# km_kenya.py
# Michael Jetsupphasuk
# October 09, 2017
# Prototype code for Frequentist RLA audits of Kenyan election

"""
Context: 

Each precinct returns a 'Form 34A' that has tally of votes.
There is also a secondary means (e.g. from video) for getting a tally.
We'll do a variant of the Kaplan-Markov comparison audit, and see if reported outcome 
is confirmed.

This routine takes as input:

   (1) A CSV file 'sites.csv' that gives, for each site (polling-site):
        (a) its 'Polling Site ID'
        (b) its 'County ID'
        (c) the number of registered voters for that pollsite 'Number voters'
        (d) the number voting for each candidate according to form 34A.

   (2) A CSV file 'audit.csv' that gives, for the sampled sites:
        (a) its 'Polling Site ID'
        (b) the number voting for each candidate in the photo/video.
       (Here "sampled" means obtaining the counts from the photo/video evidence.)

   See the associated sample data files.

This module produces a decision of whether to accept the reported winner or to run a full-hand count.

"""

import collections
import csv
import hashlib
import numpy as np
import sys
import argparse
import functools # JETS addition
import opt


##############################################################################
## global variables

candidates = []                 # list of candidates

site_ids = []                   # list of site ids

county_ids = []                 # list of county ids
number_voters = []              # list of number of voters
counts_34A = []                 # list of counts from 34As

county_id_s = {}                # county id by site
number_voters_s = {}            # number of voters by site
counts_34A_s = {}               # form 34A counts by site

counts_audit_s = {}             # audit counts by site


##############################################################################
## random integers

audit_seed = 0                  # initial value should be set from 24 decimal dice
                                # rolled at public ceremony AFTER sites.csv
                                # is created and published.

def random_int(n):
    """ Return the next random integer modulo n.
    """

    global audit_seed
    assert n > 0
    h = hashlib.sha256()
    h.update(str(audit_seed).encode())
    audit_seed = int(h.hexdigest(), 16)
    return audit_seed % n

##############################################################################

def clean(id):
    return id.strip()


def read_sites_csv(sites_file):
    """ Read the sites.csv file.
    
        This file should need to be only created once.
    """

    global candidates, site_ids, county_ids, number_voters, counts_34A
    global county_id_s, number_voters_s, counts_34A_s
    
    with open(sites_file) as file:
        reader = csv.reader(file)
        rows = [row for row in reader]
        fieldnames = rows[0]
        rows = rows[1:]

        fieldnames = [clean(fn) for fn in fieldnames]
        assert fieldnames[0] == 'Polling Site ID'
        assert fieldnames[1] == 'County ID'
        assert fieldnames[2] == 'Number voters'
        candidates = fieldnames[3:]
        assert len(candidates) > 0
        
        site_ids = [clean(row[0]) for row in rows]
        if len(site_ids) != len(set(site_ids)):
            print("Site ids are not unique!")
            c = collections.Counter(site_ids)
            for item in c:
                if c[item]>1:
                    print("ERROR: repeated site:", item, "count:", c[item])
            sys.exit()

        # site attributes as lists, in order given in file
        county_ids = [clean(row[1]) for row in rows]
        number_voters = [int(clean(row[2])) for row in rows]
        counts_34A = [[int(x) for x in row[3:]] for row in rows]

        # site attributes as dicts, mapping from site ids to values
        for i in range(len(rows)):
            county_id_s[rows[i][0]] = county_ids[i]
        for i in range(len(rows)):
            number_voters_s[rows[i][0]] = number_voters[i]
        for i in range(len(rows)):
            counts_34A_s[rows[i][0]] = counts_34A[i]

    return 


def read_audit_csv(audit_file):
    """ Read audit.csv file

        This file is appended to every time a site is audited.

        JETS NOTE: I don't understand the note above...the file is not being changed
    """

    global candidates, site_ids, county_ids, number_voters, counts_34A
    global county_id_s, number_voters_s, counts_34A_s
    global candidates_audit, site_ids_audit, counts_audit
    global counts_audit_s

    with open(audit_file) as file:
        reader = csv.reader(file)
        rows = [row for row in reader]
        fieldnames = rows[0]
        rows = rows[1:]

        fieldnames = [clean(fn) for fn in fieldnames]
        assert fieldnames[0] == 'Polling Site ID'
        candidates_audit = fieldnames[1:]
        assert len(candidates_audit) > 0
        
        site_ids_audit = [clean(row[0]) for row in rows]
        counts_audit = [[int(x) for x in row[1:]] for row in rows]
        for i in range(len(rows)):
            counts_audit_s[rows[i][0]] = counts_audit[i]

        return 


def read_files(sites_file, audit_file):
    
    global candidates, site_ids, county_ids, number_voters, counts_34A
    global county_id_s, number_voters_s, counts_34A_s
    global candidates_audit, site_ids_audit, counts_audit

    read_sites_csv(sites_file)
    print("Form 34A data:")
    print("    Candidates (34A):\n        {}".format(candidates))
    print("    Site ids:\n        {}".format(site_ids))
    print("    County ids by site:\n        {}".format(county_ids))
    print("    Number voters per site:\n        {}".format(number_voters))
    print("    Votes per candidate by site:")
    for i in range(len(counts_34A)):
        counts = counts_34A[i]
        print("       ", site_ids[i], counts)

    read_audit_csv(audit_file)
    print("Audit (video) data:")
    print("    Candidates (audit):\n        {}".format(candidates_audit))
    assert candidates == candidates_audit
    print("    Votes per candidate by audited site:")
    for i in range(len(counts_audit)):
        counts = counts_audit[i]
        print("       ", site_ids_audit[i], counts)

    return counts
        

def compute_sites_in_sample_order():
    """ Based on random number seed.

        JETS NOTE: is this necessary for Kaplan-Markov?
    """

    global site_ids
    global sites_in_sample_order

    sites_in_sample_order = []
    sites = site_ids.copy()
    while len(sites)>0:
        site = random_pick(sites)
        sites_in_sample_order.append(site)
        sites.remove(site)


def random_pick(sites):
    """ Pick and return a random site from the list 'sites'.
        The probability of picking a site is proportional to its size.

        JETS NOTE: PPEB condition - u/U
    """

    global number_voters_s

    total_size = sum([number_voters_s[site] for site in sites])
    picked_index = random_int(total_size)
    size_so_far = 0
    for i, site in enumerate(sites):
        size_so_far += number_voters_s[site]
        if picked_index < size_so_far:
            return site

        
def kaplan_markov(batches_to_sample, p_val=0.05, pad=1.03):
    """ Applies Kaplan-Markov inequality to confirm winner with p_val probability;
        otherwise, prints decision to do full hand-count. See Stark (2009).

        batches_to_sample is a list of lists of sites to sample
    
    """

    reported_margin = sum([winner-loser for winner,loser in counts_34A])
    diluted_margin = reported_margin / sum(number_voters)

    counts_34A_s_set = set(counts_34A_s)
    counts_audit_s_set = set(counts_audit_s)

    p_km = []
    for batches in batches_to_sample:
        errors = []
        upper_bounds = []
        for site in counts_34A_s_set.intersection(counts_audit_s_set, set(batches)):
            errors.append(((counts_34A_s[site][0] - counts_audit_s[site][0]) - (counts_34A_s[site][1] - counts_audit_s[site][1])) / reported_margin)
            upper_bounds.append(((counts_34A_s[site][0] - counts_34A_s[site][1]) + number_voters_s[site]) / reported_margin)
            # upper_bounds.append((2 * number_voters_s[site]) / reported_margin)

        error_normalized = sum(errors) # testing if this is >= 1
        error_max = sum(upper_bounds)

        upper_bounds = [upper_bound*pad for upper_bound in upper_bounds]
        taints = [error/upper_bound for error,upper_bound in zip(errors,upper_bounds)]

        p = []
        for batch in batches:
            p.append((1 - (1/error_max)) / (1-taints[0]))
            if len(taints) > 1:
                taints = taints[1:]
        p_km.append(functools.reduce(lambda x,y: x*y, p))

    final_p = min(p_km)

    if final_p <= p_val:
        print("Do not reject election results with probability >= {}.".format(1-p_val))
        # print("Examined {} batches by hand".format())
    else:
        print("Cannot confirm election results with probability >= {}. Full hand count required".format(1-p_val))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs simulation for Kenyan election audit.')
    parser.add_argument('--n_trials', type=int,
                        help='Specifies Number of trials to run, default 1000', required=False)
    parser.add_argument('--sites_file', type=str,
                        help="""(1) A CSV file  default 'sites.csv' that gives, for each site (polling-site):
        (a) its 'Polling Site ID'
        (b) its 'County ID'
        (c) the number of registered voters for that pollsite 'Number voters'
        (d) the number voting for each candidate according to form 34A.""", required=False)
    parser.add_argument('--audit_file', type=str,
                        help="""2 A CSV file default 'audit.csv' that gives, for the sampled sites:
        (a) its 'Polling Site ID'
        (b) the number voting for each candidate in the photo/video.
       (Here "sampled" means obtaining the counts from the photo/video evidence.)""", required=False)

    parser.set_defaults(sites_file='sites.csv')
    parser.set_defaults(audit_file= 'audit.csv')
    parser.set_defaults(n_trials=1000)
    args = parser.parse_args()
    #read_files(args.sites_file, args.audit_file)    
    #audit(args.n_trials)
    counts = read_files(args.sites_file, args.audit_file)
    
