# kenya.py
# Ronald L. Rivest
# September 28, 2017
# Prototype code for Bayesian audits of Kenyan election

"""
Context: 

Each precinct returns a 'Form 34A' that has tally of votes.
There is also a secondary means (e.g. from video) for getting a tally.
We'll do a variant of Bayesian comparison audit, and see if reported outcome 
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

This module produces an estimate of the probability (in a Bayesian sense) that
each candidate would win the election, if all sites were sampled.

"""

import collections
import csv
import hashlib
import numpy as np
import sys
import argparse
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
np.random.seed(audit_seed)

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

    global candidates, site_ids, county_ids, number_voters, counts_34A, county_site_map, site_to_county_map
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

        county_site_map = {}
        site_to_county_map = {}
        for row in rows:
            site_id = clean(row[0])
            county_id = clean(row[1])
            site_to_county_map[site_id] = county_id
            if county_id in county_site_map.keys():
                county_site_map[county_id].add(site_id)
            else:
                county_site_map[county_id] = set([site_id])
  
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
        
def audit(trials):
    """ Read data files and compute probability that each candidate would
        win (in a Bayesian sense) if all sites were to be audited.
    """

    global candidates, site_ids, number_voters, counts_34A
    global candidates_audit, site_ids_audit, counts_audit
    global sites_in_sample_order

    compute_sites_in_sample_order()
    print("Sites in sample order:")
    print(sites_in_sample_order)
    # This ordering would need to be printed out and distributed for actual use...


    printing_wanted = False
    win_count = {}
    for t in range(trials):
        if printing_wanted:
            print("\nTrial:", t)
        winner = simulate()
        if printing_wanted:
            print("Winner:", str(winner))
        win_count[winner] = win_count.get(winner,0) + 1
    print("Audit results:")
    print("    For {} trials:".format(trials))
    for winner in win_count:
        print("        {} won {} times ({:0.1f}%)"
              .format(winner, win_count[winner], 100*win_count[winner]/trials))
    risk = 100*min(win_count.values())/trials
    print("    Residual risk: {:0.1f}%".format(risk))
    return win_count
def compute_sites_in_sample_order():
    """ Based on random number seed.
    """

    global site_ids
    global sites_in_sample_order

    sites_in_sample_order = []
    sites = site_ids.copy()
    while len(sites)>0:
        site = random_pick(sites)
        sites_in_sample_order.append(site)
        sites.remove(site)

def pick( county,sites ,g=1, t =4,):
    if len(sites) < t:
        return random_pick(sites)
    else:
        if np.random.uniform(0.0, 1.0) < g:
            return random_pick(sites)
        else:
            return local_pick(county,sites)

def local_pick(county, sites):
    """
        Pick a random site in a given county. 
        @params 
            county: string -- the county to pick a site from 
            sites: list the available sites
        @returns
            None if no valid choice
            site if there is a valid choice
    """
    #county_site_map: dictionary --of all polling sites by county to set of sites
    global county_site_map
    print(county_site_map)
    if county in county_site_map.keys():
        # find sites that are in a given county and in the acceptable 'sites'
        site_options =  list(county_site_map[county].intersection(set(sites)))
        return random_pick(site_options)

def random_pick(sites):
    """ Pick and return a random site from the list 'sites'.
        The probability of picking a site is proportional to its size.
    """

    global number_voters_s

    total_size = sum([number_voters_s[site] for site in sites])
    picked_index = random_int(total_size)
    size_so_far = 0
    for i, site in enumerate(sites):
        size_so_far += number_voters_s[site]
        if picked_index < size_so_far:
            return site


def simulate(printing_wanted = False):        
    """ Run one trial of a Bayesian simulation of drawing from
        the posterior distribution. Return winner for this trial.
    """

    global candidates, site_ids, number_voters, site_to_county_map,counts_34A
    global county_id_s, number_voters_s, counts_34A_s
    global candidates_audit, site_ids_audit, counts_audit
    global counts_audit_s
    global sites_in_sample_order
    global audit_tallies

    # compute/estimate audit counts for all sites
    audit_tallies = {}
    sites_considered = []
    matrix_s = {}
    for site in sites_in_sample_order:
        if site in site_ids_audit:
            # already have it
            matrix_s[site] = opt.find_A(counts_34A_s[site],
                                        counts_audit_s[site])
            audit_tallies[site] = counts_audit_s[site]
            if printing_wanted:
                print("Did:", site, "Tally:", counts_34A_s[site], audit_tallies[site])
                print(matrix_s[site])
        else:
            # pick random x-to-y matrix (Polya urn style)
            # apply it to get estimate of audit (video) counts
            picked_site = pick(site_to_county_map[site],sites_considered)
            # no prior here -- this is 'empirical bayes'
            matrix_s[site] = matrix_s[picked_site]
            audit_tallies[site] = matrix_s[site].dot(counts_34A_s[site])
            # bound size of estimated tally for each candidate by number of voters
            audit_tallies[site] = np.array([min(at, number_voters_s[site]) for at in audit_tallies[site]])
            if printing_wanted:
                print("New:", site, "Tally:", counts_34A_s[site], audit_tallies[site])
                print(matrix_s[site])
        sites_considered.append(site)

    # compute winner from audit (video) tallies
    winner = compute_winner()

    return winner
        
def compute_winner(printing_wanted=False):
    """ Compute audit winner, based on audit tallies (may return None).

        (Wikipedia:) to win in the first round, 
        a candidate must receive over 50% of the vote and 
        25% of the vote in at least 24 counties out of the 47 counties,"

        This code needs to be reviewed for tie-handling procedures.
    """

    global candidates
    global site_ids
    global county_ids
    global audit_tallies
    global county_id_s
    
    # Computer overall tally and total_votes_cast
    overall_tally = sum([audit_tallies[site] for site in site_ids])
    total_votes_cast = sum(overall_tally)
    if printing_wanted:
        print("Overall tally:", overall_tally)
        print("Total (simulated) votes cast:", total_votes_cast)

    # Compute 'winner' as the majority vote winner
    # (more than half the cast votes)
    winner = None
    for w in range(len(candidates)):
        if overall_tally[w] > total_votes_cast/2:
            winner = w
    if winner == None:
        return winner
    
    # Check that winner has more than 25% of the vote
    # in more than half of the counties
    number_votes_c = {county:0 for county in county_ids}
    for site in site_ids:
        number_votes_c[county_id_s[site]] += sum(audit_tallies[site])
    number_good_counties = 0
    for county_id in county_ids:
        if 0.25 * number_votes_c[county_id] < \
            sum([audit_tallies[site][winner] \
                 for site in site_ids \
                 if county_id_s[site] == county_id]):
            number_good_counties += 1

    if number_good_counties > 0.5 * len(county_ids):
        return candidates[winner]
    else:
        return None

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
    win_counts = audit(args.n_trials)
    print(counts)
    print(win_counts)
    
