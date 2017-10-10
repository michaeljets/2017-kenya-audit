import kenya
def test_kenya_default():
	
    audit_file = 'audit.csv'
    sites_file = 'sites.csv'
    n_trials = 1000
    expected_counts = [400, 300]
    counts = kenya.read_files(sites_file, audit_file)
    assert(counts==expected_counts)
    expected_audit_win_counts = {'Gnu': 613, None: 387}  
    audit_win_counts = kenya.audit(n_trials)
    assert(expected_audit_win_counts == audit_win_counts)
    
def test_random_pick():
    
    audit_file = 'audit.csv'
    sites_file = 'sites.csv'
    n_trials = 1000
    expected_counts = [400, 300]
    counts = kenya.read_files(sites_file, audit_file)
    expected_pick ='site003'
    pick =  kenya.random_pick(['site002','site005','site003'])
    assert(expected_pick == pick)

def test_local_pick():
    audit_file = 'audit.csv'
    sites_file = 'sites.csv'
    n_trials = 1000
    expected_counts = [400, 300]
    counts = kenya.read_files(sites_file, audit_file)
    expected_pick ='site005'
    pick =  kenya.local_pick('county002',['site002','site005','site003'])
    assert(expected_pick == pick)

def test_pick():
    audit_file = 'audit.csv'
    sites_file = 'sites.csv'
    n_trials = 1000
    expected_counts = [400, 300]
    counts = kenya.read_files(sites_file, audit_file)

    
    
    #global case
    expected_pick ='site005'
    pick =  kenya.pick('county002',['site002','site005','site003'],  t= 1,g =1)
    assert(expected_pick == pick)

    
    #global case g = 0 but  t > candidates
    expected_pick ='site002'
    pick =  kenya.pick('county002',['site002','site005','site003'], t= 4, g =0)
    
    assert(expected_pick == pick)

    #local_case g = 0 and t < candidates
    expected_pick ='site005'
    pick =  kenya.pick('county002',['site002','site005','site003'], t= 1, g =0)
    assert(expected_pick == pick)