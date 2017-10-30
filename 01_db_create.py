import sys, os, json,subprocess
from parameter import default_param

if __name__ =='__main__' :
    param = dict([[ k.strip() for k in arg.split('=', 1)] for arg in sys.argv[1:]])
    assert 'dbname' in param, 'need "dbname".'
    
    for p in default_param :
        if p not in param :
            param[p] = default_param[p]
        elif p == 'barcode_dist':
            param[p] = [float(pp) for pp in param[p].split(',')]
        elif isinstance(default_param[p], list) :
            param[p] = param[p].split(',')
        elif isinstance(default_param[p], int) :
            param[p] = int(param[p])

    if not os.path.isdir(param['dbname']) :
        os.mkdir(param['dbname'])

    mash_db = param['mash_db'].format(**param)
    bowtie_db = param['bowtie_db'].format(**param)
    placer_db = param['placer_db'].format(**param)
    taxonomy = param['taxonomy_db'].format(**param)
    
    for folder in (mash_db, bowtie_db, placer_db, taxonomy) :
        if not os.path.isdir(folder) :
            os.makedirs(folder)
    
    json.dump(param, open(os.path.join(param['dbname'], 'dbsetting.cfg'), 'w'), indent=2, sort_keys=True)
