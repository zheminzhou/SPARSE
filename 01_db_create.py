import sys, os, msgpack

default_param = dict(
    HOME = os.path.expanduser('~'), 

    barcode_dist =    [0.1,   0.05,0.02,0.01,   0.005,0.002,0.001,   0.0005], 
    barcode_tag =     ['u',   's','r','p',      'n','m','e',         'c','a'],
    db_columns = ['index', 'barcode', 'sha256', 'size'], 
    metadata_columns = ['assembly_accession', 'taxid', 'organism_name', 'file_path', 'url_path'],
    taxa_columns = ['subspecies', 'species', 'genus','family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],

    mash_db = '{dbname}/mash_db',
    bowtie_db = '{dbname}/bowtie_db',
    taxonomy_db = '{HOME}/software/FBDA/taxonomy',
    default_MapDB = 'default', 
    
    mash = '/usr/local/bin/mash', 
    msh_param = '-k 23 -s 4000', 

    bowtie2 = '{HOME}/bin/bowtie2', 
    bowtie2_build = '{HOME}/bin/bowtie2-build',
    samtools = '{HOME}/bin/samtools',
    ipopt = '{HOME}/software/FBDA/EM/sigma-solve-model',
    
    mismatch = 0.05,
    n_thread = 20,
    representative_level = 3,
    default_ref_size = 150000000,
)


if __name__ =='__main__' :
    param = dict([[ k.strip() for k in arg.split('=', 1)] for arg in sys.argv[1:]])
    assert 'dbname' in param, 'need "dbname".'
    
    for p in default_param :
        if p not in param :
            param[p] = default_param[p]
        elif p == 'barcode_dist':
            param[p] = [float(pp) for pp in p.split(',')]
        elif isinstance(default_param[p], int) :
            param[p] = int(param[p])

    for k in param :
        if k != 'barcode_pattern' and isinstance(param[k], basestring) and '{' in param[k] :
            param[k] = param[k].format(**param)


    if not os.path.isdir(param['dbname']) :
        os.mkdir(param['dbname'])
    if not os.path.isdir(param['mash_db']) :
        os.mkdir(param['mash_db'])
    if not os.path.isdir(param['bowtie_db']) :
        os.mkdir(param['bowtie_db'])
    import json
    json.dump(param, open(os.path.join(param['dbname'], 'dbsetting.cfg'), 'w'))