import sys, os, json,subprocess

sparse_folder = os.path.realpath(__file__).replace('\\', '/')

default_param = dict(
    BIN = os.path.expanduser('~').replace('\\', '/') + '/bin', 
    SPARSE = sparse_folder,

    mash = '{BIN}/mash', 
    bowtie2 = '{BIN}/bowtie2', 
    bowtie2_build = '{BIN}/bowtie2-build',
    samtools = '{BIN}/samtools',
    ipopt = '{SPARSE}/EM/solve-model',
    
    barcode_dist =    [0.1,   0.05,0.02,0.01,   0.005,0.002,0.001,   0.0005], 
    barcode_tag =     ['u',   's' ,'r' ,'p' ,   'n'  ,'m'  ,'e'  ,   'c'    ,'a'],

    db_columns = ['index', 'deleted', 'barcode', 'sha256', 'size'], 
    metadata_columns = ['assembly_accession', 'version', 'refseq_category', 'assembly_level', 'taxid', 'organism_name', 'file_path', 'url_path'],
    taxa_columns = ['subspecies', 'species', 'genus','family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],

    mash_db = '{dbname}/mash_db',
    bowtie_db = '{dbname}/bowtie_db',
    placer_db = '{dbname}/placer_db',
    taxonomy_db = '{dbname}/taxonomy',
    default_mash = [dict( 
                            name='virus.msh',      ref=120000,       max=250000,       min=0,
                        ), dict(
                            name='default.msh',    ref=6000000,      max=10000000,     min=50000,
                        ), dict(
                            name='large.msh',      ref=12000000,     max=25000000,     min=2500000,
                        ), dict(
                            name='eukaryota.msh',  ref=120000000000, max=250000000000, min=5000000,
                    )],
    default_bowtie = [dict( 
                            MapDB='Virus',         tag='p==a',      min=0, max=100000000, superkingdom='Viroids,Viruses',
                        ), dict(
                            MapDB='representative',tag='r==a',      min=0, max=100000000, superkingdom='Bacteria,Archaea',
                        ), dict(
                            MapDB='subpopulation', tag='p!=r;p==a', min=0, max=100000000, superkingdom='Bacteria,Archaea',
                        ), dict(
                            MapDB='Eukaryota',     tag='p==a',      min=0, max=500000000, superkingdom='Eukaryota,-',
                    )],
    
    msh_param = '-k 23 -s 4000', 
    mismatch = 0.05,
    n_thread = 20,
    representative_level = 2,
    virus_size = 50000,
    default_ref_size = 150000000,
)

taxdump = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'

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