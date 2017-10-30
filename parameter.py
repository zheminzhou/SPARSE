import os

installation_param = dict(
    BIN = '/usr/local/bin/',
    
    
    mash = '{BIN}mash',
    bowtie2 = '{BIN}bowtie2',
    bowtie2_build = '{BIN}bowtie2-build',
    samtools = '{BIN}samtools',
    malt_run = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-run',
    malt_build = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-build',
)

run_param = dict(
    mismatch = 0.05,
    n_thread = 20,
    minFreq = 0.0001, # minFreq = 0.000001 reserve for ancient DNA samples
    minNum = 10,       # minNum = 5 reserve for ancient DNA samples
    HGT_prior = [[0.05, 0.99, 0.1], [0.02, 0.99, 0.2], [0.01, 0.99, 0.5]],
    UCE_prior = [487, 2000],
)

database_param = dict(
    mash_db = '{dbname}/mash_db',
    bowtie_db = '{dbname}/bowtie_db',
    placer_db = '{dbname}/placer_db',
    taxonomy_db = '{dbname}/taxonomy',

    msh_param = '-k 23 -s 4000 -S 42',

    barcode_dist =    [0.1,   0.05,0.02,0.01,   0.005,0.002,0.001,   0.0005],
    barcode_tag =     ['u',   's' ,'r' ,'p' ,   'n'  ,'m'  ,'e'  ,   'c'    ,'a'],
    default_mash = [dict(
                        name='virus.msh',      ref=120000,       max=250000,       min=0,
                    ), 
                    dict(
                        name='default.msh',    ref=6000000,      max=10000000,     min=50000,
                    ), 
                    dict(
                        name='large.msh',      ref=12000000,     max=25000000,     min=2500000,
                    ), 
                    dict(
                        name='eukaryota.msh',  ref=120000000000, max=250000000000, min=5000000,
                    )
    ],
    default_bowtie = [dict(
                        MapDB='Virus',         tag='p==a',      min=0, max=100000000, superkingdom='Viroids,Viruses',
                    ), 
                    dict(
                        MapDB='representative',tag='r==a',      min=0, max=100000000, superkingdom='Bacteria,Archaea',
                    ), 
                    dict(
                        MapDB='subpopulation', tag='p!=r;p==a', min=0, max=100000000, superkingdom='Bacteria,Archaea',
                    ), 
                    dict(
                        MapDB='Eukaryota',     tag='p==a',      min=0, max=200000000, superkingdom='Eukaryota,-',
                    )
    ],
    representative_level = 2,
)

program_param = dict(
    SPARSE = os.path.realpath(__file__).replace('\\', '/').rsplit('/', 1)[0],
    ipopt = '{SPARSE}/EM/solve-model',    
    db_columns = ['index', 'deleted', 'barcode', 'sha256', 'size'],
    metadata_columns = ['assembly_accession', 'version', 'refseq_category', 'assembly_level', 'taxid', 'organism_name', 'file_path', 'url_path'],
    taxa_columns = ['subspecies', 'species', 'genus','family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],
)

default_param = program_param
default_param.update(installation_param)
default_param.update(run_param)
default_param.update(database_param)
