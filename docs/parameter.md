All the default parameters are stored in [parameter.py](../parameter.py). Some parameters here need to be specified during installation, while others can be specified for each database or for each SPARSE run

#### parameters that need to be specified during installation
You need only point `BIN` to a folder that contains all executables, such as

* `BIN = '/usr/local/bin/'`

Alternatively, if you have all executable in your `PATH`, use 

* `BIN = ''`

You can also specify the pointer to each executable :

* `mash = '{BIN}mash',`

* `bowtie2 = '{BIN}bowtie2',`

* `bowtie2_build = '{BIN}bowtie2-build',`

* `samtools = '{BIN}samtools',`

* `malt_run = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-run',`

* `malt_build = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-build',`


#### The following parameters that can be specified on-fly. Their default values for ech database are stored in: `/path/to/sparse/database/dbsetting.cfg`


* `mismatch = 0.05,                                                       # mismatch parameter is used in the probalistic model. Given a higher value will report less references, otherwise more`


* `n_thread = 20,                                                          # number of threads for SPARSE. Higher value can accelerate the program`

* `minFreq = 0.0001,                                                       # Minimum frequencies of a strain to be reported. Use minFreq = 0.000001 for ancient DNA samples`

* `minNum = 10,                                                            # Minimum number of specific reads to report a strain. Use * minNum = 5 or less for ancient DNA samples`

* `HGT_prior = [[0.05, 0.99, 0.1], [0.02, 0.99, 0.2], [0.01, 0.99, 0.5]],  # parameters to identify core genomic regions. Suggest to use default values`

* `UCE_prior = [487, 2000],                                                # parameters to identify ultra-conserved elements. Suggest to use default values`

#### parameters to construct SPARSE databases, only for advanced uses:
`msh_param = '-k 23 -s 4000 -S 42',  # change the parameter for MASH program. reduce k and s accelerate the databsae index while allowing slightly more incorrect clusterings`

`# following three parameters are pointers to corresponding sub-folders. Change them if you want the actual data in a different folder than the database`
* `mash_db = '{dbname}/mash_db',`

* `bowtie_db = '{dbname}/bowtie_db',`

* `placer_db = '{dbname}/placer_db',`

* `taxonomy_db = '{dbname}/taxonomy',`

`# hierarchical clustering levels`

* `barcode_dist =    [0.1,   0.05,0.02,0.01,   0.005,0.002,0.001,   0.0005],`
* `barcode_tag =     ['u',   's' ,'r' ,'p' ,   'n'  ,'m'  ,'e'  ,   'c'    ,'a'],`

`# default databases`  
* `default_mash = [dict(`  
`						name='virus.msh',      ref=120000,       max=250000,       min=0,`
`					), dict(`
`						name='default.msh',    ref=6000000,      max=10000000,     min=50000,`
`					), dict(`
`						name='large.msh',      ref=12000000,     max=25000000,     min=2500000,`
`					), dict(`
`						name='eukaryota.msh',  ref=120000000000, max=250000000000, min=5000000,`
`				)],`
default_bowtie = [dict(
						MapDB='Virus',         tag='p==a',      min=0, max=100000000, superkingdom='Viroids,Viruses',
					), dict(
						MapDB='representative',tag='r==a',      min=0, max=100000000, superkingdom='Bacteria,Archaea',
					), dict(
						MapDB='subpopulation', tag='p!=r;p==a', min=0, max=100000000, superkingdom='Bacteria,Archaea',
					), dict(
						MapDB='Eukaryota',     tag='p==a',      min=0, max=200000000, superkingdom='Eukaryota,-',
				)],
representative_level = 2,`

#### These parameters are for experts, and have not been tested for varied values
`SPARSE = sparse_folder,
ipopt = '{SPARSE}/EM/solve-model',
db_columns = ['index', 'deleted', 'barcode', 'sha256', 'size'],
metadata_columns = ['assembly_accession', 'version', 'refseq_category', 'assembly_level', 'taxid', 'organism_name', 'file_path', 'url_path'],
taxa_columns = ['subspecies', 'species', 'genus','family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],`
