==================
Parameter settings
==================
All the default parameters are stored in [parameter.py](../parameter.py). Some parameters here need to be specified during installation, while others can be specified for each database or for each SPARSE run

#### parameters that need to be specified during installation
You need only point `BIN` to a folder that contains all the executables of the [dependencies](installation.md), e. g.

* `BIN = '/usr/local/bin/'`

Alternatively, if you have all executables in the system environmental parameter `$PATH`, use 

* `BIN = ''`

You can also specify a pointer for each executable file :

* `mash = '{BIN}mash',`

* `bowtie2 = '{BIN}bowtie2',`

* `bowtie2_build = '{BIN}bowtie2-build',`

* `samtools = '{BIN}samtools',`

* `malt_run = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-run',`

* `malt_build = 'xvfb-run --auto-servernum --server-num=1 {BIN}malt-build',`


#### The following parameters that can be specified on-fly. You can also specify there default values for each database in: `/path/to/sparse/database/dbsetting.cfg`


* `mismatch = 0.05`                                                       # mismatch parameter is used in the probalistic model. Given a higher value will report less bins


* `n_thread = 20`                                                          # number of threads for SPARSE. Higher value can accelerate the program

* `minFreq = 0.0001`                                                       # Minimum frequencies of a strain to be reported. Use minFreq = 0.000001 for ancient DNA samples

* `minNum = 10`                                                            # Minimum number of specific reads to report a strain. Use * minNum = 5 or less for ancient DNA samples

* `HGT_prior = [[0.05, 0.99, 0.1], [0.02, 0.99, 0.2], [0.01, 0.99, 0.5]]`  # parameters to identify core genomic regions. Suggest to use default values

* `UCE_prior = [487, 2000]`                                                # parameters to identify ultra-conserved elements. Suggest to use default values

#### parameters to construct SPARSE databases, only for advanced uses:
`msh_param = '-k 23 -s 4000 -S 42'`                                        # change the parameter for the MASH program. reduce k and s accelerate the database indexing while bring in slightly more incorrect clusterings

`# following three parameters are pointers to corresponding sub-folders. Change them if you want the actual data in a different folder than the database`

* `mash_db = '{dbname}/mash_db'`

* `bowtie_db = '{dbname}/bowtie_db'`

* `placer_db = '{dbname}/placer_db'`

* `taxonomy_db = '{dbname}/taxonomy'`

`# parameters for hierarchical clustering levels`

* `barcode_dist =    [0.1,   0.05,0.02,0.01,   0.005,0.002,0.001,   0.0005]`

* `barcode_tag =     ['u',   's' ,'r' ,'p' ,   'n'  ,'m'  ,'e'  ,   'c'    ,'a']`

* `representative_level = 2`

#### These parameters are for experts, and have not been tested for varied values

`SPARSE = sparse_folder`

`ipopt = '{SPARSE}/EM/solve-model'`

`db_columns = ['index', 'deleted', 'barcode', 'sha256', 'size']`

`metadata_columns = ['assembly_accession', 'version', 'refseq_category', 'assembly_level', 'taxid', 'organism_name', 'file_path', 'url_path']`

`taxa_columns = ['subspecies', 'species', 'genus','family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],`
