In order to do read-level taxonomic binning, representative databases need to be compiled. Four default databases were designed cover most of the genetic diversities in metagenomic samples. 
### ANI 98% database for bacteria and archaea
`python 03_db_MapDB.py dbname=refseq update=representative | python 03_db_MapDB.py dbname=refseq update=representative seqlist=stdin`

### ANI 99% database for bacteria and archaea (always use together with representative database)
`python 03_db_MapDB.py dbname=refseq update=subpopulation | python 03_db_MapDB.py dbname=refseq update=subpopulation seqlist=stdin`

### ANI 99% virus database
`python 03_db_MapDB.py dbname=refseq update=Virus | python 03_db_MapDB.py dbname=refseq update=Virus seqlist=stdin`

### ANI 99% eukaryota database (genome size <= 200MB)
`python 03_db_MapDB.py dbname=refseq update=Eukaryota | python 03_db_MapDB.py dbname=refseq update=Eukaryota seqlist=stdin`
