The refseq database from NCBI stores >100,000 complete genomes and drafts that compass all tree of life. 
### Empty framework
We firstly construct an empty database folder and assigns default control parameters for the database
`python 01_db_create.py dbname=refseq`

### Index up-to-date refseq database
A second command allows SPARSE to download all genomes in refseq on-fly and construct the database. The efficiency of index process depends on both the downloading speed and the number of assigned CPUs. When assigning 20 CPUs, you can expect the whole process finishes in about one day. 
`python 02_db_index.py dbname=refseq update=True`

### Update your database
Using the same command in an established database can incrementally add in new genomes that are in refseq but not present in your database. 
`python 02_db_index.py dbname=refseq update=True`

Be aware that the newly added genomes will not be included in the representative databases. [Run another command to update your representative databases](docs/representative.md).

### Custom your database
You can also create a custom database, or add in custom genomes to an old database. Find it in the [next section])(docs/insert.md).