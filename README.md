# apt-get packages
## These are packages for python executables
sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev
## These are essential for SPARSE
sudo apt-get install -y gfortran wget curl llvm libncurses5-dev cmake capnproto	

# additional dependencies
samtools-1.2 or later
mash
bowtie2-2.3.2 
## if you are looking for less similar mappings
malt


# python environment
SPARSE is developed on python 2.7.9 and can be run only later versions of python 2.7. It has not been tested in python 3. 
I suggest to use pyenv and pyenv-virtualenv to host a virtual environment for SPARSE. 

# python libraries
msgpack-python==0.4.8
numpy==1.13.1
pandas==0.20.3
psycopg2==2.6
pycapnp==0.6.0
ujson==1.35


# Installation
git clone https://github.com/zheminzhou/SPARSE.git
cd SPARSE/EM && make

# Create refseq database
## create an empty database framework
python 01_db_create.py dbname=refseq
## index up-to-date refseq database
python 02_db_index.py dbname=refseq update=True

# A SPARSE database can be used to predict taxonomic labels of a genome or a readset
## Compare to a genome (in fasta format) : 
python 12_query_sample.py dbname=refseq query=<query file>
# Compare to a read file (in fastq format, either gzipped or not) : 
python 12_query_sample.py dbname=refseq query=<read file> dtype=read

# In order to do read-level taxonomic binning, representative databases need to be compiled :
## default databases :
### ANI 98% database for bacteria and archaea
python 03_db_MapDB.py dbname=refseq update=representative | python 03_db_MapDB.py dbname=refseq update=representative seqlist=stdin
### ANI 99% database for bacteria and archaea (always use together with representative database)
python 03_db_MapDB.py dbname=refseq update=subpopulation | python 03_db_MapDB.py dbname=refseq update=subpopulation seqlist=stdin
### ANI 99% virus database
python 03_db_MapDB.py dbname=refseq update=Virus | python 03_db_MapDB.py dbname=refseq update=Virus seqlist=stdin
### ANI 99% eukaryota database (genome size <= 200MB)
python 03_db_MapDB.py dbname=refseq update=Eukaryota | python 03_db_MapDB.py dbname=refseq update=Eukaryota seqlist=stdin
## User-defined databases (use human reference genome as an example)
### You can query metadata of genomes in a SPARSE database using 10_query_metadata.py :
python 10_query_metadata.py dbname=refseq_20171014 assembly_accession=GCF_000001405.37 > human.tsv
The resulted file is like:
index	deleted barcode sha256	size	assembly_accession	version refseq_category assembly_level	taxid	organism_name	file_path	url_path	subspecies	species genus	family	order	class	phylum	kingdom superkingdom
107460	-	u107460.s107460.r107460.p107460.n107460.m107460.e107460.c107460.a107460 d236b7835a3f10e596f9ce3c1f988b9e897f2dea216fd3dcde880eb91963863e	3253848404	GCF_000001405.37	37	reference genome	Chromosome	9606	Homo sapiens	-	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_genomic.fna.gz	-	Homo sapiens	Homo	Hominidae	Primates	Mammalia	Chordata	Metazoa	Eukaryota
### Create a database for detecting human reads :
python 03_db_MapDB.py dbname=refseq MapDB=Human seqlist=human.tsv

# Align metagenomic reads onto databases
python 11_query_reads.py dbname=</path/to/SPARSE/database> MapDB=<comma delimited MapDB's> r1=<read_1> r2=<read_2> workspace=<workspace_name>
For example (single end reads):
python 11_query_reads.py dbname=refseq MapDB=representative,subpopulation,Virus r1=read1.fq.gz workspace=read1

# Extract reads that are specific to a set of refereces (reference with index 1 and 51):
python 21_get_specific_reads.py dbname=refseq workspace=read1 ref_id=1,51