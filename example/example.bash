echo ':::: Creating an empty database with a name "toyset"'
python 01_db_create.py dbname=toyset

echo ':::: Filling database "toyset" with 22 Salmonella complete genomes'
python 02_db_index.py dbname=toyset seqlist=example/Salmonella_toyset.txt

echo ':::: Building a mapping database named "Salmonella" in "toyset"'
python 10_query_metadata.py dbname=toyset tag=m==a | python 03_db_MapDB.py dbname=toyset MapDB=Salmonella seqlist=stdin

echo ':::: Aligning ancient Salmonella reads onto the mapping database'
python 11_query_reads.py dbname=toyset r1=example/Ragna.sample.fq.gz workspace=Ragna_toy MapDB=Salmonella

echo ':::: Obtain a standard report'
python 31_traditional_report.py Ragna_toy

echo ':::: Obtain a SPARSE-style report'
cat Ragna_toy/profile.txt

echo ':::: Retrieve all reads that are specific to the reference'
python 21_get_specific_reads.py workspace=Ragna_toy ref_id=10

