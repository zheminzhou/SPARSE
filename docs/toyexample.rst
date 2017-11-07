=====================
Toy example of SPARSE
=====================

There is a bash script in the examples folder that runs through an example of SPARSE. This example takes reads from an ancient DNA sample of 800 year old *Salmonella enterica* `from this publication <http://www.biorxiv.org/content/early/2017/02/03/105759>`_.  

The contents of the shell script download and create a Bowtie database of some *Salmonella* genomes (01-03) and then maps some of the ancient reads and displays a report.

.. code-block:: bash

    python ../01_db_create.py dbname=toyset
    python ../02_db_index.py dbname=toyset seqlist=Salmonella_toyset.txt
    python ../10_query_metadata.py dbname=toyset tag=m==a | python ../03_db_MapDB.py dbname=toyset MapDB=Salmonella seqlist=stdin
    python ../11_query_reads.py dbname=toyset r1=Ragna.sample.fq.gz workspace=Ragna_toy MapDB=Salmonella
    python ../31_traditional_report.py Ragna_toy
    cat Ragna_toy/profile.txt
    python ../21_get_specific_reads.py workspace=Ragna_toy ref_id=10
