===============
RefSeq database
===============

The refseq database from NCBI stores >100,000 complete genomes and drafts that compass all tree of life. 
We firstly construct an empty database folder and assigns default control parameters for the database.

.. code-block:: bash

    python SPARSE.py init --dbname refseq

Index refseq database or update an exising database
---------------------------------------------------
A second command allows SPARSE to download all genomes in refseq on-fly and construct the database. The efficiency of the indexing process depends on both the downloading speed and the number of assigned CPUs. When assigning 20 CPUs, you can expect the whole process to finish in about one day. 

.. code-block:: bash

    python SPARSE.py index --dbname refseq --update

Be aware that the newly added genomes are not ready for metagenomic reads. You need to run another command to update your representative databases.

We also release a pre-compiled database named "refseq_20180519", on the basis of NCBI RefSeq at 2018.05.19, at 
http://enterobase.warwick.ac.uk/sparse/

This database contains the MASH indexed master database and four default mapping databases:

.. code-block:: bash

    representative
    subpopulation
    Virus
    Eukaryota

As well as reference genomes for three important animal hosts:

.. code-block:: bash

    Human
    Swine
    Bovine


To use the database, just download and untar the package:

.. code-block:: bash

    curl -O refseq_20180531.tar.gz http://enterobase.warwick.ac.uk/sparse/refseq_20180531.tar.gz
    tar -vxzf refseq_20180531.tar.gz


Custom databases
----------------

You can also create a custom database, or add in custom genomes to an old database. 
