=====================
Toy example of SPARSE
=====================

There is a bash script in the examples folder that runs through an example of SPARSE. This example takes reads from an ancient DNA sample of 800 year old *Salmonella enterica* `from this publication <http://www.biorxiv.org/content/early/2017/02/03/105759>`_.  

The contents of the shell script download and create a Bowtie database of some *Salmonella* genomes (01-03) and then maps some of the ancient reads and displays a report.

.. code-block:: bash

    python ../SPARSE.py create --dbname toyset
    python ../SPARSE.py index --dbname toyset --seqlist Salmonella_toyset.txt
    python ../SPARSE.py query --dbname toyset --tag m==a | python ../SPARSE.py MapDB --dbname toyset --MapDB Salmonella --seqlist stdin
    python ../SPARSE.py predict --dbname toyset --r1 Ragna.sample.fq.gz --workspace Ragna_toy --MapDB Salmonella
    python ../SPARSE.py report Ragna_toy
    cat Ragna_toy/profile.txt
    python ../SPARSE.py SSR --workspace Ragna_toy --ref_id 10
