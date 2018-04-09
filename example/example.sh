#!/bin/bash
echo ':::: Creating an empty database with a name "toyset"'
    python ../SPARSE.py create --dbname toyset

echo ':::: Filling database "toyset" with 22 Salmonella complete genomes'
    python ../SPARSE.py index --dbname toyset --seqlist Salmonella_toyset.txt

echo ':::: Building a mapping database named "Salmonella" in "toyset"'
    python ../SPARSE.py query --dbname toyset --tag m==a | python ../SPARSE.py MapDB --dbname toyset --MapDB Salmonella --seqlist stdin


echo ':::: Aligning ancient Salmonella reads onto the mapping database'
    python ../SPARSE.py predict --dbname toyset --r1 Ragna.sample.fq.gz --workspace Ragna_toy --MapDB Salmonella


echo ':::: Obtain a standard report'
    python ../SPARSE.py report Ragna_toy


echo ':::: Obtain a SPARSE-style report'
    cat Ragna_toy/profile.txt


echo ':::: Retrieve all reads that are specific to the reference'
    python ../SPARSE.py SSR --workspace Ragna_toy --ref_id 10

