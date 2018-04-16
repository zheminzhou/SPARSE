#!/bin/bash
echo ':::: Creating an empty database with a name "toyset"'
    sparse init --dbname toyset

echo ':::: Filling database "toyset" with 22 Salmonella complete genomes'
    sparse index --dbname toyset --seqlist Salmonella_toyset.txt

echo ':::: Building a mapping database named "Salmonella" in "toyset"'
    sparse query --dbname toyset --tag m==a | sparse mapDB --dbname toyset --mapDB Salmonella --seqlist stdin


echo ':::: Aligning ancient Salmonella reads onto the mapping database'
    sparse predict --dbname toyset --r1 Ragna.sample.fq.gz --workspace Ragna_toy --mapDB Salmonella


echo ':::: Obtain a standard report'
    sparse report Ragna_toy


echo ':::: Obtain a SPARSE-style report'
    cat Ragna_toy/profile.txt


echo ':::: Retrieve all reads that are specific to the reference'
    sparse extract --workspace Ragna_toy --ref_id 10

