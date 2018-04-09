===================================================================
MASH based taxonomic assignment for genomic assemblies or read sets
===================================================================

SPARSE allows ultra-efficient taxonomic assignment with genomic assemblies or read sets, by using MASH to approximate average nucleotide identities (ANI). 

*genomic assembly (fasta format):*

.. code-block:: bash

    python SPARSE.py sample --dbname refseq --query <assembly file>

*Read set (in fastq format, either gzipped or not) :*

.. code-block:: bash

    python SPARSE.py sample --dbname refseq --query <read file> --read
