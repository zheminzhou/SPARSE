===================================================
Map metagenomic reads onto representative databases
===================================================

Template of how to map metagenomic reads onto a representative database::
    
    sparse predict --dbname </path/to/SPARSE/database> --mapDB <comma delimited MapDB's> --r1 <read_1> --r2 <read_2> --workspace <workspace_name>

Example (single end):

.. code-block:: bash

    sparse predict --dbname refseq --mapDB representative,subpopulation,Virus --r1 read1.fq.gz --workspace read1

The outputs consist of two files, with detailed information in the ["output" section](output.md).

Extract reference specific reads
--------------------------------

You first need to find out the indices of the interesting references in the [output files](output.md), and use the indexes to extract related reads. 

.. code-block:: bash

    sparse extract --dbname refseq --workspace read1 --ref_id <comma delimited indices>

For example, we extract all reads specific to reference id 16, which is a Vibrio cholerae genome. 

.. code-block:: bash

    sparse extract --dbname refseq --workspace read1 --ref_id 16
