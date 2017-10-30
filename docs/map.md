### Map metagenomic reads onto representative databases:
`python 11_query_reads.py dbname=</path/to/SPARSE/database> MapDB=<comma delimited MapDB's> r1=<read_1> r2=<read_2> workspace=<workspace_name>`

* Example (single end):

`python 11_query_reads.py dbname=refseq MapDB=representative,subpopulation,Virus r1=read1.fq.gz workspace=read1`

The outputs consist of two files, with detailed information in the ["output" section](output.md).

### Extract reference specific reads:
You first need to find out the indices of the interesting references in the [output files](output.md), and use the indexes to extract related reads. 

`python 21_get_specific_reads.py dbname=refseq workspace=read1 ref_id=<comma delimited indices>`

* For example, we want reference id 16 in [the example](output.md), which is a Vibrio cholerae genome. 

`python 21_get_specific_reads.py dbname=refseq workspace=read1 ref_id=16`
