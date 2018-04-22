# Strain Prediction and Analysis using Representative SEquences (SPARSE)

SPARSE indexes >100,000 reference genomes in public databases in to hierarchical clusters and uses it to predict origins of metagenomic reads. 


[![Build Status](https://travis-ci.org/zheminzhou/SPARSE.svg?branch=master)](https://travis-ci.org/zheminzhou/SPARSE)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Docs Status](https://readthedocs.org/projects/sparse/badge/)](http://sparse.readthedocs.io/en/latest/)

## Installation 

SPARSE runs on Unix and requires Python >= version 2.7

**System modules (Ubuntu 16.04) :**

* pip
* gfortran
* llvm
* libncurses5-dev
* cmake
* xvfb-run (for malt, optional)

**3rd-party software:**
* samtools (>=1.2)
* mash (>=1.1.1)
* bowtie2 (>=2.3.2)
* malt (>=0.4.0) (optional)

See [requirements.txt](requirements.txt) for python module dependencies. 

### Installation via PIP

    pip install meta-sparse

### Installation from source codes (Ubuntu) 
     
    sudo apt-get update
    sudo apt-get install gfortran llvm libncurses5-dev cmake python-pip samtools bowtie2
    git clone https://github.com/zheminzhou/SPARSE
    cd SPARSE/EM && make
    pip install -r requirements.txt 
    

### Updating SPARSE
To update SPARSE, move to installation directory and pull the latest version:  

    cd SPARSE
    git pull
    
    
## Quick Start
See http://sparse.readthedocs.io/en/latest/ for full documentation.

1. **Download reference database**

We provide a pre-compiled database based on RefSeq (dated 14.10.2017) to download at http://enterobase.warwick.ac.uk/sparse/
   Please download the complete folder refseq_20171014/ and do not change its internal folder structure. The database can be unpacked by running:
   ```
   cd refseq_20171014 && sh untar.bash
   ```
   This pre-compiled database contains four default mapping databases, which can be specified in the next step: representative, subpopulation, Virus, Eukaryota.
   
   To update the database or build a costum database, please refer to the full documentation.
   
2. **Predict read origins**

This following command will map and evaluate all reads in both fastq-files against the specified mapping databases. 
```
python SPARSE.py predict --dbname refseq_20171014 --MapDB representative,subpopulation,Virus,Eukaryota --r1 read1.fq.gz --r2 read2.fq.gz --workspace <workspace_name>
```
For single-end reads, only --r1 needs to be specified. All output files are stored in the respective workspace.

3. **Create a report**
```
python SPARSE.py report <workspace_name>
```
The report will be stored in <workspace_name>/profile.txt

4. **Extract reference specific reads**

The following command extracts all reads specific to the provided reference ids, which can be found in the output of step 2.
```
python SPARSE.py extract --dbname refseq_20171014 --workspace <workspace_name> --ref_id <comma delimited indices>
```



## Citation
SPARSE has not been formally published yet. If you use SPARSE please cite the preprint https://www.biorxiv.org/content/early/2017/11/07/215707

Zhemin Zhou, Nina Luhmann, Nabil-Fareed Alikhan, Christopher Quince, Mark Achtman, 'Accurate Reconstruction of Microbial Strains Using Representative Reference Genomes' bioRxiv 215707; doi: https://doi.org/10.1101/215707
