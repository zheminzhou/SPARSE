==================
Installation guide
==================

SPARSE runs on Unix and requires Python >= version 2.7

**hardware**

* SPASE runs best in multi-processes mode. Thus servers with at least 10-20 CPU cores were suggested.
* >= 300 GBytes of memory is required to handle over 10 million metagenomic reads.
* >= 500 GBytes of storage space.

**System modules (Ubuntu 16.04):**

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

See requirements.txt for python module dependencies.

Installation with a Conda environment (Ubuntu)
-------------------------------------
To install SPARSE and all its dependencies in an isolated `conda <https://conda.io/miniconda.html>`_ environment, you can use the environment file provided.

.. code-block:: bash

    git clone git@github.com:zheminzhou/SPARSE.git
    cd SPARSE/envs
    conda env create -f sparse_env.yml
    source activate sparse


Installation via PIP
--------------------

.. code-block:: bash

    pip install meta-sparse


Install from source files (Ubuntu)
----------------------------------

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install gfortran llvm libncurses5-dev cmake python-pip samtools bowtie2
    git clone https://github.com/zheminzhou/SPARSE
    cd SPARSE/EM && make
    pip install -r requirements.txt

Change the parameters if needed.


Updating SPARSE
---------------

To update SPARSE, move to the installation directory and pull the latest version:

.. code-block:: bash

    cd SPARSE
    git pull
