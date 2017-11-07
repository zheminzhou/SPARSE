========================
Representatives database
========================

In order to do read-level taxonomic binning, representative databases need to be compiled. Four default databases were designed cover most of the genetic diversities in metagenomic samples. 

**ANI 98% database for bacteria and archaea**

.. code-block::

    python 03_db_MapDB.py dbname=refseq update=representative | python 03_db_MapDB.py dbname=refseq update=representative seqlist=stdin

**ANI 99% database for bacteria and archaea (always use together with representative database)**

.. code-block::

    python 03_db_MapDB.py dbname=refseq update=subpopulation | python 03_db_MapDB.py dbname=refseq update=subpopulation seqlist=stdin

**ANI 99% virus database**

.. code-block::

    python 03_db_MapDB.py dbname=refseq update=Virus | python 03_db_MapDB.py dbname=refseq update=Virus seqlist=stdin

**ANI 99% eukaryota database (genome size <= 200MB)**

.. code-block::

   python 03_db_MapDB.py dbname=refseq update=Eukaryota | python 03_db_MapDB.py dbname=refseq update=Eukaryota seqlist=stdin

**Custom databases** 

In order to index a differet set of references into a representative database, see [here](custom.md)
