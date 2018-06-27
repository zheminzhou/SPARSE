========================================
Outputs
========================================

Output for 'sparse predict'
-------------------------------
The taxonomic profiling results for 'sparse query' is saved in <workspace>/profile.txt

The first three rows in 'profile.txt' summarizes the status of reads

.. code-block:: bash

  Total	<No. reads>	<No. matched reads>
  Unmatched	<% unmatched reads in total reads>	0.000
  Uncertain_match	<% unreliable matches in total reads>	<% unreliable matches in total matches>


For example

.. code-block:: bash

  Total	26726530	23783388.0
  Unmatched	11.012	0.000
  Uncertain_match	36.102	40.570



Suggests that 89% of reads found at least one match in the reference database. And 60% of all found matches are used for taxonomic predictions. 


The following lines describe the prediction for each taxonomic levels, in format:

.. code-block:: bash

  <SPARSE group>	<% in total reads>	<% in matched reads>	<taxonomic labels> (<reference IDs>)


For example

.. code-block:: bash

  ~154    2.1706  2.4392  Bacteria|-|Actinobacteria|Actinobacteria|Micrococcales|Micrococcaceae (15969,66991,66935,66915,67189,110179,40981,154,67166,67220,114405,66878,66930,82153,40861,40710,67029)
  u154    2.1701  2.4387  Bacteria|-|Actinobacteria|Actinobacteria|Micrococcales|Micrococcaceae|Rothia (15969,66991,66935,66915,67189,110179,40981,154,67166,67220,114405,66878,66930,82153,40861,40710,67029)
  s154    2.1551  2.4217  Bacteria|-|Actinobacteria|Actinobacteria|Micrococcales|Micrococcaceae|Rothia|Rothia dentocariosa (*Rothia sp. HMSC067H10/*Rothia sp. HMSC064D08/*Rothia sp. HMSC071F11/*Rothia sp. HMSC069C01) (15969,66991,66935,66915,67189,110179,40981,154,67166,67220,114405,66878,66930,82153,40861,40710,67029)
  ~613    1.4988  1.6843  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae (16778,16416,117596,16415,10931,17276,113949,60730,613)
  u613    1.4934  1.6782  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae|Veillonella (16778,16416,117596,16415,10931,17276,113949,60730,613)
  s613    1.4507  1.6302  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae|Veillonella|Veillonella parvula (*Veillonella sp. 6_1_27/*Veillonella sp. S13054-11/*Veillonella sp. 3_1_44) (16778,16416,117596,16415,10931,17276,113949,60730,613)
  r15969  0.4677  0.5256  Bacteria|-|Actinobacteria|Actinobacteria|Micrococcales|Micrococcaceae|Rothia|Rothia dentocariosa|- (15969)
  p15969  0.3907  0.4391  Bacteria|-|Actinobacteria|Actinobacteria|Micrococcales|Micrococcaceae|Rothia|Rothia dentocariosa|-|Rothia dentocariosa M567: GCF_000143585.1 (15969)
  r16416  0.1838  0.2065  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae|Veillonella|*Veillonella sp. 6_1_27|- (16416)
  p16416  0.1631  0.1833  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae|Veillonella|*Veillonella sp. 6_1_27|-|Veillonella sp. 6_1_27: GCF_000163735.1 (16416)


The SPARSE groups are internal hierarchical clustering results stored in the SPARSE database. The group label consists of two components. The prefix presents the ANI level of the cluster and the following number presents the designation of the cluster. 

For example, 's613' is a cluster '613' in 's' level (ANI 95%)
The correlation between prefix and ANI level is:

.. code-block:: bash

  ~		<90% ANI
  u		90% ANI
  s		95% ANI
  r		98% ANI
  p		99% ANI
  n		99.5% ANI
  m		99.8% ANI
  e		99.9% ANI
  c		99.95% ANI
  a		100% ANI


's' (ANI 95%) is normally treated as a 'gold standard' criterion for species definition. 

The traditional taxonomic labels for SPARSE groups are shown in format:

.. code-block:: bash

  <superkingdom>|<kingdom>|<phylum>|<class>|<order>|<family>|<genus>|<species>|<subspecies>|<reference_genome>


These taxonomic labels are summarised from the database inputed. Sometimes multiple species can be associated with one SPARSE group. For example:

.. code-block:: bash

  s613    1.4507  1.6302  Bacteria|-|Firmicutes|Negativicutes|Veillonellales|Veillonellaceae|Veillonella|Veillonella parvula (*Veillonella sp. 6_1_27/*Veillonella sp. S13054-11/*Veillonella sp. 3_1_44) (16778,16416,117596,16415,10931,17276,113949,60730,613)


Where group s613 is associated with four different species 

.. code-block:: bash

  Veillonella parvula
  *Veillonella sp. 6_1_27
  *Veillonella sp. S13054-11
  *Veillonella sp. 3_1_44
  
All these species names other than 'Veillonella parvula' starts with a prefix "*" because they are informal names. The most probable species is shown first, and followed with the other three names in a bracket. There is another bracket after the taxonomic labels. 

.. code-block:: bash

  (16778,16416,117596,16415,10931,17276,113949,60730,613)

These are the IDs of the actual reference genomes that were found in the database. They were used to extract reference specific reads using 'sparse extract' command.


Output for 'sparse report'
-------------------------------
sparse report combines multiple 'sparse predict' runs together into a tab-delimited text file, and tries to identify potential pathogens in the predictions. 

.. code-block:: bash

  #Group  #Pathogenic     ERR1659111      ERR1659110      #Species        #Taxon
  s3080   non     4.47309775569   4.84028327303   Actinomyces dentalis (*Actinomyces sp. oral taxon 414)  Bacteria|-|Actinobacteria|Actinobacteria|Actinomycetales|Actinomycetaceae|Actinomyces|Actinomyces dentalis (*Actinomyces sp. oral taxon 414)
  s1438   non     0.821962806352  3.57658189557   Desulfomicrobium orale  Bacteria|-|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfomicrobiaceae|Desulfomicrobium|Desulfomicrobium orale
  s9975   non     2.04489272864   1.85184148971   *Anaerolineaceae bacterium oral taxon 439       Bacteria|-|Chloroflexi|Anaerolineae|Anaerolineales|Anaerolineaceae|-|*Anaerolineaceae bacterium oral taxon 439
  s939    non     1.81538010098   0.712860400235  Pseudopropionibacterium propionicum     Bacteria|-|Actinobacteria|Actinobacteria|Propionibacteriales|Propionibacteriaceae|Pseudopropionibacterium|Pseudopropionibacterium propionicum
  s8820   non     1.67063037869   0.491279312566  *Ottowia sp. Marseille-P4747 (*Ottowia sp. oral taxon 894)      Bacteria|-|Proteobacteria|Betaproteobacteria|Burkholderiales|Comamonadaceae|Ottowia|*Ottowia sp. Marseille-P4747 (*Ottowia sp. oral taxon 894)
  s2215   non     1.31802856115   0.34575838713   Lautropia mirabilis     Bacteria|-|Proteobacteria|Betaproteobacteria|Burkholderiales|Burkholderiaceae|Lautropia|Lautropia mirabilis
  s2590   non     0.665641018802  0.612783437737  Actinomyces cardiffensis        Bacteria|-|Actinobacteria|Actinobacteria|Actinomycetales|Actinomycetaceae|Actinomyces|Actinomyces cardiffensis
  s2189   non     0.87220732902   0.296597041195  Corynebacterium matruchotii     Bacteria|-|Actinobacteria|Actinobacteria|Corynebacteriales|Corynebacteriaceae|Corynebacterium|Corynebacterium matruchotii
  s108979 non     0.295928369726  0.857545958706  *Actinomyces sp. oral taxon 897 Bacteria|-|Actinobacteria|Actinobacteria|Actinomycetales|Actinomycetaceae|Actinomyces|*Actinomyces sp. oral taxon 897

The first line shows the samples in the report, as well as additional annotations (starts with '#'). #Group and #Taxon are the same as in 'sparse predict' output. #Species is a simple extraction of the most probably species in #Taxon column and #Pathogenic consists of the interpretations, where

.. code-block:: bash

  non - not a pathogen
  * - commensal and normally not a pathogen
  ** - Possibly a pathogen
  *** - Pathogen
  **** - Important pathogen, possibly fatal

The numbers shows the abundances of the species in each metagenomic read set. It is normally shown in percentages, unless parameter '--absolute' is applied, which changes the numbers to be absolute read counts. 

The last row of the output is 'dark matters', which is a summary of all unknown/uncertain reads. 
