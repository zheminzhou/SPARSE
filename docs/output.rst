========================================
Outputs
========================================

Output for 'sparse query'
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

  <SPARSE tag>	<% in total reads>	<% in matched reads>	<taxonomic labels> (<reference IDs>)


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


The SPARSE tags are internal hierarchical clustering levels used in SPARSE database. It consists two components. The prefix describes the ANI level of the cluster and the following number is the designation of the cluster. 

For example, 's154' is a cluster '154' in 's' level (ANI 95%)
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
