.. currentmodule:: MAPT

****************************************
recapitulating the pPNA
****************************************

General explanation of the procedure
MAPT, the Microbiome Amplification Preference Tool, is designed to find regions of homology between the host and its microbiome. MAPT has great flexibility, as any FASTA sequence can be used as the target sequence, i.e. the sequence that should be blocked. As you can see below, for the pPNA that sequence is A. thaliana, and is designed to recapitulate the results seen in Lundberg et al., 2013. MAPT users also specify the primers (shown below as the forward and reverse primer as Fprimer and Rprimer. This helps as only organisms that have homology to these primers are amplified in silico. The microbiome sequences come from Greengenes, as this was used to design the pPNA originally (Lundberg et al., 2013). K-mer range is the size that the sequences are fractionated into, such that K-mer range of 9-13 fractionates reads into sizes 9 to 13. By doing so we can identify regions that are ideal for targeting host amplification with oligonucleotides, such as PNAs, that will bind to host DNA and block amplification.   

Required Information
--------------------------------------
Explanation of the materials needed (target sequence,
sequence file, primers, and k-mer range)

.. ipython::

    In [1]: import MAPT

    In [2]: greengenes_fasta = "data/gg_97_otus_4feb2011.fasta"

    In [3]: Athaliana_chloro="data/A_thaliana_chloroplast.fasta"

    In [4]: R926="CCGYCAATTYMTTTRAGTTT"

    In [5]: F515 = "GTGCCAGCMGCCGCGGTAA"

    In [6]: krng = (9,13)

    In [7]: result_file='data/pPNA.csv'

Now that we have all of our materials, we can begin mapping
k-mers from our greengenes data-set to the *A. thaliana* chloroplast

Generating K-mer maps
--------------------------------------

.. ipython::

    In [7]: designer_pPNA = MAPT.PNA_Designer(result_file = 'data/pPNA.csv', target_fastafile = Athaliana_chloro, sequence_file = greengenes_fasta, primer_F=F515,primer_R = R926, kmer_range = krng)

MAPT has automatically detected that some of the sequences
in our greengenes data-set failed to find **exact** primer matches. We
can view these sequences by looking at the ‘failed_amplification’ attribute
of MAPT.


With our k-mers mapped, we are able to visualize our data! But before
we do, we will quickly map the pPNA. While PNAs have a higher affinity
to the antiparallel direction, they can also bind in the parallel direction
as well. We will allow the PNA to bind in both directions for plotting

Mapping a PNA
--------------------------------------

.. ipython::

    In [9]: pPNA = "GGCTCAACCCTGGACAG"

    In [9]: mapping = designer_pPNA.map_PNA(pPNA,antiparallel_only=False)

    In [9]: import pandas#pandas is a module for working with dataframe objects in python, highly recommend!

    In [9]: pandas.DataFrame(mapping)

Now that we have data tables, use your favorite software to make
a figure! We will use a custom R script to plot our data.

