Microbiome Amplification Preference Tool (MAPT)
=======
MAPT, the Microbiome Amplification Preference Tool, is designed to find regions of homology between the host and its microbiome. MAPT has great flexibility, as any FASTA sequence can be used as the target sequence, i.e the sequence that should be blocked. As you can see below, for the gPNA that sequence is Medicago sativa. MAPT users also specify the primers (shown below as the forward and reverse primer as Fprimer and Rprimer. This helps as only organisms that have homology to these primers are amplified in silico. The microbiome sequences can come from SILVA, and can be specified by taxonomy, as well as imported as FASTA sequences. K-mer range is the size that the sequences are fractionated into, such that K-mer range of 9-13 fractionates reads into sizes 9 to 13. By doing so we can identify regions that are ideal for targeting host amplification with oligonucleotides, such as PNAs, that will bind to host DNA and block amplification.  




##Installation
Python >=3.7 is required to run MAPT
Navigate to your directory of choice and clone our github repo
`git clone https://github.com/SEpapoulis/MAPT.git`
Once MAPT has been downloaded, you can easily add it to your sys path
`import sys`
`sys.path.insert(0, 'local/path/to/MAPT)`
You can see an example jupyter notebook of this module at /Moccia_PNA/Moccia_PNA.ipynb

##Documentation
Documentation can be found under docs/_build/html/index.html
