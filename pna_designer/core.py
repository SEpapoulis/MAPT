import os
from . import InSilico_PCR
from .database.silva import silva_manager
from collections import defaultdict

#Microbiome Amplification Preference Tool (MAPT)
    
class k_mapper:
    '''
    Simple class for generating and mapping k-mers

    Nucleotide sequences can be sampled to find subsequences of 
    k length, otherwise refered to as k-mers. k-mers are generated from 
    a group of sequences and subsequently mapped to a target sequence. 
    Degeneracy/ambiguity in sequences are supported, however, poor quality
    sequences with too many ambiguous nucleotides will inflate mapping results.
    
    Parameters
    ----------
    target_sequence: str
        Nucleotide sequence that k-mers will be mapped to.
    sequences : list of str
        Nucleotide sequences that k-mers will be generated from.
    krange: tuple of ints, optional
        Size range of k-mers to be generated. Inclusive. Default: (9,14)

    Attributes
    ----------
    target : str
        Where k-mers are being mapped
    target_match : list of ints
        k-mer mapping result, number of times a k-mer mapped to a specific index in the target sequence
    target_match_unique : list of ints
        k-mer mapping result, number of times a unique k-mer mapped to a specific index in the target sequence

    '''

    def __init__(self,target_sequence,sequences=list(),krange=tuple([9,14])):
        #Inclusive range
        self.krng=list((range(krange[0],krange[1]+1)))

        self.target= target_sequence.upper()
        self._target_klocation = self._get_target_kmermap()
        #coordiantes for unique kmers that align
        self.target_match_unique = [0]*len(self.target)

        #coordiantes for all kmers that align, where kmers that are found
        #more often in seqlist will be weighted more
        self.target_match = [0]*len(self.target)

        #dictoinary mapping to defaultdicts of kmer counts
        #note, if _kdict[kmer] == 0, that means this kmer is unique
        #to the target sequence!
        self._kdict=defaultdict(int)

        for seq in sequences:
            self.add_sequence(seq)

    def _get_target_kmermap(self):
        kmer_location={}
        for k in self.krng:
            for i in range(k,len(self.target)+1):
                start=i-k
                stop =i 
                kmer = self.target[start:stop]
                if kmer in kmer_location:
                    kmer_location[kmer].append(range(start,stop))
                else:
                    kmer_location[kmer] = [range(start,stop)]
                #reverse complement
                kmer=InSilico_PCR.reverse_complement(kmer)
                if kmer in kmer_location:
                    kmer_location[kmer].append(range(start,stop))
                else:
                    kmer_location[kmer] = [range(start,stop)]
        return(kmer_location)

    def add_sequence(self,sequence):
        '''
        Add a sequence that will be used to generate k-mers

        Parameters
        ----------
            sequence: str
                Nucleotide sequence

        '''

        seq=sequence.upper()
        for i in range(0,len(seq)):
            for k in self.krng:
                if i+k <len(seq):
                    self._kdict[seq[i:i+k]]+=1
                    #reverse complement
                    self._kdict[InSilico_PCR.reverse_complement(seq[i:i+k])]+=1    

    def map_kmers(self):
        '''map kmers to target sequence
        '''
        for kmer in self._kdict:
            if kmer in self._target_klocation:
                locations = self._target_klocation[kmer]
                for location in locations:
                    #change range to list to indicies to add
                    loc = list(location)
                    for i in loc:
                        self.target_match_unique[i]+=1
                        self.target_match[i]+=self._kdict[kmer]

    def get_results(self):
        '''returns a dictionary with mapping results'''
        dat = {'Nucleotide':[],'index':[],'unique match':[],'absolute match':[]}
        for i in range(0,len(self.target_match)):
            dat['Nucleotide'].append(self.target[i])
            dat['index'].append(i)
            dat['unique match'].append(self.target_match_unique[i])
            dat['absolute match'].append(self.target_match[i])
        return(dat)

    def write_results(self,file_name):
        '''
        Generate mapping file

        A spreadsheet will be generated showing the unquie and
        absolute kmer mappings to the target sequence. The number
        of absolute mappings is relative to the number of sequences
        used to generate kmers, thus, absolute mapping will change with
        the size of the dataset. Columns are the nucleotide of the target,
        index of target sequence, number of unquie k-mer matches, and 
        the number of k-mer absolute matches


        Parameters
        ----------
        file_name: str
            Mapping file name
        
        '''

        with open(file_name,'w') as o:
            o.write('Nucleotide,index,unique match,absolute match')
            for i in range(0,len(self.target_match)):
                o.write('\n')
                o.write(','.join([self.target[i],str(i),
                str(self.target_match_unique[i]),str(self.target_match[i])]))
                
def map_PNA(target,PNA,krange=(9,13),antiparallel_only=False):
    '''
    Align a PNA to a given target sequence

    map_PNA will use the k_mapper class to generate k-mer alignments to target DNA. PNA oligomers
    can bind to DNA in either orientation, thus, the reverse of the PNA is also searched for alignments

    Parameters
    ----------
    target : str
        DNA sequence the PNA will be aligned to
    PNA : str
        PNA sequence
    krange : tuple, optional
        the k-mer range used for alignemnt. Defaults to 9-13 (inclusive)
    antiparallel_only : bool, optional
        Specifies if only the antiparallel orientation should be concidered in alignment. Defaults to False


    '''
    PNA_map = k_mapper(target)
    PNA_map.add_sequence(PNA)
    if not antiparallel_only:
        PNAr = PNA[::-1]
        PNA_map.add_sequence(PNAr)
    PNA_map.map_kmers()
    return(PNA_map.target_match)

class InputError(Exception):
    """Raised when user input arguments are incorrect"""
    pass


def _XOR(a,b):
    return(bool(a) != bool(b))

class PNA_Designer:
    '''
    A tool used to streamline the design of PNA oligos used for blocking
    amplificaiton of specific DNAs during PCR amplificaiton 
    
    PNA_designer is a module that automates key steps in designing a PNA blocker. Unfortunatly,
    no single PNA is universal in blocking DNA, thus, each requires testing in order to determine
    the sequence space being blocked during amplificaiton. This tool automates inital *in silico* 
    analysis required for development of a new PNA. To complement automatic analysis, this tool
    also makes data access more convient by utilizing RNA sequence data avalible
    at the `Silva database <https://www.arb-silva.de/>`_. Silva identifiers can be used to 
    specify either the PNA target sequence or sequences concidered in the analysis. A copy of silva
    datasets will be download and then compiled. 

    Parameters
    ----------
    results_file: str
        The name of the output file
    target_silva_accession : str, optional
        Silva accession of PNA target. Required if "target_fastafile" is not provided.
    target_fastafile: str, optional
        Fasta file with PNA target. File should only have PNA target sequence. Required if
        "target_silva_accession" is not provided.
    sequence_file: str, optional
        Fasta file used to desgin the PNA. Required if "sequence_silva_taxid" not provided.
    sequence_silva_taxid: int, optional
        Silva taxid used to desgin the PNA. Required if "sequence_silva_taxid" not provided.
    primer_F: str, optional
        Forward primer used for the in silico amplification. If not specified, the full sequence
        of the PNA target and sequences for PNA design will be used. 
    primer_R: str, optional
        Reverse primer used for the in silico amplification. If not specified, the full sequence
        of the PNA target and sequences for PNA design will be used. 
    kmer_range: tupple, optional
        Kmer range used in PNA design. Defaults range: 9-14
    silva_dataset: str, optional
        Silva Database used to fetch seqeucnes. Avalible options are **parc** (all sequences),
        **ref** (high quality sequences), or **nr** (non-redundant, clustered 
        99% identity criterion). Defaults to **nr**
    subunit: str, optional
        Specify **large** or **small** ribosomal subunit in silva database. Defaults to **small**
    database_dir: str, optional
        Local directory of silva database if "target_silva_accession" or 
        "sequence_silva_taxid" is specified. Defualt locaiton is "~/.pna_designer"
        
    Attributes
    ----------
    failed_amplification : list of str
        If primers are specified, accession's of sequences that do not have an exact match
        to the primers will be recorded

    '''

    def __init__(self,result_file=str(),target_silva_accession=str(),target_fastafile=str(),
    sequence_file=str(),sequence_silva_path=int(),primer_F=str(),primer_R=str(),kmer_range=(9,14),
    silva_dataset='ref',subunit='small',database_dir=str()):

        self.target = () #tuple of (accession,seq)
        self.primer_F = primer_F
        self.primer_R = primer_R
        self.failed_amplification=[]
        self.krange = kmer_range
        
        #raise an error if the correct sequencing data was not provided
        #both or missing or both a provided throws an error
        if not _XOR(target_silva_accession,target_fastafile):
            raise InputError("Please provide either a silva accession or a target fasta file")

        if not _XOR(sequence_file,sequence_silva_path):
            raise InputError("Please provide either a fasta file or a silva taxonomic ID")

        if _XOR(self.primer_F,self.primer_R):
            raise InputError("You must provide both Forward and Reverse primers")

        
        #initialize the silva database
        if target_silva_accession or sequence_silva_path:
            self.silva = silva_manager(location=database_dir,dataset=silva_dataset, subunit=subunit)
        
        
        #getting the target sequence
        if target_fastafile:
            for i,line in enumerate(self._iter_fasta(target_fastafile)):
                self.target = line
                if i > 0:
                    raise InputError("You can only design a PNA for a single sequence or consensus sequence")
        #no target_fasta, get silva
        else:
            dat=self.silva.get_seqs(target_silva_accession)[0]
            self.target = (dat[0],dat[1])

        #amplify target if needed
        if primer_F and primer_R:
            print("Amplifying PNA target")
            target_amplicon=InSilico_PCR.sim_amplify(self.primer_F,self.primer_R,self.target[1])
            self.target = (self.target[0],target_amplicon)

        #initialize the k-mer mapper
        self.pna = k_mapper(target_sequence=self.target[1],krange = kmer_range)
        
        #Lets starting adding sequences to make our PNA
        if sequence_silva_path:
            sequence_iterator = self._iter_silvatax(sequence_silva_path)
        else:
            sequence_iterator = self._iter_fasta(sequence_file)
        
        #if amplifying targets 
        if self.primer_F and self.primer_R:
            print("Amplifying and Collecting sequence K-mers")
            for accession,sequence in sequence_iterator:
                try:
                    seq = InSilico_PCR.sim_amplify(self.primer_F,self.primer_R,sequence)
                    self.pna.add_sequence(seq)
                except InSilico_PCR.PrimerError:
                    self.failed_amplification.append(accession)
            if self.failed_amplification:
                print('\tWARNING: {} sequences failed to amplify, please see attribute "failed_amplification"'.format(str(len(self.failed_amplification))))
        else:
            print("Collecting sequence K-mers")
            for accession,sequence in sequence_iterator:
                self.pna.add_sequence(sequence)

        #mapping kmers to target
        print("Mapping Kmers")
        self.pna.map_kmers()
        print("Mapping Complete")
        self.pna.write_results(result_file)
        print("Results have been written to {}".format(result_file))

    def _iter_silvatax(self,parent_taxpath):
        accessions = self.silva.get_accessions(parent_taxpath)
        dat = self.silva.get_seqs(accessions)
        for row in dat:
            yield tuple(row)

    def _iter_fasta(self,file):
        f=open(file,'r')
        seq=''
        key=''
        for line in f:
            if line[0] == '>':
                if key!='':
                    yield ((key,seq))
                key=line.strip()
                seq=''
            else:
                temp=line.strip()#removing end of line
                temp=temp.replace(' ',' ')#removing any possible spaces
                seq+=temp
        if key:#just incase file is empty
            yield ((key,seq))
        f.close()

    def map_PNA(self,PNA,antiparallel_only=False):
        '''
        Align a PNA to a given target sequence

        map_PNA will use the k_mapper class to generate k-mer alignments to target DNA. PNA oligomers
        can bind to DNA in either orientation, thus, the reverse of the PNA is also searched for alignments

        Parameters
        ----------
        PNA : str
            PNA sequence
        antiparallel_only : bool, optional
            Specifies if only the antiparallel orientation should be concidered in alignment. Defaults to False

        '''
        
        dat=map_PNA(self.target[1],PNA,krange=self.krange,antiparallel_only=antiparallel_only)
        results = self.pna.get_results()
        results['PNA mapping']=dat
        return(results)