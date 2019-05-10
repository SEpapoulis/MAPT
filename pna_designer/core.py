import os
from . import InSilico_PCR
from .database.silva import silva_manager
from collections import defaultdict


    
class PNA:
    def __init__(self,target_sequence=str,sequences=list(),krange=tuple()):
        #Inclusive range
        self.krng=list((range(krange[0],krange[1]+1)))

        self.target= target_sequence.upper()
        self.target_klocation = self.get_target_kmermap()
        #coordiantes for unique kmers that align
        self.target_match_unique = [0]*len(self.target)

        #coordiantes for all kmers that align, where kmers that are found
        #more often in seqlist will be weighted more
        self.target_match = [0]*len(self.target)

        #dictoinary mapping to defaultdicts of kmer counts
        #note, if kdict[kmer] == 0, that means this kmer is unique
        #to the target sequence!
        self.kdict=defaultdict(int)


    def get_target_kmermap(self):
        kmer_location={}
        for k in self.krng:
            for i in range(k,len(self.target)+1):
                start=i-k
                stop =i 
                kmer = self.target[start:stop]
                if kmer in kmer_location:
                    kmer_location[kmer].append((start,stop))
                else:
                    kmer_location[kmer] = [range(start,stop)]
        return(kmer_location)

    def add_kmers(self,seq):
        seq=seq.upper()
        for i in range(0,len(seq)):
            for k in self.krng:
                if i+k <len(seq):
                    self.kdict[seq[i:i+k]]+=1
                

    def map_kmers(self):
        for kmer in self.kdict:
            if kmer in self.target_klocation:
                locations = self.target_klocation[kmer]
                for location in locations:
                    loc = list(location)
                    for i in loc:
                        self.target_match_unique[i]+=1
                        self.target_match[i]+=self.kdict[kmer]

    def write_results(self,write_file=str()):
        if not write_file:
            write_file=os.path.join(os.getcwd(),'KmerCoordiantes.csv')
        else:
            write_file=os.path.join(write_file)
        with open(write_file,'w') as o:
            o.write('Nuclotide,index,unique match, absolute match')
            for i in range(0,len(self.target_match)):
                o.write('\n')
                o.write(','.join([self.target[i],str(i),
                str(self.target_match_unique[i]),str(self.target_match[i])]))
                

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message    
    """
    pass


def _XOR(a,b):
    return(bool(a) != bool(b))

class Designer:
    def __init__(self,target_silva_accession=str(),target_fastafile=str(),
    sequence_file=str(),sequence_silva_taxid=int(),primer_F=str(),primer_R=str(),kmer_range=(9,14),
    silva_dataset='ref',subunit='small',database_dir=str(),results_file=str()):

        self.target = () #tuple of (accession,seq)
        self.primer_F = primer_F
        self.primer_R = primer_R
        self.failed_amplification=[]

        #raise an error if the correct sequencing data was not provided
        #both or missing or both a provided throws an error
        if not _XOR(target_silva_accession,target_fastafile):
            raise InputError("Please provide either a silva accession or a target fasta file")

        if not _XOR(sequence_file,sequence_silva_taxid):
            raise InputError("Please provide either a fasta file or a silva taxonomic ID")

        if _XOR(self.primer_F,self.primer_R):
            raise InputError("You must provide both Forward and Reverse primers")

        #initialize the silva database
        if target_silva_accession or sequence_silva_taxid:
            self.silva = silva_manager(location=database_dir,dataset=silva_dataset, subunit=subunit)

        #getting the target sequence
        if target_fastafile:
            for i,line in enumerate(self.iter_fasta(target_fastafile)):
                self.target = line
                if i > 0:
                    raise InputError("You can only design a PNA for a single sequence or consensus sequence")
        #no target_fasta, get silva
        else:
            dat=self.silva.db.get_accession_dat(target_silva_accession,['accession','seq'])[0]
            self.target = (dat['accession'],dat['seq'])

        #amplify target if needed
        if primer_F and primer_R:
            target_amplicon=InSilico_PCR.sim_amplify(self.primer_F,self.primer_R,self.target[1])
            self.target = (self.target[0],target_amplicon)

        self.pna = PNA(target_sequence=self.target[1],krange = kmer_range)
        
        #Lets starting adding sequences to make our PNA
        if sequence_silva_taxid:
            sequence_iterator = self.iter_silvatax(sequence_silva_taxid)
            self.test=sequence_iterator
        else:
            sequence_iterator = self.iter_fasta(sequence_file)
        
        #if amplifying targets 
        if self.primer_F and self.primer_R:
            print("Amplifying and Collecting sequence Kmers")
            for accession,sequence in sequence_iterator:
                try:
                    seq = InSilico_PCR.sim_amplify(self.primer_F,self.primer_R,sequence)
                    self.pna.add_kmers(seq)
                except InSilico_PCR.PrimerError:
                    self.failed_amplification.append(accession)
        
        #mapping kmers to target
        print("Mapping Kmers")
        self.pna.map_kmers()
        self.pna.write_results(results_file)
                

    def iter_silvatax(self,parent_taxid):
        taxids = self.silva.db.get_descendant_taxids(parent_taxid)
        dat = self.silva.db.get_accession_dat(self.silva.db.get_accessions(taxids),['accession','seq'])
        for row in dat:
            yield tuple(row)

    def iter_fasta(self,file):
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
        
