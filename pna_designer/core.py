#ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobioal.tar.gz

#https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz

'''
The Parc datasets comprise the entire SILVA databases for the respective gene, 
whereas the Ref datasets represent a subset of the Parc comprising only 
high-quality nearly full-length sequences.


tax_slv_[ls]su_VERSION.txt
These files contain taxonomic rank designations for all taxonomic paths
used in the SILVA taxonomies. Additionally, a unique numeric identifier is 
assigned to each taxon (path). These identifiers will be mostly stable in 
upcoming SILVA releases

columns: path, taxid, rank, remark, release

taxmap_TAXNAME_[ls]su_VERSION.txt
----------------------------
mapping of each entry in the SILVA database to a taxonomic path. Different
rRNA regions of the same INSDC entry (genome) may be assigned to multiple
paths (contaminations or micro diversity among the rRNA sequences).

columns: pac, start, stop ,path ,name, taxid
'''


import os
from . import ftp,InSilico_PCR
from .database.db import silva_db
from .database import parser
from collections import defaultdict

#TODO
#API for DNA to block and lookup funcitons for DNA    
#Put the path building into its own module, it makes the core too messy
#Kmer analysis- might want to try and do a cython implimitation to get it fast
class silva:
    def __init__(self,location=str(),dataset='ref', subunit='small'):

        #Establishing the path for database files
        if location:
            self.location=location
        else:
            self.location = os.path.join(os.path.expanduser('~'),'.pna_designer')
            if not os.path.exists(self.location):
                os.makedirs(self.location)
        
        #Currently, we fetch the version number of the current release
        #TODO: add the ability to fetch different versions
        self.version_file = os.path.join(self.location,'version')
        if not os.path.exists(self.version_file):     
            self.version = self._find_release()
        else:
            f = open(self.version_file,'r')
            self.version = f.readline()
            f.close()

        #Generating the filenames for the database
        fasta_file,taxmap_file,tax_file,embl_file,tree_file = self._generate_filenames(dataset,subunit)
        self.fasta_file = os.path.join(self.location,fasta_file)
        self.taxmap_file = os.path.join(self.location,taxmap_file)
        self.tax_file = os.path.join(self.location,tax_file)
        self.embl_file = os.path.join(self.location,embl_file)
        self.tree_file = os.path.join(self.location,tree_file)


        #lets check if these files exsist, if not lets download
        host = 'ftp.arb-silva.de/current/Exports/'
        host_tax = host+'taxonomy/'
        if not os.path.exists(self.fasta_file):
            self._fetch_file(host+fasta_file)
        if not os.path.exists(self.taxmap_file):
            self._fetch_file(host_tax+taxmap_file)
        if not os.path.exists(self.tax_file):
            self._fetch_file(host_tax+tax_file)
        if not os.path.exists(self.embl_file):
            self._fetch_file(host_tax+embl_file)
        if not os.path.exists(self.tree_file):
            self._fetch_file(host_tax+tree_file)


        #we need to complie a database if there is not one already
        #TODO: add except EOFError to downloads to retry failed downloads
        #WARNING: right now, the db is compiled with all or none while silva_db can partially compile db 
        self.db_file = os.path.join(self.location,'sql.db')
        if not os.path.exists(self.db_file):
            #create new database and add tables
            self.db = silva_db(self.db_file,self.tree_file,self.fasta_file,self.taxmap_file,self.tax_file,self.embl_file)
        else:
            self.db = silva_db(self.db_file,self.tree_file)

    def _generate_filenames(self,dataset,subunit):
        # * is wildcard for silva release number
        if subunit not in ('small','large'):
            raise ValueError("Invalid subunit, must be either small (16S,18S), or large (23S,28S)")
        elif subunit == 'large' and dataset == 'nr':
            raise ValueError("There is not a silva non-redundant dataset for the large subunit,\nplease specify ref or parc")
        if dataset == 'ref': 
            fasta_file = 'SILVA_*_SSURef_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_ref_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt'       
            embl_file = 'taxmap_embl_ssu_ref_*.txt.gz'
            tree_file = 'tax_slv_ssu_*.tre'
        elif dataset == 'nr':
            fasta_file = 'SILVA_*_SSURef_Nr99_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_ref_nr_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt'
            embl_file = 'taxmap_embl_ssu_ref_nr99_*.txt.gz'
            tree_file = 'tax_slv_ssu_*.tre'
        elif dataset == 'parc':
            fasta_file = 'SILVA_*_SSUParc_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_parc_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt'
            embl_file = 'taxmap_embl_ssu_parc_*.txt.gz'  
            tree_file = 'tax_slv_ssu_*.tre'
        else:
            raise ValueError("Invalid dataset, choose one of the following: ref, nr, parc")             
        
        if subunit == 'large':
            fasta_file=fasta_file.replace('SSU','LSU')
            taxmap_file=taxmap_file.replace('ssu','lsu')
            tax_file=tax_file.replace('ssu','lsu')
            embl_file=embl_file.replace('ssu','lsu')
            tree_file=tree_file.replace('ssu','lsu')
        
        self.release = self._find_release()
        fasta_file=fasta_file.replace('*',self.release)
        taxmap_file=taxmap_file.replace('*',self.release)
        tax_file=tax_file.replace('*',self.release)
        embl_file=embl_file.replace('*',self.release)
        tree_file=tree_file.replace('*',self.release)
        return(fasta_file,taxmap_file,tax_file,embl_file,tree_file)

    def _find_release(self):
        #finding the current release
        DE = ftp.Download_Engine('ftp.arb-silva.de')
        file_list = DE.listdir('current/Exports/')
        for f in file_list:
            if '_SSURef_tax_silva.fasta.gz' in f:
                version = f.split('_')[1]
                f = open(self.version_file,'w')
                f.write(version+'\n')
                f.close()
                return(version)

    def _fetch_file(self,target_file):
        DE = ftp.Download_Engine('ftp.arb-silva.de',destination = self.location)
        DE.add_target(target_file)
        DE.download()

    
class PNA:
    def __init__(self,target_sequence=str,sequences=list(),krange=tuple()):

        self.target= target_sequence.upper()

        #coordiantes for unique kmers that align
        self.target_match_unique = [0]*len(self.target)

        #coordiantes for all kmers that align, where kmers that are found
        #more often in seqlist will be weighted more
        self.target_match = [0]*len(self.target)

        #dictoinary mapping to defaultdicts of kmer counts
        #note, if kdict[k][kmer] == 0, that means this kmer is unique
        #to the target sequence!
        self.kdict=defaultdict(int)

        #Inclusive range
        self.krng=list((range(krange[0],krange[1]+1)))

        #seqlist is a list of sequences that you DO NOT want to block
        #thus, you want the area of target that shares the least number
        #of kmers with seqlist
        self.seqlist = [seqs.upper() for seqs in sequences]

    def get_kmers(self,k,seqlist):
        kmers=[]
        for seq in seqlist:
            seq=seq.lower()
            for i in range(k,len(seq)):
                kmers.append(seq[i-k:i])
        return(set(kmers))

    def add_kmers(self,k):
        kmers = defaultdict(int)
        for seq in self.seqlist:
            for i in range(k,len(seq)):
                kmer=seq[i-k:i]
                kmers[kmer]+=1
        self.kdict[k] = kmers
                
    def add_kmerrng(self,start,end):
        print("Gathering kmer range {}-{}".format(str(start),str(end)))
        for k in range(start,end+1):
            print("\tGathering {}-mers".format(str(k)),end='\r')
            self.add_kmers(k)
        print("\t Gathering Complete!")

    def exp_add_kmers(self,seq):
        seq=seq.upper()
        for i in range(0,len(seq)):
            for k in self.krng:
                if i+k <len(seq):
                    self.kdict[seq[i:i+k]]+=1
                

    def map_kmers(self):
        for k in self.kdict:
            for i in range(k,len(self.target)+1):
                #seq_kmer is the number of kmers in seqlist
                seq_kmer = self.kdict[k][self.target[i-k:i]]
                if seq_kmer:
                    for align in list(range(i-k,i+1)):
                        self.target_match_unique[align]+=1
                        self.target_match[align]+=seq_kmer

    def write_results(self,write_file=str()):
        if not write_file:
            write_file=os.path.join(os.getcwd(),'KmerCoordiantes.csv')
        with open(write_file,'w') as o:
            o.write('Nuclotide,index,unique match, absolute match')
            for i in range(0,len(self.target_match)):
                o.write(','.join([self.target[i],str(i),self.target_match_unique[i],self.target_match[i]]))

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
    silva_dataset='ref',subunit='small',database_location=str()):

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
            self.silva = silva(location=database_location,dataset=silva_dataset, subunit=subunit)

        #getting the target sequence
        if target_fastafile:
            for i,line in enumerate(self.iter_fasta(target_fastafile)):
                self.target = line
                if i > 0:
                    raise InputError("You can only design a PNA for a single sequence or consensus sequence")
        #no target_fasta, get silva
        else:
            dat=self.silva.db.get_sequence_dat(target_silva_accession,['accession','seq'])[0]
            self.target = (dat['accession'],dat['seq'])

        #amplify target if needed
        if primer_F and primer_R:
            target_amplicon=InSilico_PCR.sim_amplify(self.primer_F,self.primer_R,self.target[1])
            self.target = (self.target[0],target_amplicon)

        self.pna = PNA(target_sequence=self.target[1],krange = kmer_range)
        
        

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
        
