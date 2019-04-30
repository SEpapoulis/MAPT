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
from . import ftp
from .database.db import silva_db
from .database import parser

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
        self.db_file = os.path.join(self.location,'sql.db')
        if not os.path.exists(self.db_file):
            #create new database and add tables
            self.db = silva_db(self.db_file,self.tree_file)
            self.db.add_fasta(self.fasta_file)
            self.db.add_taxmap(self.taxmap_file)
            self.db.add_tax(self.tax_file)
            self.db.add_embl(self.embl_file)
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

    


