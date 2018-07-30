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

class silva:
    def __init__(self,location=str(),dataset='ref', subunit='small'):
        if location:
            self.location=location
        else:
            self.location = os.path.join(os.path.expanduser('~'),'.pna_designer')
            if not os.path.exists(self.location):
                os.makedirs(self.location)
        
        # % is replaced with either parc for all sequences or ref for high quality sequences
        # * is wildcard for silva release number
        ftp = 'ftp.arb-silva.de/current/Exports/'
        ftp_tax = ftp+'taxonomy/'
                       
        if subunit not in ('small','large'):
            raise ValueError("Invalid subunit, must be either small (16S,18S), or large (23S,28S)")
        
        if dataset == 'ref': 
            fasta_file = 'SILVA_*_SSU%_tax_silva.fasta.gz'.replace('%','Ref')
            taxmap_file = 'taxmap_slv_ssu_%_*.txt.gz'.replace('%','ref')
            tax_file = 'tax_slv_ssu_*.txt'       

        elif dataset == 'nr':
            fasta_file = 'SILVA_*_SSURef_Nr99_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_ref_nr_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt'
        
        elif dataset == 'parc':
            fasta_file = 'SILVA_*_SSU%_tax_silva.fasta.gz'.replace('%','Parc')
            taxmap_file = 'taxmap_slv_ssu_%_*.txt.gz'.replace('%','parc')
            tax_file = 'tax_slv_ssu_*.txt'  
        
        else:
            raise ValueError("Invalid dataset, choose one of the following: ref, nr, parc")             
        
        if subunit == 'large':
            fasta_file=fasta_file.replace('SSU','LSU')
            taxmap_file=taxmap_file.replace('ssu','lsu')
            tax_file=tax_file.replace('ssu','lsu')

        self.release = self._find_release()
        fasta_file=fasta_file.replace('*',self.release)
        taxmap_file=taxmap_file.replace('*',self.release)
        tax_file=tax_file.replace('*',self.release)
        self.fasta_file=os.path.join(self.location,fasta_file)
        self.taxmap_file=os.path.join(self.location,taxmap_file)
        self.tax_file=os.path.join(self.location,tax_file)

        #lets check if these files exsist
        if not os.path.exists(self.fasta_file):
            self._fetch_file(ftp+fasta_file)
        if not os.path.exists(self.taxmap_file):
            self._fetch_file(ftp_tax+taxmap_file)
        if not os.path.exists(self.tax_file):
            self._fetch_file(ftp_tax+tax_file)
    



    def _find_release(self):
        #finding the current release
        DE = ftp.Download_Engine('ftp.arb-silva.de')
        file_list = DE.listdir('current/Exports/')
        for f in file_list:
            if '_SSURef_tax_silva.fasta.gz' in f:
                return(f.split('_')[1])

    def _fetch_file(self,target_file):
        DE = ftp.Download_Engine('ftp.arb-silva.de',destination = self.location)
        DE.add_target(target_file)
        DE.download()
    