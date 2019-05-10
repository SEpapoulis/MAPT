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


from .sql_core import sqldb
from .parser import gzip_fasta,gzip_table
from .TreeOfLife import load_ToL
from . import ftp
import os

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

#TODO
#API for DNA to block and lookup funcitons for DNA    
#Put the path building into its own module, it makes the core too messy
#Kmer analysis- might want to try and do a cython implimitation to get it fast


class silva_manager:
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


#TODO: taxmap and embl are redundant!!!! merge the two into one spreadsheet
#       Ideally, put ncbi_taxid as a column in taxmap

class silva_db():
    def __init__(self,dbfile,tree_file=str(),fasta_file=str(),
    taxmap_file=str(),tax_file=str(),embl_file=str()):

        self.db = sqldb(dbfile)

        #Tree of life is a dictionary of nodes
        #nodes have references to parent and child(ren) nodes

        
        #NOTE: Our accessions combine the International Nucleotide Sequence Database
        #      Collaboration (INSDC) primary accession with the start and stop sites.
        #      Thus, accessions = INSDC primary accession.start.stop      
        #      
        
        #accession is primary key
        self.taxmap_columns = (('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('taxid','INTEGER'),
                                ('ncbi_taxid','INTEGER'),
                                ('seq','TEXT'))

        #taxid primary key
        self.tax_columns = (('path', 'TEXT'), ('taxid', 'INTEGER'), 
                                ('rank','TEXT'), ('remark', 'TEXT'), 
                                ('release','INTEGER'))

        #Adding data
        if tree_file:
            self.ToL=load_ToL(tree_file)    
        if self._setup_required():
            self._compile_database(fasta_file,taxmap_file,tax_file,embl_file)

        #used as a error check for inputs
        if self.db.table_exists('tax'):
            self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])     

    def _setup_required(self):
        
        return(not self.db.table_exists('taxmap'))

    def _compile_database(self,fasta_file,taxmap_file,tax_file,embl_file):
        print("Compiling Silva database, this make take some time")
        if fasta_file:
            print('\tAdding sequences to database....', end='\r')
            self.add_fasta(fasta_file)
            print('\tAdding sequences to database....  [DONE]')
        print('\tAdding taxonmic metadata to database....', end='\r')
        if taxmap_file:
            self.add_taxmap(taxmap_file)
        if tax_file:
            self.add_tax(tax_file)
        print('\tAdding taxonmic metadata to database....  [DONE]')
        print('\tFinishing compliation of a local silva database....', end='\r')
        if embl_file:
            self.add_embl(embl_file)
            if self.db.table_exists('taxmap'):
                self.db.copycol_bypk('taxmap','embl','accession','ncbi_taxid')
                self.db.drop_table('embl')
                self.db.copycol_bypk('taxmap','seqs','accession','seq')
                self.db.drop_table('seqs')
        print('\tFinishing compliation of local silva database....  [DONE]')
        print("Compliation Complete")

    #this funciton will add our fasta data into the sql database
    #700000 sequences at a time (as to not overflow RAM)
    def add_fasta(self,fasta_file,write_size=700000,include_path=False):
        if include_path:
            seqs_columns = (('accession', 'TEXT'),
                                ('path', 'TEXT'),
                                ('seq','TEXT'))
        else:
            seqs_columns = (('accession', 'TEXT'),
                                ('seq','TEXT'))
        self.db.drop_table('seqs')
        #accession is primary key
        self.db.create_table('seqs',columns=seqs_columns,primary_key=['accession'])
        seqs_col=[el[0] for el in seqs_columns]
        i = 1
        to_db=[]
        for header,sequence in gzip_fasta(fasta_file):
            dat = header[1:].split(' ') #remove '>' and split by ' '
            entry=[dat[0]]
            if include_path:
                path = ' '.join(dat[1:]) #rebuild path
                entry.append(path)
            entry.append(sequence.replace('U','T'))
            to_db.append(entry)
            if i % write_size == 0:
                self.db.insert_rows('seqs',seqs_col,to_db)
                to_db = []
            i+=1
        #adding remaining sequences
        if to_db:
            self.db.insert_rows('seqs',seqs_col,to_db)


    def add_taxmap(self,taxmap_file):
        taxmap_columns = (('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('taxid','INTEGER'))
        self.db.drop_table('taxmap')
        self.db.create_table('taxmap',columns=taxmap_columns,
                            primary_key=['accession'])
        taxmap_cols=[el[0] for el in taxmap_columns]
        to_db=[]
        for line in gzip_table(taxmap_file):
            accession = '.'.join(line[0:3])
            to_db.append([accession]+line[3:])
        self.db.insert_rows('taxmap',taxmap_cols,to_db)
        
    
    def add_tax(self,tax_file):
        self.db.drop_table('tax')
        self.db.create_table('tax',columns=self.tax_columns,primary_key=['taxid'])
        tax_cols=[el[0] for el in self.tax_columns]
        to_db=[]
        f = open(tax_file,'r')
        for line in f:
            dat = line.strip().split('\t')
            if len(dat) == 3:
                dat.extend(['',''])
            to_db.append(dat)
        f.close()
        self.db.insert_rows('tax',tax_cols,to_db)
        #if we add tax, we must construct _ranks
        self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])
            
    def add_embl(self,embl_file):
        embl_columns = (('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('ncbi_taxid','INTEGER'))
        self.db.drop_table('embl')
        self.db.create_table('embl',columns=embl_columns,
                            primary_key=['accession'])
        embl_cols=[el[0] for el in embl_columns]
        to_db=[]
        for line in gzip_table(embl_file):
            accession = '.'.join(line[0:3])
            to_db.append([accession]+line[3:])
        self.db.insert_rows('embl',embl_cols,to_db)

    def get_organism_byname(self,name):
        name=name.lower()
        dat=self.db.get_columns(table_name='taxmap',
                                columns=['organism_name','accession','path','taxid'])
        found = []
        for row in dat:
            try:
                if name in row['organism_name'].lower():
                    found.append(row)
            except IndexError:
                print(row)
        if found:
            print('Matches:')
            for row in found:
                print('\n\tName:       {}\n\tAccession:  {}\n\tTaxon:      {}\n\ttaxid:      {}'
                    .format(row['organism_name'],row['accession'],row['path'],row['taxid']))
        else:
            print('Sorry, no matches\n')
        return(found)

    def get_taxon_byname(self,name):
        name=name.lower()
        dat=self.db.get_columns(table_name='tax',columns=['path','taxid','rank'])
        found = list()
        for row in dat:
            if name in row['path'].split(';')[-2].lower():
                found.append(row)
        if found:
            print('Matches:')
            for row in found:
                
                print('\n\tTaxon:  {}\n\ttaxid:  {}\
                    \n\trank:   {}'.format(row['path'],row['taxid'],row['rank']))
        else:
            print('Sorry, no matches\n')
        return(found)        

    #this will use the tree of life to get all descendant taxids
    def _recursive_descendant_taxids(self,taxid):
        taxids=[]
        node = self.ToL[taxid]
        for childnode in node.children:
            if childnode.children:
                taxids.extend(self.get_descendant_taxids(childnode.taxid))
            else:
                taxids.append(childnode.taxid)
        return(taxids)

    def get_descendant_taxids(self,taxid):
        taxids = self._recursive_descendant_taxids(taxid)
        taxids.append(taxid)
        return(taxids)

    def get_accessions(self,taxid):
        accessions=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn='taxid',
                                        values=taxid,columns='accession')
        return([el['accession'] for el in accessions])

    def get_accession_dat(self,accession,get_columns=None):
        if get_columns==None:
            get_columns = [el[0] for el in self.taxmap_columns]

        accessions=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn='accession',
                                        values=accession,columns=get_columns)
        return(accessions)

    def write_fasta(self,accessions,fasta_file,include_taxid=False,include_path=False,
                   include_organism_name=False,include_NCBItaxid=False):
        cols = ['accession','seq']
        taxmap_metadata=[]
        if include_organism_name:
            taxmap_metadata.append('organism_name')
        if include_path:
            taxmap_metadata.append('path')
        if include_taxid:
            taxmap_metadata.append('taxid')
        if include_NCBItaxid:
            taxmap_metadata.append('ncbi_taxid')
        cols.extend(taxmap_metadata)
        dat = self.get_accession_dat(accessions,get_columns=cols)
        with open(fasta_file,'w') as f:
            for row in dat:
                header=['>'+row['accession']]
                for metadata in taxmap_metadata:
                    header.append("[{}:{}]".format(str(metadata),str(row[metadata])))
                f.write(' '.join(header))
                seq = row['seq']
                for i in range(0,len(seq),100):
                    if i+100 <= len(dat):
                        f.write(dat[i:i+100]+'\n')
                    else:
                        f.write(dat[i:len(dat)]+'\n')
        

        


