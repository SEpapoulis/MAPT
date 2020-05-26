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
import gzip

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
#add the ability to fetch different versions
#API for DNA to block and lookup funcitons for DNA    
#Put the path building into its own module, it makes the core too messy
#Kmer analysis- might want to try and do a cython implimitation to get it fast
#add except EOFError to downloads to retry failed downloads

class silva_manager:
    '''
    Manager class used to search and retrieve Silva sequences

    The Silva manager class is used by the PNA_Designer to retrieve sequences
    under a specific taxonomic unit, however, users can access the silva_manager API
    for their own uses.
    
    Parameters
    ----------
    location: str
        Directory location for local database copy. Defualt locaiton is “~/.pna_designer”
    dataset : str
        Avalible options are **parc** (all sequences),
        **ref** (high quality sequences), or **nr** (non-redundant, clustered 
        99% identity criterion). Defaults to **nr**
    subunit: tuple of ints, optional
        Specify **large** or **small** ribosomal subunit in silva database. Defaults to **small**
    release : int, optional
        Specify a Silva release. Default = Current. Note that only release 132 or greater are avlible
    
    '''

    def __init__(self,location=str(),dataset='ref', subunit='small',release=int()):

        #Establishing the path for database files
        if location:
            self.location=location
        else:
            self.location = os.path.join(os.path.expanduser('~'),'.MAPT')
            if not os.path.exists(self.location):
                os.makedirs(self.location)
        
        #Currently, we fetch the version number of the current release
        #TODO: add the ability to fetch different versions

        #Generating the filenames for the database
        self._meta_dat = {}
        self._load_meta(release=str(release))
        self.__generate_filenames(dataset,subunit)
        if os.path.exists(self.db_file):
            self.db = silva_db(self.db_file,self.tree_file)
        else:
            self.__build_protocol()
            

    def __build_protocol(self):
        if self._meta_dat['current'] == 'True':
            host = 'ftp://ftp.arb-silva.de/current/Exports/'
        else:
            host = 'ftp://ftp.arb-silva.de/release_{}/Exports/'.format(self._meta_dat['release'])
        host_tax = host+'taxonomy/'
        if not os.path.exists(self.fasta_file):
            self._fetch_file(('fasta',host+os.path.basename(self.fasta_file)))
        if not os.path.exists(self.taxmap_file):
            self._fetch_file(('taxmap',host_tax+os.path.basename(self.taxmap_file)),False)
        if not os.path.exists(self.tax_file):
            self._fetch_file(('tax',host_tax+os.path.basename(self.tax_file)),False)
        if not os.path.exists(self.tree_file):
            self._fetch_file(('tree',host_tax+os.path.basename(self.tree_file)),False)

        self.db = silva_db(self.db_file,self.tree_file,self.fasta_file,self.taxmap_file,self.tax_file)

        #cleanup
        
        os.remove(self.fasta_file)
        os.remove(self.taxmap_file)
        os.remove(self.tax_file)
        #os.remove(self.embl_file)

    def __generate_filenames(self,dataset,subunit):
        ''' This function is not compatible with version 132, 132 did not have tax or tree
        files as gz files. for compatibility, please update
        '''
        # * is wildcard for silva release number
        if subunit not in ('small','large'):
            raise ValueError("Invalid subunit, must be either small (16S,18S), or large (23S,28S)")
        elif subunit == 'large' and dataset == 'nr':
            raise ValueError("There is not a silva non-redundant dataset for the large subunit,\nplease specify ref or parc")
        if dataset == 'ref': 
            fasta_file = 'SILVA_*_SSURef_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_ref_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt.gz'       
            tree_file = 'tax_slv_ssu_*.tre.gz'
        elif dataset == 'nr':
            fasta_file = 'SILVA_*_SSURef_Nr99_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_ref_nr_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt.gz'
            tree_file = 'tax_slv_ssu_*.tre.gz'
        elif dataset == 'parc':
            fasta_file = 'SILVA_*_SSUParc_tax_silva.fasta.gz'
            taxmap_file = 'taxmap_slv_ssu_parc_*.txt.gz'
            tax_file = 'tax_slv_ssu_*.txt.gz'
            tree_file = 'tax_slv_ssu_*.tre.gz'
        else:
            raise ValueError("Invalid dataset, choose one of the following: ref, nr, parc")             
        dbsu="ssu"
        if subunit == 'large':
            fasta_file=fasta_file.replace('SSU','LSU')
            taxmap_file=taxmap_file.replace('ssu','lsu')
            tax_file=tax_file.replace('ssu','lsu')
            tree_file=tree_file.replace('ssu','lsu')
            dbsu="lsu"

        ##self.version_file= os.path.join(self.location,'version')
        r= self._meta_dat['release']
        self.db_file = os.path.join(self.location,'silva{}_{}_{}.db'.format(r,dataset,dbsu))
        self.fasta_file=os.path.join(self.location,fasta_file.replace('*',r))
        self.taxmap_file=os.path.join(self.location,taxmap_file.replace('*',r))
        self.tax_file=os.path.join(self.location,tax_file.replace('*',r))
        self.tree_file=os.path.join(self.location,tree_file.replace('*',r))

        #posthoc changes for 132 compatibility
        if r == '132':
            self.tree_file = self.tree_file.replace('.gz','')
            self.tax_file = self.tax_file.replace('.gz','')

    def _find_release(self):
        #finding the current release  
        file_list = ftp.list_dir('ftp.arb-silva.de/current/Exports/')
        for f in file_list:
            if '_SSURef_tax_silva.fasta.gz' in f:
                release = f.split('_')[1]
                return(release)
        raise FileNotFoundError('Coult not find release version in ftp.arb-silva.de/current/Exports/')

    def _fetch_file(self,target_file,print_progress=True):
        DE = ftp.Download_Engine([target_file],destination = self.location)
        DE.Download(print_progress)

    def _write_meta(self,metafile):
        with open(metafile,'w') as f:
            for field in self._meta_dat:
                f.write(','.join([str(field),str(self._meta_dat[field])])+'\n')

    def _load_meta(self,release):
        current = None
        if release == 0:
            current = self._find_release()
            release=current
        meta_file = os.path.join(self.location,'metadat_{}.csv'.format(release))
        if os.path.exists(meta_file):
            with open(meta_file,'r') as f:
                for line in f:
                    field,dat = line.strip().split(',')
                    self._meta_dat[field]=dat
        else:
            if not current:
                current = self._find_release()
            if release == current:# ==0 for default
                self._meta_dat['release']=current
                self._meta_dat['current']='True'
            else:
                self._meta_dat['release']=release
                self._meta_dat['current']='False'

            self._write_meta(meta_file)

    def find_taxpath(self,taxon_name,print_results=True):
        '''Search for taxonomic unit in Silva database

        Finds information on the given taxonomic name provided. Spelling must be exact to taxonomic name
        to find results.

        Parameters
        ----------
        taxon_name: str
            name of taxon to be searched. Case insensitive
        print_results: bool
            All entries found will be printed with additional metadata. Default True.

        Returns
        -------
        tuple
            returns a tuple containing the db entry and the number of sequences belonging to that taxon
        '''
        dat = self.db.tax_lookup('search',taxon_name.upper(),['path','taxid','rank'])

        out = []
        while dat:
            entry = dat.pop()
            out.append((entry,len(self.get_accessions(entry['path']))))
        return(out)
            

    def find_organism(self,name,print_results=True):
        '''Search for an organism in the Silva database

        You can directly search for a organism by the name. Be aware that only eact matches to the organism name will be shown.

        Parameters
        ----------
        name: str
            name of the organism. Case sensitive. 
        print_results: bool
            All entries found will be printed with additional metadata. Default True.

        Returns
        -------
        list
            List of silva entries
        '''
        dat = self.db.taxmap_lookup('organism_name',name,['organism_name','accession','path'])
        return(dat)

    def get_accessions(self,taxpath):
        '''Retrieve silva accessions under a specific taxonomic unit

        For a given silva taxonomic path, all accessions under this taxon will be retrieved.

        Parameters
        ----------
        taxpath: str
            Desired Silva taxpath. For example, Fungi would be: "Eukaryota;Opisthokonta;Nucletmycea;Fungi;".
            Use find_taxpath() to search the database for a specific taxpath, or use the silva website

        Returns
        -------
        accessions
            List of silva accessions under the given taxonomic path
        '''
        taxid = self.db.tax_lookup('path',taxpath,retrieve=['taxid'])
        if taxid:
            taxids = self.db.get_descendant_taxids(taxid[0]['taxid'])
            accessions = self.db.taxmap_lookup('taxid',taxids,retrieve=['accession'])
            return([el['accession'] for el in accessions])
        return(list())
    
    def get_seqs(self,accessions):
        '''Retrieve sequences from input accessions

        For a given list of accessions, retrieve documented sequences

        Parameters
        ----------
        accessions: list of accessions
            list of accessions used for sequence retrieval

        Returns
        -------
        Sequence list
            A list of tuples(accession,sequence)
        '''
        dat = self.db.taxmap_lookup('accession',accessions,retrieve=['accession','seq'])
        return([tuple([el['accession'],el['seq']]) for el in dat])


#TODO: 
#dropping the tables to then add them to another is a very expensive operation,
#I should just pull them into memeory since they are going to take, at the most 100 meg

class silva_db():
    def __init__(self,dbfile,tree_file=str(),fasta_file=str(),
    taxmap_file=str(),tax_file=str()):

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
                                #('ncbi_taxid','INTEGER'), They broke the file >:(
                                ('seq','TEXT'))

        #taxid primary key
        self.tax_columns = (('path', 'TEXT'), ('taxid', 'INTEGER'), 
                                ('rank','TEXT'), ('remark', 'TEXT'), 
                                ('release','INTEGER'),
                                ('search','TEXT'))

        #Adding data
        if tree_file:
            self.ToL=load_ToL(tree_file)    
        if self._setup_required():
            self._compile_database(fasta_file,taxmap_file,tax_file)

        #used as a error check for inputs
        if self.db.table_exists('tax'):
            self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])     

    def _setup_required(self):
        
        return(not self.db.table_exists('taxmap'))

    def _compile_database(self,fasta_file,taxmap_file,tax_file):
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
        print('\tFinishing compilation of a local silva database....', end='\r')
        self.db.copycol_bypk('taxmap','seqs','accession','seq')
        self.db.drop_table('seqs')
        print('\tFinishing compilation of local silva database....  [DONE]')
        print("Compilation Complete\n")

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
        self.db.create_table('tax',columns=self.tax_columns,primary_key=['path','taxid'])
        tax_cols=[el[0] for el in self.tax_columns]
        to_db=[]
        if '.gz' in tax_file:
            f = gzip.open(tax_file,'r')
        else:
            f = open(tax_file,'r')
        for line in f:
            if '.gz' in tax_file:
                dat = str(line,'utf-8').strip().split('\t')
            else:    
                dat = line.strip().split('\t')
            if len(dat) == 3:
                dat.extend(['',''])
            try:
                dat.append(dat[0].split(';')[-2].upper())#adding to the search column
            except IndexError:
                print(dat)
                raise
            to_db.append(dat)
        if '.gz' not in f:
            f.close()
        self.db.insert_rows('tax',tax_cols,to_db)
        #if we add tax, we must construct _ranks
        self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])
    '''      
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
    '''
    '''
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
    '''
    def find_taxon(self,name):
        name=name.upper()
        dat=self.db.get_rows_bycolvalue(table_name='tax',valuecolumn='search',
                                        values=name,columns=['path','taxid','rank','search'])
        return(dat)        

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

    def tax_lookup(self,valuecolumn,values,retrieve=None):
        if retrieve==None:
            retrieve = [el[0] for el in self.tax_columns]

        dat=self.db.get_rows_bycolvalue(table_name='tax',valuecolumn=valuecolumn,
                                        values=values,columns=retrieve)
        return(dat)

    def taxmap_lookup(self,valuecolumn,values,retrieve=None):
        if retrieve==None:
            retrieve = [el[0] for el in self.taxmap_columns]

        dat=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn=valuecolumn,
                                        values=values,columns=retrieve)
        return(dat)

    def write_fasta(self,accessions,fasta_file,include_taxid=False,include_path=False,
                   include_organism_name=False):
        cols = ['accession','seq']
        taxmap_metadata=[]
        if include_organism_name:
            taxmap_metadata.append('organism_name')
        if include_path:
            taxmap_metadata.append('path')
        if include_taxid:
            taxmap_metadata.append('taxid')
        #if include_NCBItaxid:
        #    taxmap_metadata.append('ncbi_taxid')
        cols.extend(taxmap_metadata)
        dat = self.taxmap_lookup('accessions',accessions)
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
        

        


