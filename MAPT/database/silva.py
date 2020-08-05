from .sql_core import sqldb
from .parser import gzip_fasta,gzip_table,compress_gaps,gzip_fasta_batch
from .TreeOfLife import load_ToL
from . import ftp
import os
import gzip

from multiprocessing import Pool



def compression_worker(element):
    header,sequence=element
    dat = header[1:].split(' ') #remove '>' and split by ' '
    accession,encoding=compress_gaps(dat[0],sequence)
    entry=[accession]
    entry.append(sequence.replace('U','T').replace('.','').replace('-',''))
    entry.append(encoding)
    return(entry)

class silva_manager:
    '''
    Manager class used to search and retrieve Silva sequences

    The Silva manager class is used by the PNA_Designer to retrieve sequences
    under a specific taxonomic unit, however, users can access the silva_manager API
    for their own uses.
    
    Parameters
    ----------
    location: str
        Directory location for local database copy. Default location is “~/.MAPT"
    dataset : str
        Avalible options are **parc** (all sequences),
        **ref** (high quality sequences), or **nr** (non-redundant, clustered 
        99% identity criterion). Defaults to **nr**
    subunit: tuple of ints, optional
        Specify **large** or **small** ribosomal subunit in silva database. Defaults to **small**
    release : int, optional
        Specify a Silva release. Default = Current. Note that only release 132 or greater are available
    
    '''

    def __init__(self,path='',dataset='nr', subunit='small',release=int()):

        if path:
            if path[0] == '~':
                path = os.path.join(os.path.expanduser(path[0]),path[1:].strip('/'))
        else:
            path = os.path.join(os.path.expanduser('~'),'.MAPT')
        #set contact ftp for database
        if not release:
            sil_ftp = 'ftp.arb-silva.de/current/Exports/'
        else:
            sil_ftp = 'ftp://ftp.arb-silva.de/release_{}/Exports/'.format(release)

        fetch,db_name,tree_name,msa_name=self._parse_remotefiles(sil_ftp,dataset,subunit)
        self.release,self.subunit,self.dataset = db_name[:-3].split('_')[1:]
        self.db_file = os.path.join(path,db_name)#contains path
        self.tree_file = os.path.join(path,tree_name)
        self.msa_file = os.path.join(path,msa_name)
        #database outline
        self.taxmap_columns = (('accession', 'TEXT'), 
                        ('path', 'TEXT'),
                        ('organism_name', 'TEXT'), 
                        ('taxid','INTEGER'),
                        #('ncbi_taxid','INTEGER'), They broke the file >:(
                        ('seq','TEXT'))

        self.tax_columns = (('path', 'TEXT'), ('taxid', 'INTEGER'), 
                                ('rank','TEXT'), ('remark', 'TEXT'), 
                                ('release','INTEGER'),
                                ('search','TEXT'))

        self.seq_columns = (('accession', 'TEXT'),('seq','TEXT'))

        if os.path.exists(self.db_file):
            self.db = sqldb(self.db_file)
            self.ToL = load_ToL(self.tree_file)
        else:
            self.db = sqldb(self.db_file)
            self.__build_protocol(fetch)
            self.ToL = load_ToL(self.tree_file)


    def _parse_remotefiles(self,hostftp,dataset,subunit):
        ''' This function will catagorize files based off of parameters'''
      
        #Error checking
        if dataset not in ['nr','parc','ref']:
            raise ValueError("Invalid dataset, must be either parc, ref, or nr")
        if dataset == 'nr':
            dataset = 'ref_nr99'
        if subunit == 'small':
            subunit = 'ssu'
        elif subunit == 'large':
            subunit = 'lsu'
        else:
            raise ValueError("Invalid subunit, must be either small (16S,18S), or large (23S,28S)")

        filelist = []
        #get main directory
        for file in ftp.list_dir(hostftp):
            filelist.append(file)
            if 'SILVA_' in file:
                release=os.path.basename(file).split('_')[1]
        #get taxonomy
        for file in ftp.list_dir(hostftp+'taxonomy/'):
            filelist.append(file)
        #Identify the files needed
        search = {
        'fasta' : 'SILVA_{r}_{s}{d}_tax_silva_trunc.fasta.gz'.format(r=release,s=subunit,d=dataset),
        'msa' : 'SILVA_{r}_{s}{d}_tax_silva_full_align_trunc.fasta.gz'.format(r=release,s=subunit,d=dataset),
        'tree' : 'tax_slv_{s}_{r}.tre'.format(r=release,s=subunit,d=dataset), #for compatibility with 132, .gz is dropped and checked later
        'taxmap' : 'taxmap_slv_{s}_{d}_{r}.txt.gz'.format(r=release,s=subunit,d=dataset.replace('99','')),
        'tax' : 'tax_slv_{s}_{r}.txt'.format(r=release,s=subunit,d=dataset)}#for compatibility with 132, .gz is dropped and checked later
        if dataset == 'ref_nr99':
            dbname='sil_{r}_{s}_nr99.db'.format(r=release,s=subunit)    
        else:
            dbname='sil_{r}_{s}_{d}.db'.format(r=release,s=subunit,d=dataset)
        fetch = {}
        for file in filelist:
            for key in search:
                if key not in fetch:
                    if search[key].lower() == os.path.basename(file).lower():
                        fetch[key]=file
                    if key == 'tree' or key == 'tax':
                        if search[key].lower()+'.gz' == os.path.basename(file).lower():
                            fetch[key]=file
        #adding host for ftp download
        for key in fetch:
            fetch[key] = 'ftp://ftp.arb-silva.de/'+fetch[key]
        return(fetch,dbname,os.path.basename(fetch['tree']),os.path.basename(fetch['msa']))

    def __build_protocol(self,fetch):
        '''protocol to build the silva database via sqlite3'''
        path = os.path.split(self.db_file)[0]

        print('Preparing to Download silva data to {}'.format(path))
        lcl = {}
        for el in fetch:
            lcl_file = os.path.join(path,os.path.basename(fetch[el]))
            if not os.path.exists(lcl_file):
                print('Downloading {}'.format(fetch[el]))
                DE = ftp.Download_Engine([(el,fetch[el])],destination = path)
                DE.Download(print_progress=True)
            lcl[el] = lcl_file

        self.__compile(lcl)
        
        #cleanup
        os.remove(lcl['fasta'])
        os.remove(lcl['taxmap'])
        os.remove(lcl['tax'])
        

    def __compile(self,lcl):
        print("\nCompiling Silva database, this make take some time")
        print('\tAdding taxonmic metadata to database....', end='\r')
        if 'taxmap' in lcl:
            self.__add_taxmap(lcl['taxmap'])
        if 'tax' in lcl:
            self.__add_tax(lcl['tax'])
        print('\tAdding taxonmic metadata to database....  [DONE]')
        if 'fasta' in lcl:
            print('\tAdding sequences to database....',end='\r')
            self.__add_fasta(lcl['fasta'])
            print('\tAdding sequences to database....  [DONE]')
        print("Compilation Complete!\n")

    def __add_fasta_compression(self,fasta_file,tbl='seqs'):#write_size=700000,include_path=False):
        seqs_columns = (('accession', 'TEXT'),
                            ('seq','TEXT'),
                            ('gap_compression','TEXT'))
        #accession is primary key
        self.db.create_table('seqs',columns=seqs_columns,primary_key=['accession'])
        seqs_col=[el[0] for el in seqs_columns]
        #i = 1
        it=0
        to_db=[]
        a_pool = Pool()
        for dat in gzip_fasta_batch(fasta_file,batchsize=16000):
            to_db = a_pool.map(compression_worker,dat)
            self.db.insert_rows('seqs',seqs_col,to_db)
            it += len(to_db)
            print("\t{} sequences added".format(it))
            to_db=[]
            dat=[]

    def __add_fasta(self,fasta_file,write_size=700000,include_path=False):
        '''Adds SILVA_[release]_[subunit][dataset]_tax_silva_trunc.fasta.gz.txt to the database'''
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

    def __add_tax(self,tax_file):
        '''Adds tax_slv_[subunit]_[release].txt to the database'''
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
        #self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])

    def __add_taxmap(self,taxmap_file):
        '''Adds taxmap_slv_[subunit]_[dataset]_[release].txt.gz.txt to the database'''
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
        self.accession_count = len(to_db)

    def _recursive_descendant_taxids(self,taxid):
        '''For each childnode, call child nodes'''
        taxids=[]
        node = self.ToL[taxid]
        for childnode in node.children:
            if childnode.children:
                taxids.extend(self._get_descendant_taxids(childnode.taxid))
            else:
                taxids.append(childnode.taxid)
        return(taxids)

    def _get_descendant_taxids(self,taxid):
        '''inits recursive call for descendants in ToL'''
        taxids = self._recursive_descendant_taxids(taxid)
        taxids.append(taxid)
        return(taxids)

    def _find_taxon(self,name):
        name=name.upper()
        dat=self.db.get_rows_bycolvalue(table_name='tax',valuecolumn='search',
                                        values=name,columns=['path','taxid','rank','search'])
        return(dat)

    def _tax_lookup(self,valuecolumn,values,retrieve=None):
        if retrieve==None:
            retrieve = [el[0] for el in self.tax_columns]

        dat=self.db.get_rows_bycolvalue(table_name='tax',valuecolumn=valuecolumn,
                                        values=values,columns=retrieve)
        return(dat)
    def _taxmap_lookup(self,valuecolumn,values,retrieve=None):

        if retrieve==None:
            retrieve = [el[0] for el in self.taxmap_columns]

        dat=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn=valuecolumn,
                                        values=values,columns=retrieve)
        return(dat)

    def __getitem__(self,key):
            '''
            returns a entry of all assoicated metadata
            '''
            if isinstance(key,str):
                keys = [key]
            else:
                keys = list(key)
            dat = self._taxmap_lookup('accession',keys,retrieve=['accession','path','organism_name'])
            if len(dat) == 1:
                return(dat[0])
            return(dat)

    def __repr__(self):
        rprt = ["Local Silva database",'-'*20]
        rprt.append("  Path      {}".format(self.db_file))
        rprt.append("  Release   {}".format(self.release))
        rprt.append("  Dataset   {}".format(self.dataset))
        rprt.append("  Subunit   {}".format(self.subunit))
        return('\n'.join(rprt))


    #################Exposed to user################# #REQUIRES REVIEW

    def find_taxpath(self,taxon_name):
        '''Search for taxonomic unit in Silva database

        Finds information on the given taxonomic name provided. Spelling must be exact to taxonomic name
        to find results.

        Parameters
        ----------
        taxon_name: str
            name of taxon to be searched. Case insensitive

        Returns
        -------
        tuple
            returns a tuple containing the db entry and the number of sequences belonging to that taxon
        '''
        dat=[]
        if isinstance(taxon_name,str):
            dat = self._tax_lookup('search',taxon_name.upper(),['path','taxid','rank'])
        else:
            search = []
            for el in taxon_name:
                search.append(el.upper())
            dat = self._tax_lookup('search',search,['path','taxid','rank'])

        out = []
        while dat:
            entry = dat.pop()
            out.append((entry,len(self.get_accessions(entry['path']))))
        return(out)
            

    def find_organism(self,name):
        '''Search for an organism in the Silva database

        You can directly search for a organism by the name. Be aware that only exact matches to the organism name will be shown.

        Parameters
        ----------
        name: str
            name of the organism. Case sensitive. 

        Returns
        -------
        list
            List of silva entries
        '''
        dat = self._taxmap_lookup('organism_name',name,['organism_name','accession','path'])
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
        taxid = self._tax_lookup('path',taxpath,retrieve=['taxid'])
        if taxid:
            taxids = self._get_descendant_taxids(taxid[0]['taxid'])
            accessions = self._taxmap_lookup('taxid',taxids,retrieve=['accession'])
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
        dat = self.db.get_rows_bycolvalue(table_name='seqs',valuecolumn='accession',
                                            values=accessions,columns=['accession','seq'])
        return([tuple([el['accession'],el['seq']]) for el in dat])

    def get_aligned(self,accessions):
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
        if isinstance(str,accessions):
            a = set([accessions])
        else:
            a = set(accessions)
        aligned_seqs = []
        for header,sequence in gzip_fasta(self.msa_file):
            dat = header[1:].split(' ') #remove '>' and split by ' '
            if dat[0] in a:
                entry=[dat[0]]
                #if include_path:
                #    path = ' '.join(dat[1:]) #rebuild path
                #    entry.append(path)
                entry.append(sequence.replace('U','T'))
                aligned_seqs.append(entry)

    def write_aligned(self,accessions,file):
        '''Retrieve sequences from input accessions

        For a given list of accessions, retrieve documented sequences

        Parameters
        ----------
        accessions: list of accessions
            list of accessions used for sequence retrieval
        file: str
            file path to save sequences

        '''
        if isinstance(accessions,str):
            a = set([accessions])
        else:
            a = set(accessions)
        total = self.db.total_rows('seqs')
        
        j=0
        j_=0
        print("Decompressing and retrieving {} sequences from {} total sequences".format(len(accessions),total))
        with open(file,'w') as f:
            for header,sequence in gzip_fasta(self.msa_file):
                j+=1
                accession=header[1:].split(' ')[0]
                print("  Searching {0:.2f}%\t".format(j/total*100),end='\r')
                if accession in a:
                    j_+=1
                    f.write(''.join(['>',accession,'\n']))
                    sequence = sequence.replace('U','T')
                    for i in range(0,len(sequence),100):
                        if i+100 <= len(sequence):
                            f.write(sequence[i:i+100]+'\n')
                        else:
                            f.write(sequence[i:len(sequence)]+'\n')
        print("{} sequences writen to {}".format(j_,file))