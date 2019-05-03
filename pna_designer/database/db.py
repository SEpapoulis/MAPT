from .sql_core import sqldb
from .parser import gzip_fasta,gzip_table
from .TreeOfLife import load_ToL

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


#TODO: taxmap and embl are redundant!!!! merge the two into one spreadsheet
#       Ideally, put ncbi_taxid as a column in taxmap

class silva_db():
    def __init__(self,dbfile,tree_file=str(),fasta_file=str(),
    taxmap_file=str(),tax_file=str(),embl_file=str()):

        self.db = sqldb(dbfile)

        #Tree of life is a dictionary of nodes
        #nodes have references to parent and child(ren) nodes

        self.tables=('seqs','taxmap','tax','embl')
        #NOTE: accessions also have genomic location, true accessions
        #      are accessions.split('.')[0], where locations are [1] and [2]
        
        #accession is primary key
        self.seqs_columns = (('accession', 'TEXT'),
                                ('path', 'TEXT'),
                                ('seq','TEXT'))
        
        #accession is primary key
        self.taxmap_columns = (('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('taxid','INTEGER'))

        #accession is primary key
        #NOTE:ncbi_taxid is moved to taxmap and this table is dropped from the db
        self.embl_columns = (('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('ncbi_taxid','INTEGER'))
        #taxid primary key
        self.tax_columns = (('path', 'TEXT'), ('taxid', 'INTEGER'), 
                                ('rank','TEXT'), ('remark', 'TEXT'), 
                                ('release','INTEGER'))

        #Adding data
        if tree_file:
            self.ToL=load_ToL(tree_file)
        if fasta_file:
            self.add_fasta(fasta_file)
        if taxmap_file:
            self.add_taxmap(taxmap_file)
        if tax_file:
            self.add_tax(tax_file)
        if embl_file:
            self.add_embl(embl_file)
            if self.db.table_exists('taxmap'):
                self.db.copycol_bypk('taxmap','embl','accession','ncbi_taxid')
                self.db.drop_table('embl')
        #used as a error check for inputs
        if self.db.table_exists('tax'):
            self._ranks = set([el[0] for el in self.db.get_columns('tax','rank')])     

    def _col_names(self,table_columns):
        return([el[0] for el in table_columns])

    #this funciton will add our fasta data into the sql database
    #700000 sequences at a time (as to not overflow RAM)
    def add_fasta(self,fasta_file,write_size=700000):
        self.db.drop_table('seqs')
        self.db.create_table('seqs',columns=self.seqs_columns,primary_key=['accession'])
        seqs_col=self._col_names(self.seqs_columns)
        print('Adding sequences to database....', end='\r')
        i = 1
        to_db=[]
        for header,sequence in gzip_fasta(fasta_file):
            dat = header[1:].split(' ') #remove '>' and split by ' '
            path = ' '.join(dat[1:]) #rebuild path
            dat=[dat[0],path,sequence.replace('U','T')]
            to_db.append(dat)
            if i % write_size == 0:
                self.db.insert_rows('seqs',seqs_col,to_db)
                to_db = []
            i+=1
        #adding remaining sequences
        if i % write_size != 0:
            self.db.insert_rows('seqs',seqs_col,to_db)
        
        print('Adding sequences to database....  [DONE]')


    def add_taxmap(self,taxmap_file):
        self.db.drop_table('taxmap')
        self.db.create_table('taxmap',columns=self.taxmap_columns,
                            primary_key=['accession'])
        taxmap_cols=self._col_names(self.taxmap_columns)
        to_db=[]
        for line in gzip_table(taxmap_file):
            accession = '.'.join(line[0:3])
            to_db.append([accession]+line[3:])
        self.db.insert_rows('taxmap',taxmap_cols,to_db)
        
    
    def add_tax(self,tax_file):
        self.db.drop_table('tax')
        self.db.create_table('tax',columns=self.tax_columns,primary_key=['taxid'])
        tax_cols=self._col_names(self.tax_columns)
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
        self.db.drop_table('embl')
        self.db.create_table('embl',columns=self.embl_columns,
                            primary_key=['accession'])
        embl_cols=self._col_names(self.embl_columns)
        to_db=[]
        for line in gzip_table(embl_file):
            accession = '.'.join(line[0:3])
            to_db.append([accession]+line[3:])
        self.db.insert_rows('embl',embl_cols,to_db)
        
    def get_seqs(self,accessions):
        seqs=self.db.get_rows_bycolvalue(table_name='seqs',valuecolumn='accession',
                                        values=accessions,columns=['accession','seq'])
        return(seqs)


    def get_organism_byname(self,name):
        name=name.lower()
        dat=self.db.get_columns(table_name='taxmap',columns=['organism_name','accession','path','taxid'])
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
    def get_descendant_taxids(self,taxid):
        taxids=[]
        node = self.ToL[taxid]
        for childnode in node.children:
            if childnode.children:
                taxids.extend(self.get_descendant_taxids(childnode.taxid))
            else:
                taxids.append(childnode.taxid)
        return(taxids)
    
    def get_accessions(self,taxid):
        accessions=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn='taxid',
                                        values=taxid,columns='accession')
        return(accessions)

    def get_accession_dat(self,accession,get_columns=None):
        if get_columns==None:
            get_columns = [el[0] for el in self.taxmap_columns]

        accessions=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn='accession',
                                        values=accession)
        return(accessions)

    def get_sequence_dat(self,accession,get_columns=None):
        if get_columns==None:
            get_columns = [el[0] for el in self.seqs_columns]
        seq_dat=self.db.get_rows_bycolvalue(table_name='seqs',valuecolumn='accession',
                                        values=accession,columns=get_columns)
        return(seq_dat)
    
    def write_fasta(self,accessions,fasta_file,include_taxid=False,include_path=False,
                   include_organism_name=False,include_NCBItaxid=False):
        taxmap_metadata=[]
        if include_path:
            taxmap_metadata.append('path')
        if include_taxid:
            taxmap_metadata.append('taxid')
        if include_organism_name:
            taxmap_metadata.append('organism_name')
        if taxmap_metadata:
            metadata=self.db.get_rows_bycolvalue(table_name='embl',valuecolumn='accession',
                                            values='accession',columns=taxmap_metadata)
        
        if include_NCBItaxid:
            temp=self.db.get_rows_bycolvalue(table_name='embl',valuecolumn='accession',
                                        values='accession',columns=['accession,ncbi_taxid'])
            accession_NCBI = {el['accession']:el['ncbi_taxid'] for el in temp}
        else:
            accession_NCBI={}

        


