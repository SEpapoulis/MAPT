from .sql_core import sqldb
from .parser import gzip_fasta,gzip_table
from .TreeOfLife import load_ToL

class silva_db():
    def __init__(self,dbfile,tree_file=str()):

        self.db = sqldb(dbfile)

        #Tree of life is a dictionary of nodes
        #nodes have references to parent and child(ren) nodes
        if tree_file:
            self.ToL=load_ToL(tree_file)


        self.tables=['seqs','taxmap','tax','embl']

        #accession is primary key
        self.seqs_columns = [('accession', 'TEXT'),
                                ('path', 'TEXT'),
                                ('seq','TEXT')]
        
        #accession is primary key
        self.taxmap_columns = [('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('taxid','INTEGER')]

        #accession is primary key
        self.embl_columns = [('accession', 'TEXT'), 
                                ('path', 'TEXT'),
                                ('organism_name', 'TEXT'), 
                                ('ncbi_taxid','INTEGER')]
        #taxid primary key
        self.tax_columns = [('path', 'TEXT'), ('taxid', 'INTEGER'), 
                                ('rank','TEXT'), ('remark', 'TEXT'), 
                                ('release','INTEGER')]

        #used as a error check for inputs
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

    def get_name(self,identifier):
        at=self.db.get_rows_bycolvalue(table_name='taxmap',valuecolumn='taxid',
                                        values=rank.lower(),columns=['path','taxid'])

    #Using rank to limit searchspace, we look for the name among 
    def get_taxon_byname(self,name,rank):
        if rank not in self._ranks:
            raise ValueError("{} is not a phylogenetic rank, please specifyone of the following:\n{}".format(
                rank,'\n'.join(self._ranks)))

        dat=self.db.get_rows_bycolvalue(table_name='tax',valuecolumn='rank',
                                        values=rank.lower(),columns=['path','taxid'])
        names = {el[0].split(';')[-2].lower():el[1] for el in dat}
        if name.lower() in names.keys():
            taxid=names[name.lower()]
            print("{} taxid = {}".format(name,str(taxid)))
            return(taxid)
        else:
            print("{} was not found in the taxanomic rank of {}".format(name,rank))
        

    def get_descendant_accessions(self,taxid):
        taxids=[]
        node = self.ToL[taxid]
        for childnode in node.children:
            if childnode.children:
                taxids.extend(self.get_descendant_accessions(childnode.taxid))
            else:
                taxids.append(childnode.taxid)
        return(taxids)

