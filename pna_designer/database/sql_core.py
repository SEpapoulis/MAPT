import sqlite3


class sqldb():
    def __init__(self,db_file=str()):
        self.conn,self.c = self.connect(db_file)
        

    def connect(self,db):
        conn = sqlite3.connect(db)
        #this allows row elements to be called by columns
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        return(conn,c)

    def close(self):
        self.conn.close()

    def create_table(self,table_name,primary_key=list(),columns=list()):
        column_field=[]
        for field in columns:
            column_field.append(' '.join(field))
        if primary_key:
            pk = ', PRIMARY KEY('+','.join(primary_key)+')'
        else:
            pk=''
        self.c.execute('CREATE TABLE IF NOT EXISTS {tn} ({fields}{p})'.format(
            tn=table_name, fields=','.join(column_field),p=pk))
        self.conn.commit()

    def drop_table(self,table_name):
        self.c.execute('DROP TABLE IF EXISTS {tn}'.format(tn=table_name))
        self.conn.commit()

    def total_rows(self,table_name, print_out=False):
        """ Returns the total number of rows in the database """
        self.c.execute('SELECT COUNT(*) FROM {}'.format(table_name))
        count = self.c.fetchall()
        if print_out:
            print('\nTotal rows: {}'.format(count[0][0]))
        return (count[0][0])

    def table_col_info(self, table_name, print_out=False):
        """ Returns a list of tuples with column informations:
        (id, name, type, notnull, default_value, primary_key)
        """
        self.c.execute('PRAGMA TABLE_INFO({})'.format(table_name))
        info = self.c.fetchall()

        if print_out:
            print("\nColumn Info:\nID, Name, Type, NotNull, DefaultVal, PrimaryKey")
            for col in info:
                print(col)
        return info


    def insert_rows(self,table_name,column_names,rows):
        qm = ','.join(['?']*len(column_names))
        self.c.executemany("INSERT INTO {tm} ({c}) VALUES ({q})".format(
            tm=table_name,c=','.join(column_names),q=qm),rows)
        self.conn.commit()

    def get_columns(self,table_name,columns='*'):
        if type(columns) == type(list()):
            col = ", ".join(columns)
        else:
            col = columns
        
        self.c.execute('SELECT {c} FROM {t} '.format(c=col, t=table_name))
        return(self.c.fetchall())
    

    #get values from specified columns. all columns are 
    #returned if nothing is passed
    def get_rows_bycolvalue(self,table_name,valuecolumn,values,columns=['*']):
        if type(columns) == str:
            columns=[columns]
        if type(values) == type(list()):
            values = [str(el) for el in values]
            vals = "'"+"', '".join(values)+"'"
        else:
            vals = "'"+str(values)+"'"
        self.c.execute("SELECT {c} FROM {tn} WHERE {vc} IN ({v})".format(
            tn=table_name,c=','.join(columns),vc=valuecolumn,v = vals)
        )
        return(self.c.fetchall())
    
    def table_exists(self,tablename):
        self.c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='{tn}'".format(tn=tablename))
        fetchlist=self.c.fetchall()
        if fetchlist:
            return(True)
        else:
            return(False)

    def itertable(self,tablename):
        self.c.execute('SELECT * FROM {t}'.format(t=tablename))
        for row in self.c:
            yield(row)
    
    def add_column(self,table,column):
        col,coltype=column
        self.c.execute("ALTER TABLE {t} ADD COLUMN {c} {ct}".format(t=table,c=col,ct=coltype))

    def copycol_bypk(self,table,tablesource,keycol,datcol):
        cols=self.table_col_info(tablesource)
        cols={el['name']:el for el in cols}
        if keycol not in cols:
            raise Exception("keycol not found in {}".format(tablesource))
        if datcol not in cols:
            raise Exception("datcol not found in {}".format(tablesource))

        #creating a new column in table named datcol with datcol type from tablesource
        self.add_column(table,(datcol,cols[datcol]['type']))
        
        #populate column where keycol table == keycol tablesource
        self.c.execute("UPDATE {t} SET {dc} = (SELECT {dc} FROM {ts} WHERE {kc} = {t}.{kc})".format(
            t= table , ts = tablesource, kc=keycol, dc=datcol
        ))
        self.conn.commit()

    #This is a work in progress, still unclear how to use extentions
    #def spelling(self,table,column):
    #    self.c.execute("CREATE VIRTUAL TABLE demo USING spellfix1")
    #    self.c.execute("INSERT INTO demo(word) SELECT {} FROM {};".format(column,table))
    #    self.conn.commit()
    
    
    
        '''this is a sample funciton from a tutorial
    def values_in_col(self, table_name, print_out=True):
        """ Returns a dictionary with columns as keys
        and the number of not-null entries as associated values.
        """
        self.c.execute('PRAGMA TABLE_INFO({})'.format(table_name))
        info = self.c.fetchall()
        col_dict = dict()
        for col in info:
            col_dict[col[1]] = 0
        for col in col_dict:
            self.c.execute('SELECT ({0}) FROM {1} '
                    'WHERE {0} IS NOT NULL'.format(col, table_name))
            # In my case this approach resulted in a
            # better performance than using COUNT
            number_rows = len(self.c.fetchall())
            col_dict[col] = number_rows
        if print_out:
            print("\nNumber of entries per column:")
            for i in col_dict.items():
                print('{}: {}'.format(i[0], i[1]))
        return col_dict
    '''