import os
import gzip


amb_baseset={'R':set(['A','G']),
            'Y':set(['C','T']),
            'S':set(['G','C']),
            'W':set(['A','T']),
            'K':set(['G','T']),
            'M':set(['A','C']),
            'B':set(['C','T','G']),
            'D':set(['A','G','T']),
            'H':set(['A','C','T']),
            'V':set(['A','C','G']),
            'N':set(['A','C','T','G']),
            }

nt_complement={
    'A':'T','T':'A','G':'C','C':'G', 
    'N':'N','R':'Y','Y':'R','S':'S',
    'W':'W','K':'M','M':'K','B':'V',
    'V':'B','D':'H','H':'D'}

def chunkify(fname,size=1024*1024*100):
    fileEnd = os.path.getsize(fname)
    chunkEnd=0
    chunkStart=0
    with open(fname,'rb') as f:
        while chunkStart < fileEnd:
            chunkStart = chunkEnd
            if chunkStart + size > fileEnd:
                size = fileEnd -chunkStart
            f.seek(size,1)
            f.readline() #gets pointer to the end of a line!
            chunkEnd = f.tell()#endchunk is now the end of a line
            yield (chunkStart, chunkEnd - chunkStart)

def gzip_table(filename,skip_first_row=True):
    with gzip.open(filename,'r') as f:
        if skip_first_row:
            f.readline()
        for line in f:
            dat = str(line,'utf-8').strip()
            dat = dat.split('\t')
            yield (dat)

def gzip_fasta(fasta_file):
    with gzip.open(fasta_file, 'r') as f:
        seq=""
        header=""
        for line in f:
            line = str(line,'utf-8')
            if line[0] == '>':
                if header!="":
                    yield (header,seq)
                header=line.strip()
                seq=""
            else:
                temp=line.strip()#removing end of line
                temp=temp.replace(" ","")#removing any possible spaces
                seq+=temp
        #for the last sequence
        if header:
            yield(header,seq)


#returns a list of tuples, where each tuple is (decendants,taxid)
def split_newick(s):
    counter = 0
    nodedat = tuple()
    descendants = []
    taxid = str()
    for i,char in enumerate(s):
        if char == '(':
            if counter == 0:
                start = i+1
            counter +=1
        elif char ==')':
            counter-=1
            if counter == 0:
                stop=i
                nodedat = (start,stop)
        elif counter < 1 and char.isdigit():
            taxid = taxid+char
        elif counter < 1 and char ==',':
            if taxid:
                taxid=int(taxid)
            else:
                taxid=None
            if nodedat:
                descendants.append((s[nodedat[0]:nodedat[1]],taxid))
            else:
                descendants.append((None,taxid))
            taxid = str()
            nodedat = tuple()
    #add before exit
    if taxid:
        taxid=int(taxid)
    else:
        taxid=None
    if nodedat:
        descendants.append((s[nodedat[0]:nodedat[1]],taxid))
    else:
        descendants.append((None,taxid))
    return(descendants)