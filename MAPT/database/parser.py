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

def compress_gaps(accession,seq):
    ''' this function produces an encoding of seqs without gaps'''
    encoding = False
    current_encoding = [str(),str(),str()]
    compression = []
    i=0 #used to count all non - or .
    e=int() #used to count sequence of - or .
    for char in seq:
        #only count when not . or -
        if encoding:
            if not (char == '.' or char == '-'):
                current_encoding[2] = str(e) #insert e of char
                compression.append(','.join(current_encoding))
                encoding=False
                current_encoding = [str(),str(),str()]
                i+=1
            else:
                e+=1
        else:
            if char == '.' or char == '-':
                encoding = True
                e = 1 #number of char to encode
                current_encoding[0] = char #what character
                current_encoding[1] = str(i) #start position relative to compression
            else: 
                i+=1
    if current_encoding[0]:
        current_encoding[2] = str(e) #end position
        compression.append(','.join(current_encoding))
    return(accession,';'.join(compression))

def inflate(encoding):
    char,pos,numchar = encoding.split(',')
    return(char*int(numchar),int(pos))

def inflate_gaps(accession,seq,compression):
    gapseq=[]
    last_encoding_position=int()
    for encoding in compression.split(';'):
        gaps,pos = inflate(encoding)
        if pos: #add nt from gapless sequence
            #print([last_encoding_position,pos])
            gapseq.append(seq[last_encoding_position:pos])
        gapseq.append(gaps)
        last_encoding_position = pos
    return(accession,''.join(gapseq))

def read_fasta(fasta_file):
    accession_seq = {}                                                                                                        
    with open(fasta_file, 'r') as f:
        seq=""
        header=""
        for line in f:
            if line[0] == '>':
                if header!="":
                    accession = header[1:].split(' ')[0]
                    accession_seq[accession]=seq
                header=line.strip()
                seq=""
            else:
                temp=line.strip()#removing end of line
                temp=temp.replace(" ","")#removing any possible spaces
                seq+=temp
        #for the last sequence
        if header:
            accession = header[1:].split(' ')[0]
            accession_seq[accession]=seq
    return(accession_seq)

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

def gzip_fasta_batch(fasta_file,batchsize=5000):                                                                                                        
    with gzip.open(fasta_file, 'r') as f:
        seq=""
        header=""
        dat=[]
        seq = []
        for line in f:
            line = str(line,'utf-8')
            if line[0] == '>':
                if header!="":
                    dat.append((header,''.join(seq)))
                header=line.strip()
                seq=[]
            else:
                temp=line.strip()#removing end of line
                temp=temp.replace(" ","")#removing any possible spaces
                seq.append(temp)
            if len(dat) == batchsize:
                yield dat
                dat = []
        #for the last sequence
        yield dat

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