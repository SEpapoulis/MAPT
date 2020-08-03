
nt_complement={
        'A':'T','T':'A','G':'C','C':'G', 
        'N':'N','R':'Y','Y':'R','S':'S',
        'W':'W','K':'M','M':'K','B':'V',
        'V':'B','D':'H','H':'D'}

def complement(DNA):
    DNA=DNA.upper()
    complement = [nt_complement[nt] for nt in DNA]
    return(''.join(complement))
    
def reverse_complement(DNA):
    DNA=DNA.upper()
    complement = [nt_complement[nt] for nt in DNA]
    return(''.join(complement[::-1]))


def get_primerSS(primer):
    amb={'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],
        'B':['C','G','T'],'V':['A','C','G'],'D':['A','G','T'],'H':['A','C','T'],'N':['A','T','G','C']}
    Pset=['']
    for char in primer:
        if char in amb: #for each sequence in Fset, duplicate and add each possilbe amb to child
            nts = amb[char]
            _temp=[]
            for seq in Pset:
                for nt in nts:
                    _temp.append(seq+nt)
            Pset=_temp
        else:
            for i in range(0,len(Pset)):
                Pset[i]=Pset[i]+char
    return(set(Pset))

def search_primer(seq,primer,direction):
    primerlen=len(primer)
    primer_set=get_primerSS(primer)
    for i in range(primerlen,len(seq)):
        end=i
        start=i-primerlen
        window=seq[start:end]
        if window in primer_set:
            if direction == 'F':
                return(start)
            else:
                return(end+1)



class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class PrimerError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message    
    """
    pass

def sim_amplify(Forward_primer,Reverse_primer,seq):
    #for searchers, we need the reverse complement when scanning DNA
    Reverse_primer_RC = reverse_complement(Reverse_primer)
    start=search_primer(seq,Forward_primer,'F')
    end=search_primer(seq,Reverse_primer_RC,'R')
    if start and end:
        if end<start:
            end=search_primer(seq[start:],Reverse_primer_RC,'R')+start
    #retry with reverse complement
    if start == None or end == None:
        start=search_primer(reverse_complement(seq),Forward_primer,'F')
        end=search_primer(reverse_complement(seq),Reverse_primer_RC,'R')
        if start and end:
            if end<start:
                end=search_primer(reverse_complement(seq)[start:],Reverse_primer_RC,'R')+start
    if start == None or end == None:
        raise PrimerError("Could not find forward or reverse primer in sequence")
    return(seq[start:end])

