from .parser import split_newick
import gzip

#Base datastructre for tree
class Node:
    def __init__(self,taxid):
        self.taxid = taxid
        self.parent = None
        self.children = []

    def get_taxid(self):
        return(self.taxid)

    def set_name(self,name):
        self.name = name

    def get_children(self):
        return(self.children)

    def add_child(self,child):
        self.children.append(child)
    
    def get_parent(self):
        return(self.parent)
    def set_parent(self,parent):
        self.parent = parent

#TreeOfLife is a list of Nodes, taxid is an int, children is [(raw_newick,taxid) ...], 
# and parent is a node
def add_node(taxid_node,taxid,children,parent_node):
    #create current node and add it as a child to parent
    current_node = Node(taxid)
    if parent_node != None:
        parent_node.add_child(current_node)
    #add current node to dicitonary for automatic lookup
    taxid_node[taxid] = current_node
    #make sure current node can reference the parent
    current_node.set_parent(parent_node)
    #we will now check all children to see if they are leaf
    for child in children:
        cnode,childtax = child
        #cnode == None when there is no more newick to parse
        if cnode != None:
            #more newick to parse, add node recursivly 
            add_node(taxid_node,childtax,split_newick(cnode),current_node)
        else:
            child_node = Node(childtax)
            child_node.set_parent(current_node)
            taxid_node[childtax] = child_node
            current_node.add_child(child_node)

def load_ToL(file):
    taxid_node={}
    if '.gz' in file:
        f =gzip.open(file,'r')
        dat = f.readlines()[0]
        dat = str(dat,'utf-8')
    else:
        with open(file, 'r') as f:
            dat=f.readlines()[0]
    dat=dat.strip(';\n')
    newick_raw,roottax=split_newick(dat)[0]
    children = split_newick(newick_raw)
    add_node(taxid_node,roottax,children,None)
    return(taxid_node)        
