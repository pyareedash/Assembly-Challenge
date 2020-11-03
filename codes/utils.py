import numpy as np
#from de_bruijn import *
#from correct_errors import *

def load_reads(filename):
    f = open(filename, 'r')
    reads = []
    for line in f:
        if not line.startswith('>'):
            reads.append(line.strip())
    return reads

def reverse_complement(read):
    rc = {"A":"T","T":"A", "C":"G", "G":"C"}
    read_rc = [rc[c] for c in reversed(read)]
    read_rc = ''.join(read_rc)
    return read_rc

def count_k_mers(reads, k):
    #TODO it shall be a dictionary
    k_mer_dict = {}
    reads_rc = [reverse_complement(r) for r in reads]
    for r in reads+reads_rc:
        for i in range(len(r)-k+1):
            v = r[i:i+k]
            if v in k_mer_dict:
                k_mer_dict[v]+=1
            else:
                k_mer_dict[v]=1
    return k_mer_dict

def print_list(l):
    l_str = ''
    for x in l:
        l_str += x + ', '
    return l_str[:-2]

def query_de_bruijn(k_mer_dict, k_mer, alphabet = ["A","C","G","T"], prnt = False):
    if k_mer not in k_mer_dict:
        print("NOT IN GRAPH")
        return
    incoming = []
    outgoing = []
    for letter in alphabet:
        in_node = letter + k_mer[:-1]
        out_node = k_mer[1:] +letter
        if in_node in k_mer_dict:
            incoming.append(in_node)
        if out_node in k_mer_dict:
            outgoing.append(out_node)
    if prnt:
        print("Incoming k-mers - "+print_list(incoming))
        print("Outgoing k-mers - "+print_list(outgoing))
        
    #TODO
    return [incoming, outgoing]

def write_fasta(de_bruijn,output_name):
    f = open(output_name, 'w')
    for i,node in enumerate(de_bruijn.nodes):
        f.write('>%i\n'%(i+1))
        f.write(node.label+'\n')
    f.close()
