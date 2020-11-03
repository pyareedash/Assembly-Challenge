import numpy as np
from de_bruijn import *
from correct_errors import *


def assemble_data(input_file, k, k_mer, output_file):
    #TODO read text file
    reads = load_reads(input_file)
    #TODO construct k_mer_dict
    k_mer_dict = count_k_mers(reads, k)
    #TODO query the graph
    query_de_bruijn(k_mer_dict, k_mer, prnt = True)
    
    #TODO extra: correct error in reads before collapsing the nodes
    #TODO collapse nodes of de bruijn
    graph = DeBruijnGraph(k,k_mer_dict)
    #print('Number of nodes: ',len(graph.nodes))
    #print('Number of edges: ',len(graph.edges))
    #TODO correct errors
    #print('remove low coverage')
    remove_low_coverage(graph,2)
    #print('Number of nodes: ',len(graph.nodes))
    #print('Number of edges: ',len(graph.edges))
    #print('remove tips')
    remove_tips(graph,k)
    #print([n.label for n in graph.nodes])
    #print(len([n.label for n in graph.nodes]))
    #print('Number of nodes: ',len(graph.nodes))
    #print('Number of edges: ',len(graph.edges))
    #print('remove bubbles')
    remove_bubbles(graph)
    #print('Number of nodes: ',len(graph.nodes))
    #print('Number of edges: ',len(graph.edges))
    
    write_fasta(graph, output_file)
    
    
    return graph
