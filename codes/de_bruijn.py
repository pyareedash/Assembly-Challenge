import numpy as np
from correct_errors import *
from utils import *


class Node:
    def __init__(self, label, coverage):
        self.label = label
        self.coverage = coverage
        self.incoming = set()
        self.outgoing = set()
    
    def add_outgoing_edge(self, node_in):
        self.outgoing.add(node_in)
        node_in.incoming.add(self)

class DeBruijnGraph:
    def __init__(self,k,k_mer_dict, alphabet = ["A","C","G","T"]):
        self.k = k
        self.collapse_nodes(k_mer_dict,alphabet)
    
    def collapse_nodes(self,k_mer_dict, alphabet = ["A","C","G","T"]):
        #start points
        starts = set()
        for k_mer in k_mer_dict:
            incoming, outgoing = query_de_bruijn(k_mer_dict,k_mer,alphabet = alphabet)
            if len(incoming) != 1:
                starts.add(k_mer)
            if len(outgoing) > 1:
                for out_node in outgoing:
                    starts.add(out_node)
        self.nodes = []
        #collapse nodes
        for node in starts:
            new_node = node
            node_count = k_mer_dict[node]
            incoming, outgoing = query_de_bruijn(k_mer_dict,node, alphabet = alphabet)
            while len(outgoing) == 1:
                next_node = outgoing[0]
                next_incoming, next_outgoing = query_de_bruijn(k_mer_dict,next_node,alphabet = alphabet)
                if len(next_incoming) == 1:
                    new_node = new_node + next_node[-1]
                    node_count += k_mer_dict[next_node]
                    outgoing = next_outgoing
                else:
                    break
            node_avg = node_count / (len(new_node)-self.k+1)
            new_node = Node(new_node,node_avg)
            self.nodes.append(new_node)        
        #edges in graph
        self.nodes = set(self.nodes)
        self.create_edges()
        return
        
    
    def create_edges(self):
        #print([node.label for node in self.nodes])
        starts = {}
        ends = {}
        for node in self.nodes:
            start = node.label[:self.k-1]
            end = node.label[-self.k+1:]
            for key, dct in [[start,starts], [end,ends]]:
                if key in dct:
                    dct[key].append(node)
                else:
                    dct[key] = [node]
        self.edges = []
        for key in starts:
            if key in ends:
                for node_out in ends[key]:
                    for node_in in starts[key]:
                        node_out.add_outgoing_edge(node_in)
                        self.edges.append((node_out,node_in)) 
        self.edges = set(self.edges)
    
    def remove_node(self,node):
        #print('removing',node.label)
        self.nodes.remove(node)
        for in_node in node.incoming:
            #if node not in in_node.outgoing:
            #    print(in_node.label,in_node,'removing',node.label,node)
            in_node.outgoing.remove(node)
            self.edges.remove((in_node, node))
        for out_node in node.outgoing:
            #if node not in out_node.incoming:
            #    print(node.label,out_node.label)
            out_node.incoming.remove(node)
            self.edges.remove((node, out_node))
        changed_nodes = node.incoming.union(node.outgoing)
        return changed_nodes
      
    
    def simplify(self, changed_nodes):
        contract_in = [(list(node.incoming)[0],node)
                       for node in changed_nodes 
                           if len(node.incoming)==1 
                               and len(list(node.incoming)[0].outgoing)==1
                      ]
        contract_out = [(node,list(node.outgoing)[0]) 
                        for node in changed_nodes 
                            if len(node.outgoing)==1
                                and len(list(node.outgoing)[0].incoming)==1
                       ]
        contract_edges = list(set(contract_in + contract_out))
        deleted_nodes = set()
        new_nodes = set()
        while len(contract_edges)>0:
            (out_node,in_node) = contract_edges.pop(0)
            if (out_node in deleted_nodes or in_node in deleted_nodes):
                continue
            if (out_node.label == in_node.label):
                #print("equal")
                continue
            new_node = self.merge(out_node,in_node)
            new_nodes.add(new_node)
            deleted_nodes.add(out_node)
            deleted_nodes.add(in_node)
            #self.changed_nodes += [out_node,in_node]
            #add new edges
            if len(new_node.incoming) ==1 and len(list(new_node.incoming)[0].outgoing) == 1:
                contract_edges.append((list(new_node.incoming)[0],new_node))
            if len(new_node.outgoing) ==1 and len(list(new_node.outgoing)[0].incoming) == 1:
                contract_edges.append((new_node,list(new_node.outgoing)[0]))
        return new_nodes - deleted_nodes
    
    
    def merge(self, out_node,in_node):
        new_label = out_node.label + in_node.label[self.k-1:]
        w_out = len(out_node.label)-self.k+1
        w_in = len(in_node.label)-self.k+1
        new_coverage = (out_node.coverage*w_out + in_node.coverage*w_in)/(w_out+w_in)
        new_node = Node(new_label,new_coverage)
        new_node.incoming = out_node.incoming.copy()
        new_node.outgoing = in_node.outgoing.copy()
        #print('merge\n',out_node.label,'\n',in_node.label)
        if new_label[:self.k-1] == new_label[-self.k+1:]:
            new_node.incoming.add(new_node)
            new_node.outgoing.add(new_node)
            out_node.outgoing.add(new_node)
            in_node.incoming.add(new_node)
            new_node.incoming.add(in_node)
            new_node.outgoing.add(out_node)
            new_node.incoming.add(out_node)
            new_node.outgoing.add(in_node)
            #print("equal after merge")
            #print(out_node.label,in_node.label)
        for new_in in new_node.incoming:
            new_in.outgoing.add(new_node)
            self.edges.add((new_in, new_node))
        for new_out in new_node.outgoing:
            new_out.incoming.add(new_node)
            self.edges.add((new_node, new_out))
        #new_node.before_merge = [out_node,in_node]
        #print('new',new_node.label)
        self.nodes.add(new_node)
        self.remove_node(out_node)
        self.remove_node(in_node)
        return new_node
    
    def collapse_bubble(self,node1,node2):
        node2.incoming = node1.incoming.union(node2.incoming)
        node2.outgoing = node1.outgoing.union(node2.outgoing)
        for in_node in node1.incoming:
            in_node.outgoing.add(node2)
        for out_node in node1.outgoing:
            out_node.incoming.add(node2)
        node2.coverage = node1.coverage+node2.coverage
        changed_nodes = self.remove_node(node1)
        new_nodes = self.simplify(changed_nodes)
        return new_nodes
