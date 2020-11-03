import numpy as np
#from utils import *
#from de_bruijn import *

def remove_low_coverage(graph, below):
    nodes_to_remove = [node for node in graph.nodes if node.coverage<below]
    l = len(nodes_to_remove)
    #print('nodes to remove: ',l)
    #i=0
    for node in nodes_to_remove:
        if node in graph.nodes:
            changed_nodes = graph.remove_node(node)
            _ = graph.simplify(changed_nodes)
        #i+=1
        #if i%10000 == 0: print('remaining:',l-i)

def remove_tips(graph,k):
    tip_candidates = order_nodes_by_coverage([node for node in graph.nodes if testtip(node,k)])
    #print(len(tip_candidates))
    #i = 0
    while len(tip_candidates)>0:
        #i+=1
        #if i %100 == 0:
            #print('tips remaining: ',len(tip_candidates))
        node = tip_candidates.pop(0)
        if node not in graph.nodes:
            continue
        #print('remove tip', node.coverage)
        changed_nodes = graph.remove_node(node)
        new_nodes = graph.simplify(changed_nodes)
        for new_n in new_nodes:
            if testtip(new_n,k):
                tip_pos = bisect.bisect([n.coverage for n in tip_candidates],new_n.coverage)
                tip_candidates.insert(tip_pos, new_n)
    
def testtip(node,k):
    if len(node.incoming) == 0 or len(node.outgoing) ==0:
        if len(node.label) < 2*k:
            return True
    return False

def remove_bubbles(graph, threshold =0.1):
    bubble_candidates = []
    sorted_nodes = order_nodes_by_coverage(graph.nodes)
    for i,node in enumerate(sorted_nodes):
        bubble_candidates += find_bubble_pairs(node,sorted_nodes[i+1:])
    #print(len(bubble_candidates))
    while len(bubble_candidates)>0:
        (node1,node2)= bubble_candidates.pop(0)
        #print(len(bubble_candidates))
        if (node1 not in graph.nodes) or node2 not in (graph.nodes):
            continue
        #print(node1.label, node2.label)
        if hybrid_similarity(node1.label, node2.label, threshold):#is_similar(node1.label, node2.label,threshold):
            #node2.coverage is always larger or equal than node1.coverage
            new_nodes = graph.collapse_bubble(node1,node2)
            bubbles_to_add = []
            for new_node in new_nodes:
                bubbles_to_add += find_bubble_pairs(new_node,graph.nodes)
            bubble_candidates = bubble_candidates + bubbles_to_add
            bubble_coverages = [(node1.coverage,node2.coverage) for (node1,node2) in bubble_candidates]
            bubble_order = np.argsort(np.array(bubble_coverages,dtype=np.dtype([('x', float), ('y', float)])))
            bubble_candidates = [bubble_candidates[i] for i in bubble_order]

def is_similar(s,t,threshold = 0.1):
    distance = edit_distance(s,t)
    length = len(s)
    if distance/length < threshold:
        return True
    return False
            
            
def find_bubble_pairs(node,other_nodes):
    bubble_pairs = []
    for other in other_nodes:
        if node != other:
            common_in_nodes = node.incoming.intersection(other.incoming)
            common_out_nodes = node.outgoing.intersection(other.outgoing)
            if len(common_in_nodes) != 0 and len(common_out_nodes) != 0:
                bubble_pairs.append(tuple(order_nodes_by_coverage([node, other])))
    return bubble_pairs
    
def order_nodes_by_coverage(nodes):
    l_nodes = list(nodes)
    ordering = np.argsort([node.coverage for node in l_nodes])
    ordered_nodes = [l_nodes[i] for i in ordering]
    return ordered_nodes

def hybrid_similarity(s, t, threshold = 0.1,window = 50):
    ls, lt = len(s), len(t)
    minl = min(ls,lt)
    if abs(ls-lt) > threshold*minl:
        #print('rejecting by length')
        return False
    char_diff = character_differences(s,t)
    if char_diff > 2*threshold*minl:
        #print('rejecting by character differences')
        return False
    avgdiff = compare_random_slices(s,t, window)
    if avgdiff > 2*threshold:
        #print('rejecting by random samples')
        return False
    #print('computing edit distance on random slices')
    avgdiff = compute_edit_distance_on_slices(s,t,window)
    #print('end: computing edit distance')
    if avgdiff < threshold:
        #print('accepting by edit_distance')
        return True
    else:
        #print('rejecting by edit distance')
        return False
    
def compare_random_slices(s,t, window = 50):
    avgdiff = 0
    #window = int(size*max(len(s),len(t)))
    minl = min(len(s),len(t))
    if minl>window:
        start_points = np.random.choice(range(minl-window), int(minl*0.5))
        for idx in start_points:
            diff = character_differences(s[idx:idx+window], t[idx:idx+window])
            avgdiff += diff
        return avgdiff/(len(start_points)*window)
    else:
        return character_differences(s,t)/minl

def character_differences(s,t):
    letters = set(s+t)
    dict_s = {l:0 for l in letters}
    dict_t = {l:0 for l in letters}
    for l in s:
        dict_s[l]+=1
    for l in t:
        dict_t[l]+=1
    diff = 0
    for l in letters:
        diff += abs(dict_s[l]-dict_t[l])
    return diff
    
def compute_edit_distance_on_slices(s,t, window = 50):
    minl = min(len(s),len(t))
    start_points = np.random.choice(range(int(minl)), int(minl*0.01))
    if len(start_points) == 0:
        return edit_distance(s, t)/minl
    avgdiff = 0
    for idx in start_points:
        diff = edit_distance(s[idx:idx+window], t[idx:idx+window])
        avgdiff += diff
    return avgdiff/(len(start_points)*window)

def edit_distance(s, t):
    m, n = len(s), len(t)
    T = [] # create empty matrix
    T.append([i for i in range(m+1)]) # init row 0
    for i in range(1, n +1):
        T.append([None] * (m +1)) # add empty row i
        T[i][0] = i # init 0 th column of row i
        for j in range(1, m +1):
            T[i][j] = min(
                T[i-1][j-1]+(0 if s[j-1]== t[i-1] else 1) ,
                T[i-1][j] + 1,
                T[i][j-1] + 1
            )
    return T[n][m]
