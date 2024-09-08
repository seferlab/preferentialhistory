#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 07:57:51 2018

@author: nus
"""
import copy, sys, math
def clean(x):
    if 'e' in x:
        c = int(float(x))
    else:
	    c = int(math.ceil(float(x)))
    assert(c == 0 or c == 1)
    return c

def process_anchor(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, node = [int(x) for x in li[0].lstrip('NodeAnchor[').rstrip(']').split(',')]
        soln[1][graph] = node
    return soln
    
def process_nodex(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, node = [int(x) for x in li[0].lstrip('NodeX[').rstrip(']').split(',')]
        soln[2][graph].append(node)
    return soln
    
def process_phantom(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, ni, nj = [int(x) for x in li[0].lstrip('Phantom[').rstrip(']').split(',')]
        soln[0][graph][nj] = ni #nj in graph, ni in prev graph
    return soln

def read_file(fname, soln):
    with open(fname) as f:
        for line in f:
            if line.startswith('NodeAnchor'):
                # anchor node
                soln = process_anchor(line, soln)
            elif line.startswith('NodeX'):
                # duplicated nodes
                soln = process_nodex(line, soln)
            elif line.startswith('Phantom'):
                soln = process_phantom(line, soln)
    return soln
                
def create_ds(n):
    #n: number of nodes in extant
    #for each graph, we create a dict that maps nodes through phantom edge
    dictlist = []
    for i in range(n+1): #dummy 0,1,2 -> no graphs
        phantommap = {}
        for j in range(1,i+1):
            phantommap[j] = -1 #will be assigned when we read the solution file
        dictlist.append(copy.deepcopy(phantommap))
    #for each graph we have a dict of anchor and X nodes
    anchordict, Xdict = {}, {}
    for i in range(2, n):
        anchordict[i] = -1
    for i in range(3, n+1):
        Xdict[i] = []
        
    soln = (dictlist, anchordict, Xdict)
    return soln

def get_anchor_list(n, f, soln):
    """\
        Find (anchor,remove) node list in an reversed order (from last to the first), 
        node indices all of the last graph
        n: 
            last graph size
        f:
            first graph size
        soln:
            solution
    """
    # n graph number 
    ARlist = []
    Rlist = []

    # n th solution
    # 0 phantom edge, 1 anchor node, 2 duplicated node
    # ar a list of length 2, duplicated nodes
    ar = soln[2][n]
    # smaller one or larger one, user chosen 
    # ar.reverse()
    # inverse order from large to small, ARlist anchor remove list, node index in last graph
    ARlist.append(ar)   
    Rlist.append(ar[-1])
    # recursively find the pair before
    #get current mapping for (n-1)-st graph
    mapping = {} 

    for j in range(1,n+1):
        # except the removed node, find the mapping from previous graph to last graph
        if j not in Rlist:
            # n = 7, for each node j in 7th graph, find the corresponding previous node in the 6th graph
            # mapping[pre_g] -> g, phantom used as indication
            mapping[soln[0][n][j]] = j
    #traverse the graphs in reverse order, until the first graph
    for g in range(n-1, f, -1):
        # e.g. g = 6, duplicated nodes in 6, map to nodes in 7
        ar = [mapping[x] for x in soln[2][g]] 
        # delete
        # ar.reverse()
        ARlist.append(ar)
        Rlist.append(ar[-1])
        # find new mapping from graph to last graph, such that the anchor remove pair are all of the index of the last graph
        newmap = {}
        for j in range(1,g+1):
            if mapping[j] not in Rlist:
                newmap[soln[0][g][j]] = mapping[j]
        mapping = newmap
    
    return ARlist


def get_all_networks(extant, anchor_rem_list):
    """\
        Calculate the graph for each level, key is to get the ordered anchr_rem_list, from the get_anchor_list function
        return a list of graph obj
    """
    nx_list =[extant]
    curr_nx = extant.copy()
    # a r are anchor and duplicated(remove) nodes indices in the extant (last) graph
    # the order from the last to the first
    for a,r in anchor_rem_list:
        # remove r node
        Nr = curr_nx.neighbors(r)
        Na = curr_nx.neighbors(a)
        newN = set(Nr) - set(Na)
        newN -= set([a])
        newN -= set([r]) #self loop in r
        curr_nx.remove_node(r)
        # add in new connection from r's neighbors to anchor node
        for n in newN:
            curr_nx.add_edge(a,n) #adds n in graph if it does not exist
        nx_list.append(curr_nx.copy())

        # we don't confine the number of starting nodes to be 2, greedy with 1
        # if curr_nx.number_of_nodes() == 2:
        #     break
    print("number of graphs:", len(nx_list), ", number of anchor remove pairs:", len(anchor_rem_list))
    assert(len(nx_list) == len(anchor_rem_list) + 1)
    return nx_list


def compute_totlkl(nx_list, anchor_rem_list, q_mod, q_con):
    """\
    compute likelihood from reconstructed graph list and ARlist
    """
    tot_lkl, tot_logLKL = 1.0, 0.0
    for i, curr_g in enumerate(nx_list):
        # curr_g the current graph object, with last graph index
        # corresponding (anchor, remove) pair for the graph(before removed)
        u, v = anchor_rem_list[i]
        # calculate the likelihood for the pair
        lkl, logLKL = compute_lkl(u, v, curr_g, q_mod, q_con)
        tot_lkl *= lkl
        tot_logLKL += logLKL
    return tot_lkl, tot_logLKL


def compute_lkl(u, v, curr_g, q_mod, q_con, constants = False):
    """\
    compute the likelihood for each anchor-remove pair
    """
    U = set(curr_g.neighbors(u)) - set([v])
    V = set(curr_g.neighbors(v)) - set([u])
    
    gamma = q_con if v in curr_g.neighbors(u) else 1-q_con
    Cuv = len(U & V)  #intersection: Common neighbors
    Suv = len(U) + len(V) - 2*Cuv  #symm diff: Single neighbors
    
    if constants:
        LKL = (1-q_mod)**(Cuv) * (q_mod/2)**(Suv) * gamma
        logLKL = (math.log(1-q_mod)*Cuv) + ((math.log(q_mod)-math.log(2))*Suv) + math.log(gamma)
    else:
        LKL = (1-q_mod)**(Cuv) * (q_mod)**(Suv) * gamma
        logLKL = (math.log(1-q_mod)*Cuv) + (math.log(q_mod)*Suv) + math.log(gamma)

    return LKL, logLKL