#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 11:49:13 2019

@author: nus
"""
import networkx as nx
import os
import sys

def read_ARlist(fname):
    anchor_rem_list = []
    with open(fname) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            anchor_rem_list.append((li[0],li[1]))
    return anchor_rem_list

def read_input(extantFNAME):    
    input_edges = []
    with open(extantFNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            input_edges.append((li[0],li[1]))
    
    extant = nx.Graph()
    extant.add_edges_from(input_edges)
    return extant

def get_all_networks(extant, anchor_rem_list):
    nx_list =[extant]
    curr_nx = extant.copy()
    #print(extant.edges())
    for a,r in anchor_rem_list:
        #print('A R #nodes',a,r,curr_nx.number_of_nodes())
        Nr = curr_nx.neighbors(r)
        Na = curr_nx.neighbors(a)
        newN = set(Nr) - set(Na)
        newN -= set([a])
        curr_nx.remove_node(r)
        for n in newN:
            curr_nx.add_edge(a,n)
        #print(curr_nx.edges())
        nx_list.append(curr_nx.copy())
        if curr_nx.number_of_nodes() == 3:
            break
    #print(len(nx_list),len(anchor_rem_list))
    #assert(len(nx_list) == len(anchor_rem_list)-2)
    return nx_list

def write_files(nx_list, fname, N):
    flist = []
    for i, gr in enumerate(nx_list[1:],1):
        with open(fname+'_G'+str(N-i), 'w') as f:
            for j,k in list(gr.edges):
                f.write(str(j)+' '+str(k)+'\n')
            flist.append(fname+'_G'+str(N-i))
    return flist

def extract_params(fname):
    li = fname.split('_')
    for x in li:
        if x.startswith('qmod'):
            qmod = x.split('=')[1]
        if x.startswith('qcon'):
            qcon = x.split('=')[1]
        if x.startswith('perc'):
            perc = x.split('=')[1]
        if x.startswith('nG'):
            nG = x.split('=')[1]
    return float(qmod), float(qcon), float(perc), int(nG)

n = int(sys.argv[1])
ARfname = sys.argv[2]
extantFNAME = sys.argv[3]

extant = read_input(extantFNAME)

if os.path.isfile(ARfname):
    print(ARfname)
    anchor_rem_list = read_ARlist(ARfname)
    nx_list = get_all_networks(extant, anchor_rem_list) 
    fname_prefix = ARfname + '_recon'
    flist = write_files(nx_list, fname_prefix, n)
