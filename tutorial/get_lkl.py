# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 13:02:53 2017
@author: disvr
traceback.py
"""

import networkx as nx
import math
import os.path

def compute_lkl(u, v, curr_g, q_mod, q_con, constants = False):
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

def read_ARlist(fname):
    anchor_rem_list = []
    with open(fname) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            anchor_rem_list.append((li[0],li[1]))
    return anchor_rem_list

def read_input(extantFNAME, greedyFNAME):    
    input_edges = []
    with open(extantFNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            input_edges.append((li[0],li[1]))
    
    extant = nx.Graph()
    extant.add_edges_from(input_edges)
    
    greedy_ARlist = read_ARlist(greedyFNAME)
            
    return extant, greedy_ARlist

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
                
def compute_totlkl(nx_list, anchor_rem_list, q_mod, q_con):
    tot_lkl, tot_logLKL = 1.0, 0.0
    for i, curr_g in enumerate(nx_list):
        u, v = anchor_rem_list[i]
        #print(i, u, v)
        #print(curr_g.edges())
        lkl, logLKL = compute_lkl(u, v, curr_g, q_mod, q_con)
        #print('LKL:', lkl, logLKL)
        tot_lkl *= lkl
        tot_logLKL += logLKL
    return tot_lkl, tot_logLKL

def getILPlkl(ilpfname):
    with open(ilpfname) as f:
        for line in f:
            if line.startswith('Obj'):
                lkl = float(line.split()[-1].strip())
                break
    return lkl
        

statsfname = 'lkl_stats.csv'
g = open(statsfname,"a")
diff = []

q_mod = 0.9
q_con = 0.3
folder = 'illus3_n8'
fname = 'nx_run=1_qmod=0.9_qcon=0.3_n=8'

extantFNAME = folder + '/' + 'extant'
greedyFNAME = folder + '/' + 'greedy_' + fname
extant, greedy_ARlist = read_input(extantFNAME, greedyFNAME)
nx_list = get_all_networks(extant, greedy_ARlist)
gtot_lkl, gtot_logLKL = compute_totlkl(nx_list, greedy_ARlist, q_mod, q_con)
gtot_logLKL = round(gtot_logLKL,3)
g.write("{},{}\n".format(greedyFNAME,gtot_logLKL))

for numsol in range(30):
    ilpFNAME = folder + '/' + 'ILP_' + fname + '_solution' + str(numsol) + '.txt'
    if os.path.isfile(ilpFNAME):
        ilptot_logLKL = getILPlkl(ilpFNAME)
        ilptot_logLKL = round(ilptot_logLKL,3)
	g.write("{},{}\n".format(ilpFNAME,ilptot_logLKL))
g.close()
