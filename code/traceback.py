# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 13:02:53 2017
@author: disvr
traceback.py
"""

"""
input1: file with (anchor, node removed) at each step of evolution
input2: file with extant graph in "edge per line" format
input3,4: q_con, q_mod

output1: loglikelihood
output2: graphs at each step, each in a file
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
        

#q_mod, q_con = 0.1, 0.1
#run = 1
#folder = 'run'+str(run)
#
#extantFNAME = folder+'/'+'nx_run=1_qmod=0.1_qcon=0.1_n=4'
#greedyFNAME = folder+'/'+'greedy_nx_run=1_qmod=0.1_qcon=0.1_n=4'
#ilpFNAME = folder+'/'+'ILP_nx_run=1_qmod=0.1_qcon=0.1_n=4_solution.txt'

statsfname = 'stats.csv'
g = open(statsfname,"w")
g.write("run,q_con,q_mod,n,LKL_greedy,LKL_ILP\n")
diff = []

probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]    
count1, count2 = 0,0
for run in range(1,6):
    folder = 'run'+str(run)
    for q_con in probL:
        for q_mod in probL:
            for n in range(4,7):#11):
                fname = 'nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)
                extantFNAME = folder +'/'+ fname
                greedyFNAME = folder +'/greedy_'+ fname
                ilpFNAME = folder +'/ILP_'+ fname + '_solution.txt'
                print(fname)
                
                if os.path.isfile(ilpFNAME) :
                    extant, greedy_ARlist = read_input(extantFNAME, greedyFNAME)
                    nx_list = get_all_networks(extant, greedy_ARlist)
                    gtot_lkl, gtot_logLKL = compute_totlkl(nx_list, greedy_ARlist, q_mod, q_con)
                    ilptot_logLKL = getILPlkl(ilpFNAME)
                    
                    ilptot_logLKL = round(ilptot_logLKL,3)
                    gtot_logLKL = round(gtot_logLKL,3)

		    g.write("{},{},{},{},{},{}\n".format(run,q_con,q_mod,n,gtot_logLKL,ilptot_logLKL))

                    if ilptot_logLKL > gtot_logLKL:
                        print('ILP better', ilptot_logLKL, gtot_logLKL)
                        count1 += 1
			diff.append(ilptot_logLKL-gtot_logLKL)
                    elif ilptot_logLKL == gtot_logLKL:
                        print('ILP equal', ilptot_logLKL, gtot_logLKL)
                        count2 += 1
                    else:
                        print('Greedy better. How?', ilptot_logLKL, gtot_logLKL)
                        boo
print count1, count2
g.close()
import numpy as np
print len(diff), np.mean(diff), np.std(diff), np.max(diff), np.min(diff), np.median(diff)
        

