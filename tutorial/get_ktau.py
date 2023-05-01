# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 13:02:53 2017
@author: disvr
traceback.py
"""
import sys, os
import scipy.stats as stats

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

def get_anchor_list(fname):
    anchorlist = []
    with open(fname) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            anchorlist.append(int(li[1]))
    return anchorlist

statsfname = 'ktau_stats.csv'
g = open(statsfname,"w")

nG = 8
q_mod = 0.9
q_con = 0.3
folder = 'illus3_n8'
fname = 'nx_run=1_qmod=0.9_qcon=0.3_n=8'

extantFNAME = folder + '/' + 'extant'
greedyFNAME = folder + '/' + 'greedy_' + fname

truth = range(nG,0,-1)

anchorlist = get_anchor_list(greedyFNAME)
diff = set(truth) - set(anchorlist)
if len(diff) != 0:
    assert(len(diff) == 2)
    anchorlist.extend(diff)
#print(anchorlist)
#print(truth)
tau_greedy, pval_greedy = stats.kendalltau(anchorlist, truth)
#print(tau_ILP, pval_ILP)	
g.write("{},{},{}\n".format(greedyFNAME, round(tau_greedy,3), round(pval_greedy,3)))
	
for numsol in range(30):
    ilpARlist = folder + '/ARlist_ILP_' + fname + '_solution' + str(numsol) + '.txt'
    if os.path.isfile(ilpARlist):
	anchorlist = get_anchor_list(ilpARlist)
	diff = set(truth) - set(anchorlist)
	#print(diff)
	if len(diff) != 0:
	    assert(len(diff) == 2)
	    anchorlist.extend(diff)

	#print(anchorlist)
	#print(truth)
	tau_ILP, pval_ILP = stats.kendalltau(anchorlist, truth)
	#print(tau_ILP, pval_ILP)	
        g.write("{},{},{}\n".format(ilpARlist, round(tau_ILP,3), round(pval_ILP,3)))
g.close()
