"""
Input: 
1. List of Reconstructed networks (filenames)
2. List of True networks (filenames)
Both lists must be matching #nodes

- read nw from list of files
- compare i-th nw from G3 to Gn-1: convert networkx to grakel and get kernel 
- output median or mean normalized by n
"""
import numpy as np
from grakel import GraphKernel
from copy import deepcopy
import os

def relabel(el, nodes):
    nodeli = sorted(nodes)
    #print(nodeli)
    labelmap = {}
    for i in range(1,len(nodes)+1):
        labelmap[nodeli[i-1]] = i
    #print(labelmap)
    newEL = []
    for i,j in el:
        newEL.append((labelmap[i],labelmap[j]))
    return newEL

def get_el(fname):
    el = []
    nodes = []
    #MAX, MIN = 0, 1
    with open(fname) as f:
        for line in f:
            #print(line)
            li = [int(x.strip()) for x in line.split()]
            assert(len(li) == 2)
            i, j = li
            el.append((i,j))

            #if i > MAX:
            #    MAX = i
            #if j > MAX:
            #    MAX = j
            #assert(i >= MIN)
            #assert(j >= MIN)
            
            if i not in nodes:
                nodes.append(i)
            if j not in nodes:
                nodes.append(j)

    #print(nodes)
    #print(MAX)
    #print(el)
    #assert(len(nodes) == MAX)
    return relabel(el, nodes), len(nodes)

def create_adj(N, edge_list):
    adj = np.zeros((N,N),dtype=int)
    for i,j in edge_list:
        adj[i-1,j-1] = 1
    return adj.tolist()

def compute_sim(N1,N2):
    wl_kernel = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": 5}, {"name": "subtree_wl"}], normalize=True)
    wl_kernel.fit(N1)
    sim = wl_kernel.transform(N2)
    return sim[0][0] 

def get_kernel_stats(fnli1, fnli2):
    Nli1, Nli2 = [], []
    for fname in fnli1:
        el, N = get_el(fname)
        L = create_adj(N, el)
        labels = {}
        for i in range(1,N+1):
            labels[i-1] = str(i)    
        Nli1.append([[L, deepcopy(labels)]])
    
    for fname in fnli2:
        el, N = get_el(fname)
        L = create_adj(N, el)
        labels = {}
        for i in range(1,N+1):
            labels[i-1] = str(i)    
        Nli2.append([[L, deepcopy(labels)]])
    
    assert(len(Nli1)==len(Nli2))
    
    simLi = []
    for i in range(len(Nli1)):
        sim = compute_sim(Nli1[i],Nli2[i])
        simLi.append(sim)

    mu, med, std = np.mean(simLi), np.median(simLi), np.std(simLi)
    print(mu,med,std)
    return mu, med, std, simLi

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


def get_true_fnames(folder,run,q_mod,q_con,nG):
    fnli = []
    for i in range(3,nG):
        fname = folder + '/nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(i)
        #fname = folder + '/nx_run='+str(run)+'_nG='+str(nG)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_perc='+str(perc)+'_n='+str(i)
        fnli.append(fname)
    return fnli

statsfname = 'gr_stats.csv'
g = open(statsfname,"w")

nG = 8
q_mod = 0.9
q_con = 0.3
run = 1
folder = 'illus3_n8'
fname = 'nx_run=1_qmod=0.9_qcon=0.3_n=8'

extantFNAME = folder + '/' + 'extant'
greedyFNAME = folder + '/' + 'greedy_' + fname

fnli1 = get_true_fnames(folder,run,q_mod,q_con,nG)
fnli2 = []
for i in range(3,nG):
    fnli2.append(greedyFNAME + '_recon_G'+str(i))
assert(len(fnli1)==len(fnli2))
#print(fnli1)
#print(fnli2)
mu_greedy, med_greedy, std_greedy, li_greedy = get_kernel_stats(fnli1, fnli2)
g.write("{},{},{},{},{}\n".format(greedyFNAME,mu_greedy,med_greedy,std_greedy,li_greedy))

for numsol in range(30):
    ARfname = folder + '/ARlist_ILP_' + fname + '_solution' + str(numsol) + '.txt'
    if os.path.isfile(ARfname):
        print(ARfname)
        fnli2 = []
        for i in range(3,nG):
            fnli2.append(ARfname + '_recon_G'+str(i))
        assert(len(fnli1)==len(fnli2))
        #print(fnli1)
        #print(fnli2)
        mu_ILP, med_ILP, std_ILP, li_ILP = get_kernel_stats(fnli1, fnli2)
	g.write("{},{},{},{},{}\n".format(ARfname, mu_ILP,med_ILP,std_ILP,li_ILP))
g.close()


"""
probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
count1, count2 = 0,0
for run in range(1,6):
    folder = 'run'+str(run)
    for q_con in probL:
        for q_mod in probL:
            for n in range(4,11):
                fname = 'nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)
                #extantFNAME = folder +'/'+ fname
                fnli1 = get_true_fnames(folder,run,q_mod,q_con,n)
                
                ARfname = folder + '/ARlist_'+fname
                if os.path.isfile(ARfname):
                    print(ARfname)
                    fnli2 = []
                    for i in range(3,n):
                        fnli2.append(ARfname + '_recon_G'+str(i))
                    assert(len(fnli1)==len(fnli2))
                    #print(fnli1)
                    #print(fnli2)
                    mu_ILP, med_ILP, std_ILP, li_ILP = get_kernel_stats(fnli1, fnli2)
                    
                    greedyFNAME = folder +'/greedy_'+ fname
                    fnli2 = []
                    for i in range(3,n):
                        fnli2.append(greedyFNAME + '_recon_G'+str(i))
                    assert(len(fnli1)==len(fnli2))
                    #print(fnli1)
                    #print(fnli2)
                    mu_greedy, med_greedy, std_greedy, li_greedy = get_kernel_stats(fnli1, fnli2)
                    #raw_input()
                    g.write("{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(run,q_con,q_mod,n,0.0,mu_greedy,med_greedy,std_greedy,mu_ILP,med_ILP,std_ILP,li_greedy,li_ILP))
g.close()
"""



#Network FILES AT EACH STEP - FROM TRACEBACK.PY
#fnli1 = ['nx_run=1_qmod=0.9_qcon=0.9_n=9']
#fnli2 = ['nx_run=1_qmod=0.9_qcon=0.9_n=9']



"""
N = 7
labels = {}
for i in range(1,N+1):
    labels[i-1] = str(i)

el1 = [(1,2),(1,3),(2,5),(2,4),(3,4),(3,7),(5,6),(6,7)]
L1 = create_adj(N, el1)
#print(L1)
N1 = [[L1,labels]]
#print(N1)

el2 = [(1,2),(1,3),(2,3),(2,4),(3,4),(3,7),(4,6),(5,6),(6,7)]
L2 = create_adj(N, el2)
#print(L2)
N2 = [[L2,labels]]
#print(N2)

sim = compute_sim(N1,N2)
print(sim)
"""
