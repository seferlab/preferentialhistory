
import subprocess
import os

"""
ARlist_ILP_nx_run\=1_qmod\=0.4_qcon\=0.7_n\=9_solution9.txt
ARlist_ILP_nx_run\=1_qmod\=0.4_qcon\=0.7_n\=9_solution2.txt
greedy_nx_run\=1_qmod\=0.4_qcon\=0.7_n\=9
"""

def read_ARlist(fname):
    anchor_rem_list = []
    with open(fname) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            anchor_rem_list.append((li[0],li[1]))
    return anchor_rem_list

nG = 8
folder = 'illus3_n8'
"""
fname = 'greedy_nx_run=1_qmod=0.9_qcon=0.3_n=8'
ofprefix = 'greedy_G'
exofname = 'extant_greedy.png'
nodecolor = 'lightblue'
fname = 'ARlist_ILP_nx_run=1_qmod=0.9_qcon=0.3_n=8_solution0.txt'
ofprefix = 'ILP_0_G'
exofname = 'extant_ILP0.png'
nodecolor = 'lightgreen'
"""
fname = 'ARlist_ILP_nx_run=1_qmod=0.9_qcon=0.3_n=8_solution22.txt'
ofprefix = 'ILP_22_G'
exofname = 'extant_ILP22.png'
nodecolor = 'lightgreen'
anchor_rem_list = read_ARlist(folder + '/' + fname)

fnli = {}
for i in range(3,nG):
    fnli[i] = folder + '/' + fname + '_recon_G'+str(i)

j = 0
gfname = folder + '/'+ 'extant'
anchor, dup = anchor_rem_list[j]
j += 1
if os.path.isfile(gfname):
    cmd = ['python','draw_graph.py',gfname,exofname,str(nG),nodecolor,anchor,dup]
    subprocess.call(cmd)
for i in range(nG-1,2,-1):
    gfname = fnli[i]
    anchor, dup = anchor_rem_list[j]
    j += 1
    if os.path.isfile(gfname):
	ofname = ofprefix + str(i) +'.png'
        cmd = ['python','draw_graph.py',gfname,ofname,str(nG),nodecolor,anchor,dup]
        subprocess.call(cmd)



