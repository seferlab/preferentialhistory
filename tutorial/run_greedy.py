# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 09:56:45 2017

@author: disvr
"""
from subprocess import check_output
import sys

def write_graph(filename, elist):
    with open(filename, 'w') as f:
        for x,y in elist:
            f.write(str(x)+' '+str(y)+'\n')

fname = sys.argv[1]
q_con = sys.argv[2]
q_mod = sys.argv[3]
ofname = sys.argv[4]

cmd = ['python','archaeology.py','-m','dmc',fname,str(q_mod),str(q_con)]
output = check_output(cmd)
edges = []
#print(output)
for s in output.splitlines()[1:]:
    #print(s)
    li = s.split()
    assert(len(li) == 2)
    edges.append([int(x) for x in li])
                    
write_graph(ofname,edges)
