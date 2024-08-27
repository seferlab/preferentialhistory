# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 09:56:45 2017

@author: disvr
"""
from subprocess import check_output


def write_graph(filename, elist):
    with open(filename, 'w') as f:
        for x,y in elist:
            f.write(str(x)+' '+str(y)+'\n')

probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]    

for run in range(1,6):
    folder = 'run'+str(run)
    for q_con in probL:
        for q_mod in probL:
            for n in range(4,11):
                fname = 'nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)
                name = folder +'/'+ fname
                print(name)
                cmd = ['python','archaeology.py','-m','dmc',name,str(q_mod),str(q_con)]
                output = check_output(cmd)
                edges = []
                #print(output)
                for s in output.splitlines()[1:]:
                    #print(s)
                    li = s.split()
                    assert(len(li) == 2)
                    edges.append([int(x) for x in li])
                    
                oname = folder+'/greedy_'+fname
                write_graph(oname,edges)