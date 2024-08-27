#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:14:28 2019

@author: nus
"""
from subprocess import check_output
import os

qconli = [0.3]
qmodli = [0.9]
n = 8
folder = 'illus3_n8'

for numsol in range(30):
    for q_con in qconli:
        for q_mod in qmodli:
	    FNAME = 'ILP_nx_run=1_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)+'_solution'+str(numsol)+'.txt'
	    ilpFNAME = folder +'/' + FNAME

	    if os.path.isfile(ilpFNAME):
	        try:
	    	    cmd = ['python','get_anchor_list2.py',str(n),ilpFNAME]
		    output = check_output(cmd)
		    ofname = folder+'/ARlist_'+FNAME
		    with open(ofname,'w') as f:
		        f.write(output)
		except:
		    print(ilpFNAME)
		    #pass
