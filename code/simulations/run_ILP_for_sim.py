# -*- coding: utf-8 -*-
from subprocess import check_output

probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]    
nrange = range(7,11)
runrange = range(1,6)

for n in nrange:
    for run in runrange:
        folder = 'run'+str(run)
        q_con = random.choice(probL) #for q_con in probL:
        q_mod = random.choice(probL) #for q_mod in probL:
        #for q_con in probL:
	#    for q_mod in probL:
        fname = 'nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)
        extantFNAME = folder +'/'+ fname
        greedyFNAME = folder +'/greedy_'+ fname
        outputFNAME = folder +'/ILP_'+ fname
        print(fname)
        run_ILP(extantFNAME, q_con, q_mod, greedyFNAME, outputFNAME, timeLimit, numThreads)
        cmd = ['python','ILP.py',str(q_mod),str(q_con),extantFNAME, greedyFNAME, outputFNAME]
        output = check_output(cmd)

