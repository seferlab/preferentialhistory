# -*- coding: utf-8 -*-
from subprocess import check_output

probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]    
nrange = range(7,11)
runrange = range(1,6)

for n in nrange:
    for run in runrange:
        folder = 'run'+str(run)
        p = random.choice(probL) #for p in probL:
        fname = 'nx_run='+str(run)+'p='+str(p)+'_n='+str(n)
        extantFNAME = folder +'/'+ fname
        outputFNAME = folder +'/ILP_'+ fname
        print(fname)
        run_ILP(extantFNAME,outputFNAME, p, timeLimit ,numThreads)
        cmd = ['python','ILP-PA.py',str(p),extantFNAME, outputFNAME]
        output = check_output(cmd)