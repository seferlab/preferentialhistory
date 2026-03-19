# -*- coding: utf-8 -*-
# from ilp_pa import run_ILP
from ilp_pa_1 import run_ILP

probL = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
nrange = range(6,11)
runrange = range(1,34)

for run in runrange:
    for n in nrange: 
        folder = 'run'+str(run)
        for p in probL:
            fname = 'nx_run='+str(run)+'_p='+str(p)+'_n='+str(n)
            extantFNAME = folder + '/' + fname
            outputFNAME = folder +'/ILP_'+ fname
            run_ILP(extantFNAME,outputFNAME, p, timeLimit = 2 ,numSolutions=10,
                    poolMode=1,focus=1,numThreads = 0, memusage=10)
            
# for file in ['bzip', 'commander']:
#     folder = 'data/' + file
#     fname = file
#     extantFNAME = folder + '/'+ fname
#     outputFNAME = folder + '/'+ 'ILP_'+ fname
#     run_ILP(extantFNAME,outputFNAME, p=0.5,numSolutions=10,
#             poolMode=1,focus=1,numThreads = 0, memusage=10)


# for file in ['bzip', 'commander']:
#     folder = 'data/' + file
#     fname = file 
#     for p in probL:
#         extantFNAME = folder + '/'+ fname
#         outputFNAME = folder + '/'+ 'ILP_'+ fname + '_p='+str(p)
#         run_ILP(extantFNAME, outputFNAME, p ,numSolutions=10,
#                 poolMode=1,focus=1,numThreads = 0, memusage=10)
        
folder = "figure-3"
fname = "nx_p=0.5_n=6"
extantFNAME = folder + '/'+ fname
outputFNAME = folder + '/ILP_'+ fname + '.txt' 
run_ILP(extantFNAME,outputFNAME, p=0.5, timeLimit = 2 ,numSolutions=10,
                    poolMode=1,focus=1,numThreads = 0, memusage=10)