# -*- coding: utf-8 -*-
from ilp_pa import run_ILP
import glob, os

p = 0.5

for n in [50, 100, 250, 500, 1000]:
    folder = f'clusters_for_table_2/nx_p=0.5_n={n}/nx_p=0.5_n={n}_1'
    for fpath in sorted(glob.glob(os.path.join(folder, 'cluster-*.edges'))):
        fname = os.path.splitext(os.path.basename(fpath))[0]   
        extantFNAME = fpath                                    
        outputFNAME = folder + '/ILP_' + fname                 

        run_ILP(extantFNAME, outputFNAME, p, timeLimit = 4
                , numSolutions=10,
                poolMode=1, focus=1, numThreads=0, memusage=10)
  


for n in [12, 20, 50, 100]:
    for t in [1, 2, 4, 8, 16, 32, 64]:
        folder = f'clusters_for_table_3/nx_p=0.5_n={n}/nx_p=0.5_n={n}_1'
        for fpath in sorted(glob.glob(os.path.join(folder, 'cluster-*.edges'))):
            fname = os.path.splitext(os.path.basename(fpath))[0]   
            extantFNAME = fpath                                    
            outputFNAME = folder + '/ILP_' + fname                 
    
            run_ILP(extantFNAME, outputFNAME, p, timeLimit = t
                    , numSolutions=1,
                    poolMode=1, focus=1, numThreads=0, memusage=10)

for n in [200]:
    for t in [1.7, 17, 51, 68, 85, 136]:
        folder = f'clusters_for_table_4/nx_p=0.5_n={n}/nx_p=0.5_n={n}_1'
        for fpath in sorted(glob.glob(os.path.join(folder, 'cluster-*.edges'))):
            fname = os.path.splitext(os.path.basename(fpath))[0]   
            extantFNAME = fpath                                    
            outputFNAME = folder + '/ILP_' + fname                 
    
            run_ILP(extantFNAME, outputFNAME, p, timeLimit = t
                    , numSolutions=10,
                    poolMode=1, focus=1, numThreads=0, memusage=10)
            

