# -*- coding: utf-8 -*-
# from subprocess import check_output
from ilp_pa import run_ILP
import glob, os

for_t4 = [200]
time_limit_list = [1.7, 17, 51, 68, 85, 136]
p = 0.5

for n in for_t4:
    for t in time_limit_list:
        folder = f'clusters_for_table_4/nx_p=0.5_n={n}/nx_p=0.5_n={n}_1'
        for fpath in sorted(glob.glob(os.path.join(folder, 'cluster-*.edges'))):
            fname = os.path.splitext(os.path.basename(fpath))[0]   
            extantFNAME = fpath                                    
            outputFNAME = folder + '/ILP_' + fname                 
    
            run_ILP(extantFNAME, outputFNAME, p, timeLimit = t
                    , numSolutions=10,
                    poolMode=1, focus=1, numThreads=0, memusage=10)