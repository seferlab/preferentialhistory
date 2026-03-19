# -*- coding: utf-8 -*-
# from subprocess import check_output
from ilp_pa import run_ILP
import glob, os

for_t2 = [50, 100, 250, 500, 1000]
p = 0.5

for n in for_t2:
    folder = f'clusters_for_table_2/nx_p=0.5_n={n}/nx_p=0.5_n={n}_1'
    for fpath in sorted(glob.glob(os.path.join(folder, 'cluster-*.edges'))):
        fname = os.path.splitext(os.path.basename(fpath))[0]   
        extantFNAME = fpath                                    
        outputFNAME = folder + '/ILP_' + fname                 

        run_ILP(extantFNAME, outputFNAME, p, timeLimit = 4
                , numSolutions=10,
                poolMode=1, focus=1, numThreads=0, memusage=10)
  