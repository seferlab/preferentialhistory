
import subprocess
import os

qconli = [0.3]
qmodli = [0.9]
n = 8
folder = 'illus3_n8'

extantFNAME = folder + '/extant'

for numsol in range(30):
    for q_con in qconli:
        for q_mod in qmodli:
            ARfname = folder + '/ARlist_ILP_nx_run=1_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)+'_solution'+str(numsol)+'.txt'

            if os.path.isfile(ARfname):
                cmd = ['python','ARlist2GraphList.py',str(n),ARfname,extantFNAME]
                subprocess.call(cmd)
