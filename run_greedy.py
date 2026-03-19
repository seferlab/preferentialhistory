# -*- coding: utf-8 -*-
from subprocess import check_output
def write_graph(filename):
    with open(filename, 'w') as f:
        for s in output.splitlines()[1:]:
            s = s.decode().strip()
            if s:
                f.write(s + '\n')
probL = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

runrange = range(1,34)
nrange = range(6,11)

for run in runrange:
    for n in nrange: 
        for p in probL:
            folder = 'run'+str(run)
            fname = 'nx_run='+str(run)+'_p='+str(p)+'_n='+str(n)
            extantFNAME = folder + '/' + fname
            outputFNAME = folder +'/'+'/ILP_greedy_'+ fname + '.txt'
            cmd = ['python', 'archaeology.py', '-m', 'pa', extantFNAME, str(p)] 
            output = check_output(cmd)
            write_graph(outputFNAME)


for_t2 = [50,100,250,500,1000]
folder = "clusters_for_table_2"
for n in for_t2:
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    fname = f'nx_p=0.5_n={n}'
    extantFNAME = folder_new + '/'+ fname
    outputFNAME = folder_new + '/ILP_greedy_'+ fname + '.txt' 
    cmd = ['python', 'archaeology.py', '-m', 'pa', extantFNAME, str(0.5)] 
    output = check_output(cmd)
    write_graph(outputFNAME)
     
for_t3 = [12, 20, 50, 100]
folder = "clusters_for_table_3"
for n in for_t3:
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    fname = f'nx_p=0.5_n={n}'
    extantFNAME = folder_new + '/'+ fname
    outputFNAME = folder_new + '/ILP_greedy_'+ fname + '.txt' 
    cmd = ['python', 'archaeology.py', '-m', 'pa', extantFNAME, str(0.5)] 
    output = check_output(cmd)
    write_graph(outputFNAME)

for_t4 = [200]
folder = "clusters_for_table_4"
for n in for_t4:
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    fname = f'nx_p=0.5_n={n}'
    extantFNAME = folder_new + '/'+ fname
    outputFNAME = folder_new + '/ILP_greedy_'+ fname + '.txt' 
    cmd = ['python', 'archaeology.py', '-m', 'pa', extantFNAME, str(0.5)] 
    output = check_output(cmd)
    write_graph(outputFNAME)
    
for file in ['bzip', 'commander']:
    folder = 'data/' + file
    fname = file
    extantFNAME = folder + '/'+ fname
    outputFNAME = folder + '/ILP_greedy_'+ fname + '.txt' 
    cmd = ['python', 'archaeology.py', '-m', 'pa', extantFNAME, str(0.5)] 
    output = check_output(cmd)
    write_graph(outputFNAME)
     