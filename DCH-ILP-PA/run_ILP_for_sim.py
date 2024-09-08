try:
    from subprocess import check_output
except:
    print("cannot import subprocess")

#change depending on compute power
numcores = 4 #num of cores
memusage = 32 #max memory usage in GB

numsol = 1 #number of solutions

run = 1
q_con = 0.4
q_mod = 0.7
# time limit in hour
lim = 6

folder = './Example_data/nG50'

# totally 8 subgraphs, run each sub-graph one by one, or run in parallel given enough computational resources
sub_graph = [0, 1, 2, 3, 4, 5, 6, 7]
for i in sub_graph:
    fname = 'group'+ str(i)
    extantFNAME = folder +'/'+ fname + '.json'
    outputFNAME = folder +'/ILP_nogreedy_'+ fname + '_tlim='+str(lim)
    cmd = ['python','./ILP.py','-e',extantFNAME,'-o',outputFNAME,'-c',str(q_con),'-m',str(q_mod),'-t',str(lim),'-s',str(numsol),'-n',str(numcores), '-u', str(memusage)]
    print("start: " + fname)
    output = check_output(cmd)