from ilp_pa_1 import run_ILP

for file in ['bzip', 'commander']:
    folder = 'data/' + file
    fname = file
    extantFNAME = folder + '/'+ fname
    outputFNAME = folder + '/'+ 'ILP_'+ fname
    run_ILP(extantFNAME,outputFNAME, p=0.5,numSolutions=10,
            poolMode=1,focus=1,numThreads = 0, memusage=10)


for file in ['bzip', 'commander']:
    folder = 'data/' + file
    fname = file 
    for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        extantFNAME = folder + '/'+ fname
        outputFNAME = folder + '/'+ 'ILP_'+ fname + '_p='+str(p)
        run_ILP(extantFNAME, outputFNAME, p ,numSolutions=10,
                poolMode=1,focus=1,numThreads = 0, memusage=10)