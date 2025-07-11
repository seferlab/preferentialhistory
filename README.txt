For a tutorial on usage see the folder 'Tutorial'

------
ILP-PA.py
------

Program arguments:
usage: ILP-PA.py [-h] -e EXTANT -o OUTPUT  [-p P_VALUE] [-t TIMELIMIT]
              [-s NUMSOLUTIONS] [-m POOLMODE] [-f FOCUS] 
              [-n NUMCORES]  [u - memusage]
ILP to reconstruct maximum likelihood within the Preferential Attachment model. Uses Gurobi solver, see
http://www.gurobi.com/documentation/8.0/refman/index.html

Required arguments:
  -e EXTANT, --extant EXTANT
                        Extant filename. Format: one edge per line
  -o OUTPUT, --output OUTPUT
                        Output Filename
Optional arguments:
  -h, --help            show this help message and exit
  -p P_VALUE, --p_value P_VALUE
                        ILP-PA Model Parameter p (default: 0.5)
  -t TIMELIMIT, --timelimit TIMELIMIT
                        Time Limit in hours (default: 2)
  -n NUMCORES, --numcores NUMCORES
                        Number of cores. 0 uses all available cores (default:
                        0)
  -s NUMSOLUTIONS, --numsolutions NUMSOLUTIONS
                        Number of ILP Solutions (default: 30)
  -m POOLMODE, --poolmode POOLMODE
                        Pool mode, possible values 0,1,2 (default: 1). 0:
                        finds one optimal solution. 1: find multiple solutions
                        not necessarily the best. 2: find n best multiple
                        solutions.
  -f FOCUS, --focus FOCUS
                        MIP focus, possible values 1,2,3 (default: 1). 1:
                        finds feasible solutions quickly. 2: to prove
                        optimality, if good quality solutions can be found
                        easily. 3: if the best objective bound is moving very
                        slowly or not at all.
  -u MEMUSAGE --memusage MEMUSAGE 
                        Memory Usage Limit in GBytes (default: 10)"