For a tutorial on usage see the folder 'Tutorial'

------
ILP.py
------

Program arguments:
usage: ILP.py [-h] -e EXTANT -g GREEDY -o OUTPUT [-c Q_CON] [-m Q_MOD]
              [-t TIMELIMIT] [-n NUMCORES] [-s NUMSOLUTIONS] [-p POOLMODE]
              [-f FOCUS] [-r HEURISTICS]

ILP to reconstruct maximum likelihood network evolution history. Uses Gurobi solver, see
http://www.gurobi.com/documentation/8.0/refman/index.html

Required arguments:
  -e EXTANT, --extant EXTANT
                        Extant filename. Format: one edge per line
  -g GREEDY, --greedy GREEDY
                        Filename of Solution from Greedy Approach ReverseDMC
  -o OUTPUT, --output OUTPUT
                        Output Filename

Optional arguments:
  -h, --help            show this help message and exit
  -c Q_CON, --q_con Q_CON
                        DMC Model Parameter q_con (default: 0.7)
  -m Q_MOD, --q_mod Q_MOD
                        DMC Model Parameter q_mod (default: 0.4)
  -t TIMELIMIT, --timelimit TIMELIMIT
                        Time Limit in hours (default: 24)
  -n NUMCORES, --numcores NUMCORES
                        Number of cores. 0 uses all available cores (default:
                        0)
  -s NUMSOLUTIONS, --numsolutions NUMSOLUTIONS
                        Number of ILP Solutions (default: 30)
  -p POOLMODE, --poolmode POOLMODE
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
  -r HEURISTICS, --heuristics HEURISTICS
                        Fraction of time spent on MIP heuristics, between 0
                        and 1 (default: 0.5). Larger values produce more and
                        better feasible solutions, at a cost of slower
                        progress in the best bound.

