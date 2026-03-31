import math
import networkx as nx
from gurobipy import Model, GRB, quicksum, LinExpr
import re
import numpy as np
import argparse
rng = np.random.default_rng()


def read_input(extantFNAME):    
    extant_edges = []
    with open(extantFNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            extant_edges.append((int(li[0]),int(li[1])))
    
    extant = nx.Graph()
    extant.add_edges_from(extant_edges)
    
    return extant_edges, extant

RE_X = re.compile(r"\[(\d+),\s*(\d+)\]")



def run_ILP(extantFNAME,outputFNAME, p ,numSolutions, poolMode, focus, numThreads, memusage):
    
    extant_edges, extant = read_input(extantFNAME)
    V = sorted(extant.nodes()) 
    E = list(extant.edges())
    d = {node: set(extant.neighbors(node)) for node in V}
    if len(E) == 1:
        with open(f"{outputFNAME}_order{0}.txt", "w") as f:               
            edges_str = [f"{u} {v}" for (u, v) in E]
            for edge in edges_str:
                f.write(edge + "\n")
        return None
    # ==========Model==========
    
    m = Model("ilp-pa")
   
    edges = len(E)
    time = len(E)
    n = len(V)
    degree = n -1

    m.Params.Threads = numThreads
    m.Params.PoolSearchMode = poolMode
    m.Params.PoolSolutions = numSolutions
    m.Params.NodefileStart = memusage
    # No time limit - solve until optimal or interrupted
    # ========Variables========

    range_for_x = [(e, t) for e in range(edges) for t in range(time)]
    x = m.addVars(range_for_x, vtype=GRB.BINARY, name="x+f'{e}_{t}'")

    range_for_y = [(i, k, t) for i in range(n) for k in range(degree) for t in range(time)]
    y = m.addVars(range_for_y, vtype=GRB.BINARY, name="y+f'{i}_{k}_{t}'")

    range_for_z = [(e, k, t) for e in range(edges) for k in range(degree) for t in range(time)]
    z = m.addVars(range_for_z, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="z+f'{e}_{k}_{t}'")

    print(E)
    # ========Objective========

    # Equation (19)
    obj_expr = LinExpr();
    for t in range(1, time - 1):
        for e in range(edges):
            i = E[e][0] - 1
            j = E[e][1] - 1
            for k in range(degree):
                obj_expr += z[e, k, t+1] * math.log(((p * k) / (2 * t)) + ((1 - p) / (t + 1)))

    m.setObjective(obj_expr, sense = GRB.MAXIMIZE)
    
    # ========Constraints========

    # Equation (6) 
    for t in range(time):
        m.addConstr(quicksum(x[e, t] for e in range(edges)) == 1, name="sum (x+f'{e}_{t}') == 1")

    # Equation (7)
    for e in range(edges):
        m.addConstr(quicksum(x[e, t] for t in range(time)) == 1, name="sum (x+f'{e}_{t}') == 1")    

    # Equation (8)
    for t in range(time):
        for i in range(n):
            m.addConstr(quicksum(y[i, k, t] for k in range(degree)) == 1, name="sum (y+f'{i}_{k}_{t}') == 1")
    
    # Equation (9, 10)
    for t in range(time - 1):    
        for e in range(edges):
            i = E[e][0] - 1
            j = E[e][1] - 1
            m.addConstr(x[e, t+1] <= y[i, 0, t] + y[j, 0, t] , name="x+f'{e}_{t+1}' <= y+f'{i}_0_{t}' + y+f'{j}_0_{t}'")
            m.addConstr(x[e, t+1] <= (1 - y[i, 0, t]) + (1 - y[j, 0, t]) , name="x+f'{e}_{t+1}' <= (1 - y+f'{i}_0_{t}') + (1 - y+f'{j}_0_{t}')")

    # Equation (11)
    for i in range(n):
        neighbours = d[i + 1]
        sum_cons_11 = 0
        for e in neighbours:
            edge = (e , i + 1) if (e , i + 1) in E else (i + 1,e)
            ind = E.index(edge) 
            sum_cons_11 += x[ind,1]
        m.addConstr(y[i, 1, 1] == sum_cons_11, name="y+f'{i}_1_1' == sum (x+f'{e}_1' for e in neighbours)")

    # Equation (12, 13, 14, 15, 16)
    for t in range(time - 1):
        for k in range(t + 1):
            for i in range(n):
                neighbours = d[i + 1]
                sum__cons_12_13_14 = 0
                for e in neighbours:
                    edge = (e , i + 1) if (e , i + 1) in E else (i + 1,e)
                    ind = E.index(edge)
                    sum__cons_12_13_14 += x [ind,t] 
                if k < degree:
                    m.addConstr(y[i, k, t + 1] <= 2 - (y[i, k, t] + sum__cons_12_13_14), name="y+f'{i}_{k}_{t+1}' <= 2 - (y+f'{i}_{k}_{t}' + sum (x+f'{e}_{t}' for e in neighbours))")
                    m.addConstr(y[i, k, t + 1] >= y[i, k, t] - sum__cons_12_13_14, name="y+f'{i}_{k}_{t+1}' >= y+f'{i}_{k}_{t}' - sum (x+f'{e}_{t}' for e in neighbours)")
                    m.addConstr(y[i, k, t + 1] <= y[i, k, t] + sum__cons_12_13_14, name="y+f'{i}_{k}_{t+1}' <= y+f'{i}_{k}_{t}' + sum (x+f'{e}_{t}' for e in neighbours)")
                    if k > 0:
                        m.addConstr(y[i, k, t + 1] <= y[i, k, t] + y[i, k - 1, t], name="y+f'{i}_{k}_{t+1}' <= y+f'{i}_{k}_{t}' + y+f'{i}_{k-1}_{t}")
                        m.addConstr(y[i, k, t + 1] >= y[i, k - 1, t] + sum__cons_12_13_14 - quicksum(1 for e in neighbours), name="y+f'{i}_{k}_{t+1}' >= y+f'{i}_{k-1}_{t}' + sum (x+f'{e}_{t}' - 1 for e in neighbours)")
    
    # Equation (20, 21)
    for t in range(time - 1):
        for e in range(edges):
            i = E[e][0] - 1
            j = E[e][1] - 1
            for k in range(1, degree):
                m.addConstr(z[e, k, t + 1] <= x[e, t + 1], name="z+f'{e}_{k}_{t+1}' <= x+f'{e}_{t+1}'")
                m.addConstr(z[e, k, t + 1] <= y[i, k, t] + y[j, k, t], name="z+f'{e}_{k}_{t+1}' <= y+f'{i}_{k}_{t}' + y+f'{j}_{k}_{t}'")
                # m.addConstr(z[e, k, t + 1] == x[e, t + 1]*(y[i, k, t] + y[j, k, t]),
                #                     name="z+f'{e}_{k}_{t+1}' == x+f'{e}_{t+1}* (y+f'{i}_{k}_{t}' + y+f'{j}_{k}_{t})'")
        
    m.optimize()
    if m.status == GRB.Status.OPTIMAL:
        write_selected_edges(m, E, outputFNAME)
    else:
        if m.Status in [GRB.INF_OR_UNBD, GRB.INFEASIBLE]:
            print("Model is infeasible or unbounded")
            m.computeIIS()
            for c in m.getConstrs():
                if c.IISConstr:
                    print(f"Infeasible constraint: {c.ConstrName}")
            m.write('model_ilp_debug.mps')
        else:
            print(f"Optimization ended with status {m.Status}")
    return m

def write_selected_edges(m, candidate_edges, outputFNAME):
    nSolutions = getattr(m, "SolCount", 0)
    for e in range(nSolutions):
        m.setParam(GRB.Param.SolutionNumber, e)
        picked = []
        for v in m.getVars():
            name = v.VarName
            if not name or name[0] != 'x':
                continue
            if v.Xn < 0.5:
                continue
            mobj = RE_X.search(name)
            if not mobj:
                continue
            ei = int(mobj.group(1))
            t  = int(mobj.group(2))
            if 0 <= ei < len(candidate_edges):
                picked.append((t, candidate_edges[ei]))
        picked.sort(key=lambda z: z[0])
        outp = f"{outputFNAME}_order{e}.txt"
        with open(outp, "w", newline="\n") as f:
            for _, (u, v) in picked:
                f.write(f"{u} {v}\n")

def calculate_edge_history(m, candidate_edges, sol_idx=0):
    if m is None or getattr(m, "SolCount", 0) == 0:
        return {}
    m.setParam(GRB.Param.SolutionNumber, sol_idx)
    hist = {}
    for v in m.getVars():
        name = v.VarName
        if not name or name[0] != 'x':
            continue
        if v.Xn < 0.5:
            continue
        mobj = RE_X.search(name)
        if not mobj:
            continue
        ei = int(mobj.group(1))
        t  = int(mobj.group(2))
        if 0 <= ei < len(candidate_edges):
            hist[t] = candidate_edges[ei]
    return dict(sorted(hist.items()))


def main():

    parser = argparse.ArgumentParser(description="ILP for PA")

    parser.add_argument("-e","--extant",type=str,required=True,help="Extant filename. Format: one edge per line")
    parser.add_argument("-o","--output",type=str,required=True,help="Output Filename")
    parser.add_argument("-p","--p_value",type=float,default=0.5,help="p value for the model (default: %(default)s)")
    # parser.add_argument("-t","--timelimit",type=int,default=2,help="Time Limit in hours (default: %(default)s)")
    parser.add_argument("-u","--memusage",type=float,default=10,help="Memory Usage Limit in GBytes (default: %(default)s)")
    parser.add_argument("-n","--numcores",type=int,default=0,help="Number of cores. 0 uses all available cores (default: %(default)s)")
    parser.add_argument("-s","--numsolutions",type=int,default=10,help="Number of ILP Solutions (default: %(default)s)")
    parser.add_argument("-m","--poolmode",type=int,default=1,help="Pool mode, possible values 0,1,2 (default: %(default)s). 0: finds one optimal solution. 1: find multiple solutions not necessarily the best. 2: find n best multiple solutions.")
    parser.add_argument("-f","--focus",type=int,default=1,help="MIP focus, possible values 1,2,3 (default: %(default)s). 1: finds feasible solutions quickly. 2: to prove optimality, if good quality solutions can be found easily. 3: if the best objective bound is moving very slowly or not at all.")

    args,_= parser.parse_known_args()

    extantFNAME = args.extant 
    outputFNAME = args.output 
    p = args.p_value
    # timeLimit = args.timelimit * 3600  # Convert hours to seconds
    numSolutions = args.numsolutions
    poolMode = args.poolmode
    focus = args.focus
    numThreads = args.numcores
    memusage = args.memusage
    model= run_ILP(extantFNAME
                    , outputFNAME
                    , p
                    # , timeLimit
                    ,numSolutions
                    , poolMode
                    , focus 
                    , numThreads
                    , memusage

    )
    candidate_edges, _ = read_input(extantFNAME)

    edge_history = calculate_edge_history(model, candidate_edges, sol_idx=0)
    return edge_history

if __name__ == "__main__": 
    main()
    
    
        