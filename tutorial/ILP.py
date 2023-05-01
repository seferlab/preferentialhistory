# -*- coding: utf-8 -*-
import math, sys, random
import networkx as nx
#import matplotlib.pyplot as plt
from gurobipy import *
import argparse

def plot(first_edges,m,nG):
    #plot    
    myG = nx.Graph()
    myG.add_edges_from(first_edges)
    nx.draw(myG, with_labels=True)
    plt.show()
    
    for g in range(3,nG):
        phantom = []
        prevg = g-1                    
        for i in range(1,prevg+1):
            for j in range(1,g+1):
                namestr = 'Phantom['+str(g)+','+str(i)+','+str(j)+']'
                val = m.getVarByName(namestr).x
                if val == 1:
                    print(namestr, val)
                    phantom.append((i,'-->',j))
                
        print('Phantom',phantom)
     
        myG = nx.Graph()
        edges = []
        for i in range(1,g+1):
            for j in range(i+1,g+1):
                namestr = 'Edge['+str(g)+','+str(i)+','+str(j)+']'
                val = m.getVarByName(namestr).x
                #print(namestr, val)
                if val == 1:
                    edges.append((i,j))
        myG.add_edges_from(edges)
        nx.draw(myG, with_labels=True)
        plt.show()
    
    phantom = []           
    prevg, g = nG-1, nG
    
    for i in range(1,prevg+1):
        for j in range(1,g+1):
            namestr = 'Phantom['+str(g)+','+str(i)+','+str(j)+']'
            val = m.getVarByName(namestr).x
            if val == 1:
                print(namestr, val)
                phantom.append((i,'-->', j))
    print('Phantom',phantom)
    
    myG = nx.Graph()
    myG.add_edges_from(extant_edges)
    nx.draw(myG, with_labels=True)
    plt.show()


def get_all_networks(extant, anchor_rem_list):
    nx_list =[extant]
    curr_nx = extant.copy()
    for a,r in anchor_rem_list:
        #print(a,r,curr_nx.number_of_nodes())
        Nr = curr_nx.neighbors(r)
        Na = curr_nx.neighbors(a)
        newN = set(Nr) - set(Na)
        newN -= set([a])
        newN -= set([r]) #self loop in r
        curr_nx.remove_node(r)
        for n in newN:
            curr_nx.add_edge(a,n) #adds n in graph if it does not exist
        nx_list.append(curr_nx.copy())
        if curr_nx.number_of_nodes() == 2:
            break
    print(len(nx_list),len(anchor_rem_list))
    assert(len(nx_list) == len(anchor_rem_list)-1)
    return nx_list

def read_input(extantFNAME, greedyFNAME):    
    extant_edges = []
    with open(extantFNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            extant_edges.append((int(li[0]),int(li[1])))
    
    extant = nx.Graph()
    extant.add_edges_from(extant_edges)
    
    anchor_rem_list = []
    with open(greedyFNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            anchor_rem_list.append((int(li[0]),int(li[1])))
            
    return extant_edges, extant, anchor_rem_list

def run_ILP(extantFNAME, q_con, q_mod, greedyFNAME, outputFNAME, timeLimit, numThreads, numSolutions=30, poolMode=1, focus=1, heuristicstime=0.5):
    #####################################################################################################
    #Create a new model
    m = Model("ReconstructNetwork")
    
    #Model parameters
    m.Params.MIPFocus = focus
    m.Params.Heuristics = heuristicstime
    #m.Params.ImproveStartGap = 2
    m.Params.Threads = numThreads
    m.Params.TimeLimit = timeLimit
    m.Params.PoolSearchMode = poolMode
    m.Params.PoolSolutions = numSolutions
    #m.Params.IntFeasTol = 1e-9
    #m.Params.DualReductions = 0
    #m.Params.Presolve = 0
    #m.Params.ConcurrentMIP = numThreads if numThreads != 0 else 1
    #m.Params.SubMIPNodes = int(sys.maxsize/(2**40))
    
    ##################################################################################################### 
    # Input
    #q_con = 0.7
    #q_mod = 0.4
    #folder = 'C:/Users/DISVR/Desktop/networks - xiuwei/baseline methods/netarch/src/'
    #input_file1, input_file2 = folder+'extant4.txt', folder+'AR4.txt'
    
    extant_edges, extant, anchor_rem_list = read_input(extantFNAME, greedyFNAME)
    nG = len(extant.nodes())
    
    #extant_edges = [(1,2),(2,3),(1,4)]
    #extant_edges = [(1,2),(1,5),(2,5),(3,5),(4,5),(3,4),(5,6)]
    #nG = 4
    #nG = 6
    
    LG = []
    for i in range(nG):
        LG.append([])
        for j in range(nG):
            LG[i].append(0)
    for i,j in extant_edges:
        LG[i-1][j-1] = 1
        LG[j-1][i-1] = 1
    #print(LG)
    
    #LG = [[0,1,0,0,1],[1,0,0,0,1],[0,0,0,0,1],[0,0,0,0,1],[1,1,1,1,0]] #known edges --> extant graph
    FG = [[0,1],[1,0]] #known edges --> first graph
    first_edges = [(1,2)]
    
    #####################################################################################################
    #Construct ILP
    
    ##Variables
    
    #e_ijg edge variables for each graph except extant and first graph
    edges = []
    for g in range(3,nG):
        for i in range(1, g+1):
            for j in range(1, g+1):
                if i != j:
                    edges.append((g,i,j))
    E = m.addVars(edges, vtype=GRB.BINARY, name='Edge')
    #print(E)
    
    #d_ig nodes representing anchors -- to meet path constraint for duplicating nodes -- not for extant
    anchor = []
    for g in range(2,nG):
        for i in range(1,g+1):
            anchor.append((g,i))
    Na = m.addVars(anchor, vtype=GRB.BINARY, name='NodeAnchor')
    #print(Na)
    
    #x_ig
    nodes = []
    for g in range(3,nG+1):
        for i in range(1, g+1):
            nodes.append((g,i))
    Nx = m.addVars(nodes, vtype=GRB.BINARY, name='NodeX')
    #print(Nx)
    
    #y_kg
    nodes = []
    for g in range(3,nG+1):
        for k in range(1, g+1):
            nodes.append((g,k))
    Ny = m.addVars(nodes, vtype=GRB.BINARY, name='NodeY')
    #print(Ny)
    
    #z_lg
    nodes, temp = [], []
    neighbor = []
    for g in range(3,nG+1):
        for l in range(1, g+1):
            nodes.append((g,l))
            temp.append((g,l))
            neighbor.append((g,l))
    Nz = m.addVars(nodes, vtype=GRB.BINARY, name='NodeZ')
    Nw = m.addVars(temp, vtype=GRB.BINARY, name='NodeW')
    Nn = m.addVars(neighbor, vtype=GRB.BINARY, name='NodeN')
    #print(Nz)
    
    #phantom edges
    phantom = []
    for prevg in range(2, nG):
        g = prevg+1
        for i in range(1, prevg+1):
            for j in range(1,g+1):
                #print('phantom',g,i,j)
                phantom.append((g,i,j))
    P = m.addVars(phantom, vtype=GRB.BINARY, name='Phantom')
    
    ###############################################################################
    ##dummy nodes for intermediate values
    
    Dummya, Dummyb, Dummyc = [],[],[]
    Dummycn, Dummyexy, Dummyexz = [], [], []
    DummyP1, DummyP2 = [],[]
    DummyP2a, DummyP2b = [],[]
    DummyEP1, DummyEP2, DummyTP1, DummyTP2 = [],[],[],[]
    for g in range(3,nG+1):
        for i in range(1, g+1):
            for j in range(1, g+1):
                Dummya.append((g,i,j))
                Dummyb.append((g,i,j))
                Dummyc.append((g,i,j))
                Dummycn.append((g,i,j))
                Dummyexy.append((g,i,j))
                Dummyexz.append((g,i,j))
                DummyP1.append((g,i,j))
                DummyP2.append((g,i,j))
                DummyP2a.append((g,i,j))
                DummyP2b.append((g,i,j))
                DummyEP1.append((g,i,j))
                DummyEP2.append((g,i,j))
                DummyTP1.append((g,i,j))
                DummyTP2.append((g,i,j))
    Da = m.addVars(Dummya, vtype=GRB.BINARY, name='DummyA')
    Db = m.addVars(Dummyb, vtype=GRB.BINARY, name='DummyB')
    Dc = m.addVars(Dummyc, vtype=GRB.BINARY, name='DummyC')
    Dcn = m.addVars(Dummycn, vtype=GRB.BINARY, name='DummyCN')
    EXY = m.addVars(Dummyexy, vtype=GRB.BINARY, name='DummyEXY')
    EXZ = m.addVars(Dummyexz, vtype=GRB.BINARY, name='DummyEXZ')
    P1 = m.addVars(DummyP1, vtype=GRB.BINARY, name='DummyP1')
    P2 = m.addVars(DummyP2, vtype=GRB.BINARY, name='DummyP2')
    P2a = m.addVars(DummyP2a, vtype=GRB.BINARY, name='DummyP2a')
    P2b = m.addVars(DummyP2b, vtype=GRB.BINARY, name='DummyP2b')
    EP1 = m.addVars(DummyEP2, vtype=GRB.BINARY, name='DummyEP1')
    EP2 = m.addVars(DummyEP2, vtype=GRB.BINARY, name='DummyEP2')
    TP1 = m.addVars(DummyTP1, vtype=GRB.BINARY, name='DummyTP1')
    TP2 = m.addVars(DummyTP2, vtype=GRB.BINARY, name='DummyTP2')
    
    DummyZ1, DummyZ2 = [],[]
    for g in range(3,nG+1):
        for i in range(1, g+1): #x_i: anchor/del
            for j in range(1, g+1): #x_j: anchor/del
                for k in range(1, g+1):# z_k
                    if i != j and i != k and j != k:
                        DummyZ1.append((g,i,j,k))
                        DummyZ2.append((g,i,j,k))
    Dz1 = m.addVars(DummyZ1,vtype=GRB.BINARY, name='DummyZ1')
    Dz2 = m.addVars(DummyZ2,vtype=GRB.BINARY, name='DummyZ2')
    
    DummyP0 = []
    for prevg in range(2, nG):
        g = prevg+1
        for k in range(1,prevg+1):
            for i in range(1,g+1):
                for j in range(1,g+1):
                    if i != j:
                        DummyP0.append((g,k,i,j))
    S1term = m.addVars(DummyP0, vtype=GRB.BINARY, name='S1term')
    
    DummyP2, DummyP3, DummyP4 = [], [], []
    for prevg in range(2, nG):
        g = prevg+1
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:      
                    for l in range(1, prevg+1):
                        for k in range(1, prevg+1):
                            if l != k:
                                DummyP2.append((g,l,k,i,j))
                                DummyP3.append((g,l,k,i,j))
                                DummyP4.append((g,l,k,i,j))
    S2aterm = m.addVars(DummyP2, vtype=GRB.BINARY, name='S2aterm')
    S2bterm = m.addVars(DummyP3, vtype=GRB.BINARY, name='S2bterm')
    S3term = m.addVars(DummyP4, vtype=GRB.BINARY, name='S3term')
    
    print('All variables defined')
    
    #####################################################################################################
    # Set objective
    
    obj = LinExpr();
    for g in range(3,nG+1):
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:
                    obj += 0.5 * math.log(q_con) * Da[g,i,j] #0.5 bcs we count i,j twice
                    obj += 0.5 * math.log(1 - q_con) * Db[g,i,j]
        obj += sum(Ny[g,k]*math.log(1-q_mod) for k in range(1,g+1))
        obj += sum(Nz[g,l]*math.log(q_mod) for l in range(1,g+1))
    # + constraints for extant graph    
    m.setObjective(obj, GRB.MAXIMIZE)
    
    #####################################################################################################
    # Constraints:
    
    ###e_ij = e_ji
    for g in range(3,nG):
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:
                    ID = str(g)+'-'+str(i)+'-'+str(j)
                    #print('EdgeEquality'+ID)
                    m.addConstr(E[g,i,j] == E[g,j,i], 'EdgeEquality'+ID)
    
    ###############################################################################
    #print('Linearize Objective Function:')
    
    ###e_ij x_i x_j = Da_ij for each g
    for g in range(3,nG):
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:
                    m.addConstr(Da[g,i,j] <= E[g,i,j],'ObjDa1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Da[g,i,j] <= Nx[g,i],'ObjDa2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Da[g,i,j] <= Nx[g,j],'ObjDa3-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Da[g,i,j] >= E[g,i,j]+Nx[g,i]+Nx[g,j]-2,'ObjDa4-'+str(g)+'-'+str(i)+'-'+str(j))
                    #print('ObjDa',g,i,j)
    g = nG
    for i in range(1,g+1):
        for j in range(1,g+1):
            if i != j:
                m.addConstr(Da[g,i,j] == LG[i-1][j-1] * Dc[g,i,j],'ObjDa-LG-'+str(g)+'-'+str(i)+'-'+str(j))
                #print('ObjDa',g,i,j)
                
               
    ### (1-e_ij)x_ix_j = Db_ij   
    for g in range(3,nG+1):
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:
                    m.addConstr(Db[g,i,j] == Dc[g,i,j] - Da[g,i,j],'ObjDbc1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Dc[g,i,j] <= Nx[g,i],'ObjDc2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Dc[g,i,j] <= Nx[g,j],'ObjDc3-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(Dc[g,i,j] >= Nx[g,i]+Nx[g,j]-1,'ObjDc4-'+str(g)+'-'+str(i)+'-'+str(j))
                    #print('ObjDbc',g,i,j)
                    
    #print('--------------')
    
    ###############################################################################
                
    #print('Only 2 dup nodes per graph:')
    for g in range(3,nG+1):
        expr = LinExpr()
        for i in range(1,g+1):
            expr += Nx[g,i]
            #print(g,'d',i)
        m.addConstr(expr == 2, 'anchordup'+str(g))
        #print('anchordup'+str(g))
    #print('--------------')
    
    ##if q == 2, then y = 1, q = \sum
    ##0 <= y <= 1, y >= q-1, 2y <= q
    #print('Define Common Neighbor: 2 y_k <= \sum_i e_ik x_i; y_k >= \sum_i e_ik x_i -1')
    for g in range(3,nG):
        for k in range(1, g+1): #x_i: anchor/del
            expr = LinExpr()
            for i in range(1, g+1): #x_j: anchor/del
                if k != i:
                    expr += Dcn[g,i,k]
                    m.addConstr(Dcn[g,i,k] <= E[g,k,i], 'CommNeigh1-'+str(g)+'-'+str(i)+'-'+str(k))
                    m.addConstr(Dcn[g,i,k] <= Nx[g,i], 'CommNeigh2-'+str(g)+'-'+str(i)+'-'+str(k))
                    m.addConstr(Dcn[g,i,k] >= E[g,k,i]+Nx[g,i]-1, 'CommNeigh3-'+str(g)+'-'+str(i)+'-'+str(k))
            m.addConstr(expr >= 2*Ny[g,k], 'CommNeighI-'+str(g)+'-'+str(k))
            m.addConstr(expr <= Ny[g,k] + 1, 'CommNeighII-'+str(g)+'-'+str(k))
            #print('CommNeigh-'+str(g)+'-'+str(k))
            
    g = nG
    for k in range(1, g+1): #x_i: anchor/del
        expr = LinExpr()
        for i in range(1, g+1): #x_j: anchor/del
            if k != i:
                expr += LG[i-1][k-1]*Nx[g,i]
        m.addConstr(expr >= 2*Ny[g,k], 'CommNeighI-'+str(g)+'-'+str(k))
        m.addConstr(expr <= Ny[g,k] + 1, 'CommNeighII-'+str(g)+'-'+str(k))
    
        #print('CommNeigh-'+str(g)+'-'+str(k))
    
    
    #print('--------------')
    
    ##if q == 1 and x_k = 0, then z_k = 1, q = \sum
    ##0 <= z <= 1, z+2y = q
    #print('Define Neighbors in Symm Diff 2z_k + y_k = \sum_i e_ik x_i')
    
    for g in range(3,nG):
        for k in range(1, g+1): #x_i: anchor/del
            expr = LinExpr()
            for i in range(1, g+1): #x_j: anchor/del
                if k != i:
                    expr += Dcn[g,i,k]
            m.addConstr(expr == 2*Ny[g,k] + Nw[g,k], 'SymmDiffI-'+str(g)+'-'+str(k))
            m.addConstr(Nz[g,k] >= Nw[g,k] - Nx[g,k], 'SymmDiffII-'+str(g)+'-'+str(k))
            m.addConstr(Nz[g,k] <= Nw[g,k], 'SymmDiffIII-'+str(g)+'-'+str(k))
            m.addConstr(Nz[g,k] <= 1 - Nx[g,k], 'SymmDiffIV-'+str(g)+'-'+str(k))
            #print('SymmDiff-'+str(g)+'-'+str(k))
            
    g = nG
    for k in range(1, g+1): #x_i: anchor/del
        expr = LinExpr()
        for i in range(1, g+1): #x_j: anchor/del
            if k != i:
                expr += LG[i-1][k-1]*Nx[g,i]
        m.addConstr(expr == 2*Ny[g,k] + Nw[g,k], 'SymmDiffI-'+str(g)+'-'+str(k))
        m.addConstr(Nz[g,k] >= Nw[g,k] - Nx[g,k], 'SymmDiffII-'+str(g)+'-'+str(k))
        m.addConstr(Nz[g,k] <= Nw[g,k], 'SymmDiffIII-'+str(g)+'-'+str(k))
        m.addConstr(Nz[g,k] <= 1 - Nx[g,k], 'SymmDiffIV-'+str(g)+'-'+str(k))
        #print('SymmDiff-'+str(g)+'-'+str(k))
    
    
    #print('--------------')
    
    #print('Define Neigbor: n_k = y_k + z_k')
    for g in range(3,nG+1):
        for k in range(1,g+1):
            m.addConstr(Nn[g,k] == Ny[g,k]+Nz[g,k],'Neighbor-'+str(g)+'-'+str(k))
            #print('Neighbor-'+str(g)+'-'+str(k))
    
    ###############################################################################
    #Phantom Edges
    
    #print('At least one outgoing phantom edge:')
    for prevg in range(2, nG):
        g = prevg+1
        for i in range(1, prevg+1):
            expr = LinExpr()
            for j in range(1, g+1):
                #print(g,'<<out',i,j)
                expr += P[g,i,j]
            #print(g,'!')
            m.addConstr(expr >= 1, 'LBPhantomOutG'+str(g)+'-'+str(i))
            #print('LBPhantomOutG'+str(g)+'-'+str(i))
    #print('--------------')  
    
    #print('At most two outgoing phantom edge')
    for prevg in range(2, nG):
        g = prevg+1
        for i in range(1, prevg+1):
            expr = LinExpr()
            for j in range(1, g+1):
                #print(g,'>>out',i,j)
                expr += P[g,i,j]
            #print(g,'!')
            m.addConstr(expr <= 2, 'UBPhantomOutG'+str(g)+'-'+str(i))
            #print('UBPhantomOutG'+str(g)+'-'+str(i))
    #print('--------------')  
    
    #print('One incoming phantom edge:')
    for prevg in range(2, nG):
        g = prevg+1
        for j in range(1, g+1):
            expr = LinExpr()
            for i in range(1, prevg+1):
                #print(g,'<<in',i,j)
                expr += P[g,i,j]
            m.addConstr(expr == 1, 'PhantomInG'+str(g)+'-'+str(j))
            #print('PhantomInG'+str(g)+'-'+str(j))
    #print('--------------')  
    
    
    ###############################################################################
    
    #print('Only one anchor node per graph -- not in extant graph:')
    for g in range(2, nG):
        expr = LinExpr()
        for i in range(1,g+1):
            expr += Na[g,i]
            #print(g,'anchorNode',i)
        m.addConstr(expr == 1, 'anchor'+str(g))
        #print('anchor'+str(g))
    #print('--------------')  
    
    ###############################################################################
    #Phantom Edge Path
    #if x_ix_j = 1, then sum_k d_k P_ki P_kj = 1
    #elif x_n_j = 1, then  \sum_lk P_kj P_li d_l (1 - d_k) e_lk= 1
    # else 1 + \sum_lk P_kj P_li e_lk - e_ij = 1
    
    for prevg in range(3, nG-1):
        g = prevg+1
        for i in range(1,g+1):
            for j in range(1,g+1):
                if i != j:
                    
                    S1 = LinExpr()
                    for k in range(1,prevg+1):
                        S1 += S1term[g,k,i,j]  #P[g,k,i]*Na[prevg,k]*P[g,k,j]
                        m.addConstr(S1term[g,k,i,j] <= Na[prevg,k],'S1X1-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                        m.addConstr(S1term[g,k,i,j] <= P[g,k,i],'S1X2-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                        m.addConstr(S1term[g,k,i,j] <= P[g,k,j],'S1X3-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                        m.addConstr(S1term[g,k,i,j] >= Na[prevg,k]+P[g,k,i]+P[g,k,j]-2,'S1X4-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                        #print('S1X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                        
                    #l: anchor, i: x, j: neighbor, k, non-anchor
                    #P[g,l,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,k,j]
                    S2a = LinExpr()
                    for l in range(1, prevg+1):
                        for k in range(1, prevg+1):
                            if l != k:
                                S2a += S2aterm[g,l,k,i,j]  
                                m.addConstr(S2aterm[g,l,k,i,j] <= P[g,l,i],'S2aX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2aterm[g,l,k,i,j] <= P[g,k,j],'S2aX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2aterm[g,l,k,i,j] <= Na[prevg,l],'S2aX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2aterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2aX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2aterm[g,l,k,i,j] <= E[prevg,l,k],'S2aX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2aterm[g,l,k,i,j] >= P[g,l,i]+P[g,k,j]+Na[prevg,l]-Na[prevg,k]+E[prevg,l,k]-3,'S2aX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                
                    #l: anchor, i: neighbor, j: x, k, non-anchor
                    #P[g,k,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,l,j]
                    S2b = LinExpr()
                    for l in range(1, prevg+1):
                        for k in range(1, prevg+1):
                            if l != k:
                                S2b += S2bterm[g,l,k,i,j]  
                                m.addConstr(S2bterm[g,l,k,i,j] <= P[g,k,i],'S2bX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2bterm[g,l,k,i,j] <= P[g,l,j],'S2bX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2bterm[g,l,k,i,j] <= Na[prevg,l],'S2bX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2bterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2bX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2bterm[g,l,k,i,j] <= E[prevg,l,k],'S2bX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S2bterm[g,l,k,i,j] >= P[g,k,i]+P[g,l,j]+Na[prevg,l]-Na[prevg,k]+E[prevg,l,k]-3,'S2bX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                
                    T = LinExpr()
                    for l in range(1, prevg+1):
                        for k in range(1, prevg+1):
                            if l != k:
                                T += S3term[g,l,k,i,j]  #P[g,l,i]*E[prevg,l,k]*P[g,k,j] 
                                m.addConstr(S3term[g,l,k,i,j] <= P[g,l,i],'S3X1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S3term[g,l,k,i,j] <= E[prevg,l,k],'S3X2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S3term[g,l,k,i,j] <= P[g,k,j],'S3X3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                m.addConstr(S3term[g,l,k,i,j] >= P[g,l,i]+E[prevg,l,k]+P[g,k,j]-2,'S3X4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                                #print('S3X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                    #S3 = T + 1 - LG[i-1][j-1]
                                
                    #P1 = Dc[g,i,j]
                    m.addConstr(P2a[g,i,j] <= Nx[g,i],'P2aX1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2a[g,i,j] <= Nn[g,j],'P2aX2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2a[g,i,j] >= Nx[g,i]+Nn[g,j]-1,'P2aX3-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    m.addConstr(P2b[g,i,j] <= Nn[g,i],'P2bX1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2b[g,i,j] <= Nx[g,j],'P2bX2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2b[g,i,j] >= Nn[g,i]+Nx[g,j]-1,'P2bX3-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    m.addConstr(S2a == P2a[g,i,j],'S2a-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(S2b == P2b[g,i,j],'S2b-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    m.addConstr(S1 == Dc[g,i,j],'S1-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    #P2 = x_i n_j OR n_i x_j
                    m.addConstr(P2[g,i,j] <= P2a[g,i,j] + P2b[g,i,j],'P2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2[g,i,j] >= P2a[g,i,j],'P2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(P2[g,i,j] >= P2b[g,i,j],'P2X3-'+str(g)+'-'+str(i)+'-'+str(j))
        
                    m.addConstr(T >= P2[g,i,j],'TX1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(T <= 1 - Dc[g,i,j] + P2[g,i,j],'TX2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(T <= 1,'TUB-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(T >= 0,'TLB-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    #e_ij * (1-P1), e_ij*(1-P2)
                    m.addConstr(EP1[g,i,j] <= E[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(EP1[g,i,j] <= 1-Dc[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(EP1[g,i,j] >= E[g,i,j]-Dc[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    m.addConstr(EP2[g,i,j] <= E[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(EP2[g,i,j] <= 1-P2[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(EP2[g,i,j] >= E[g,i,j]-P2[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                   
                    #T*(1-P1), T*(1-P2)
                    m.addConstr(TP1[g,i,j] <= T,'TP1X1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(TP1[g,i,j] <= 1 - Dc[g,i,j],'TP1X2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(TP1[g,i,j] >= T - Dc[g,i,j],'TP1X3-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(TP2[g,i,j] <= T,'TP2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(TP2[g,i,j] <= 1 - P2[g,i,j],'TP2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(TP2[g,i,j] >= T - P2[g,i,j],'TP2X3-'+str(g)+'-'+str(i)+'-'+str(j))
                    
                    m.addConstr(EP1[g,i,j] <= TP1[g,i,j],'EPTP1-'+str(g)+'-'+str(i)+'-'+str(j))
                    m.addConstr(EP2[g,i,j] >= TP2[g,i,j],'EPTP2-'+str(g)+'-'+str(i)+'-'+str(j))
                    
    prevg, g = nG-1,nG
    for i in range(1,g+1):
        for j in range(1,g+1):
            if i != j:
                S1 = LinExpr()
                for k in range(1,prevg+1):
                    S1 += S1term[g,k,i,j]  #P[g,k,i]*Na[prevg,k]*P[g,k,j]
                    m.addConstr(S1term[g,k,i,j] <= Na[prevg,k],'S1X1-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] <= P[g,k,i],'S1X2-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] <= P[g,k,j],'S1X3-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] >= Na[prevg,k]+P[g,k,i]+P[g,k,j]-2,'S1X4-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    #print('S1X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                    
                #l: anchor, i: x, j: neighbor, k, non-anchor
                #P[g,l,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,k,j]
                S2a = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            S2a += S2aterm[g,l,k,i,j]  
                            m.addConstr(S2aterm[g,l,k,i,j] <= P[g,l,i],'S2aX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= P[g,k,j],'S2aX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= Na[prevg,l],'S2aX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2aX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= E[prevg,l,k],'S2aX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] >= P[g,l,i]+P[g,k,j]+Na[prevg,l]-Na[prevg,k]+E[prevg,l,k]-3,'S2aX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            
                #l: anchor, i: neighbor, j: x, k, non-anchor
                #P[g,k,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,l,j]
                S2b = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            S2b += S2bterm[g,l,k,i,j]  
                            m.addConstr(S2bterm[g,l,k,i,j] <= P[g,k,i],'S2bX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= P[g,l,j],'S2bX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= Na[prevg,l],'S2bX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2bX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= E[prevg,l,k],'S2bX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] >= P[g,k,i]+P[g,l,j]+Na[prevg,l]-Na[prevg,k]+E[prevg,l,k]-3,'S2bX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            
                T = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            T += S3term[g,l,k,i,j]  #P[g,l,i]*E[prevg,l,k]*P[g,k,j] 
                            m.addConstr(S3term[g,l,k,i,j] <= P[g,l,i],'S3X1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S3term[g,l,k,i,j] <= E[prevg,l,k],'S3X2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S3term[g,l,k,i,j] <= P[g,k,j],'S3X3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S3term[g,l,k,i,j] >= P[g,l,i]+E[prevg,l,k]+P[g,k,j]-2,'S3X4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S3X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                #S3 = T + 1 - LG[i-1][j-1]
                            
                #P1 = Dc[g,i,j]
                #x_i n_j
                m.addConstr(P2a[g,i,j] <= Nx[g,i],'P2aX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2a[g,i,j] <= Nn[g,j],'P2aX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2a[g,i,j] >= Nx[g,i]+Nn[g,j]-1,'P2aX3-'+str(g)+'-'+str(i)+'-'+str(j))
                #n_i x_j
                m.addConstr(P2b[g,i,j] <= Nn[g,i],'P2bX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2b[g,i,j] <= Nx[g,j],'P2bX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2b[g,i,j] >= Nn[g,i]+Nx[g,j]-1,'P2bX3-'+str(g)+'-'+str(i)+'-'+str(j))
                #S2 = P2 constraints
                m.addConstr(S2a == P2a[g,i,j],'S2a-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(S2b == P2b[g,i,j],'S2b-'+str(g)+'-'+str(i)+'-'+str(j))
                #S1 = P1
                m.addConstr(S1 == Dc[g,i,j],'S1-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #P2 = x_i n_j OR n_i x_j
                m.addConstr(P2[g,i,j] <= P2a[g,i,j] + P2b[g,i,j],'P2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2[g,i,j] >= P2a[g,i,j],'P2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2[g,i,j] >= P2b[g,i,j],'P2X3-'+str(g)+'-'+str(i)+'-'+str(j))
                #T >= P2, T <= 1-P1
                m.addConstr(T >= P2[g,i,j],'TX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T <= 1 - Dc[g,i,j] + P2[g,i,j],'TX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T <= 1,'TUB-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T >= 0,'TLB-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #e_ij * (1-P1), e_ij*(1-P2)
                m.addConstr(EP1[g,i,j] == LG[i-1][j-1] - LG[i-1][j-1]*Dc[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP2[g,i,j] == LG[i-1][j-1] - LG[i-1][j-1]*P2[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #T*(1-P1), T*(1-P2)
                m.addConstr(TP1[g,i,j] <= T,'TP1X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP1[g,i,j] <= 1 - Dc[g,i,j],'TP1X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP1[g,i,j] >= T - Dc[g,i,j],'TP1X3-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] <= T,'TP2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] <= 1 - P2[g,i,j],'TP2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] >= T - P2[g,i,j],'TP2X3-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #e_ij(1-P1) <= T(1-P1)
                #e_ij(1-P2) >= T(1-P2)
                m.addConstr(EP1[g,i,j] <= TP1[g,i,j],'EPTP1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP2[g,i,j] >= TP2[g,i,j],'EPTP2-'+str(g)+'-'+str(i)+'-'+str(j))
    
                
    prevg, g = 2,3
    for i in range(1,g+1):
        for j in range(1,g+1):
            if i != j:
                S1 = LinExpr()
                for k in range(1,prevg+1):
                    S1 += S1term[g,k,i,j]  #P[g,k,i]*Na[prevg,k]*P[g,k,j]
                    m.addConstr(S1term[g,k,i,j] <= Na[prevg,k],'S1X1-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] <= P[g,k,i],'S1X2-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] <= P[g,k,j],'S1X3-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    m.addConstr(S1term[g,k,i,j] >= Na[prevg,k]+P[g,k,i]+P[g,k,j]-2,'S1X4-'+str(g)+'-'+str(i)+'-'+str(j)+'-'+str(k))
                    #print('S1X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                    
                #l: anchor, i: x, j: neighbor, k, non-anchor
                #P[g,l,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,k,j]
                S2a = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            S2a += FG[l-1][k-1] * S2aterm[g,l,k,i,j]  
                            m.addConstr(S2aterm[g,l,k,i,j] <= P[g,l,i],'S2aX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= P[g,k,j],'S2aX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= Na[prevg,l],'S2aX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2aX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #m.addConstr(S2aterm[g,l,k,i,j] <= E[prevg,l,k],'S2aX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2aterm[g,l,k,i,j] >= P[g,l,i]+P[g,k,j]+Na[prevg,l]-Na[prevg,k]-2,'S2aX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                           
                            
                #l: anchor, i: neighbor, j: x, k, non-anchor
                #P[g,k,i]*Na[prevg,l]*(1-Na[prevg,k)*E[prevg,l,k]*P[g,l,j]
                S2b = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            S2b += FG[l-1][k-1] * S2bterm[g,l,k,i,j]  
                            m.addConstr(S2bterm[g,l,k,i,j] <= P[g,k,i],'S2bX1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= P[g,l,j],'S2bX2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= Na[prevg,l],'S2bX3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] <= 1 - Na[prevg,k],'S2bX4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #m.addConstr(S2bterm[g,l,k,i,j] <= E[prevg,l,k],'S2bX5-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S2bterm[g,l,k,i,j] >= P[g,k,i]+P[g,l,j]+Na[prevg,l]-Na[prevg,k]-2,'S2bX6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S2X1-6-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            
                T = LinExpr()
                for l in range(1, prevg+1):
                    for k in range(1, prevg+1):
                        if l != k:
                            T += FG[l-1][k-1] * S3term[g,l,k,i,j]  #P[g,l,i]*E[prevg,l,k]*P[g,k,j] 
                            m.addConstr(S3term[g,l,k,i,j] <= P[g,l,i],'S3X1-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #m.addConstr(S3term[g,l,k,i,j] <= E[prevg,l,k],'S3X2-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S3term[g,l,k,i,j] <= P[g,k,j],'S3X3-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            m.addConstr(S3term[g,l,k,i,j] >= P[g,l,i]+P[g,k,j]-1,'S3X4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                            #print('S3X1-4-'+str(g)+'-'+str(l)+'-'+str(k)+'-'+str(i)+'-'+str(j))
                #S3 = T + 1 - LG[i-1][j-1]
                            
                #P1 = Dc[g,i,j]
                m.addConstr(P2a[g,i,j] <= Nx[g,i],'P2aX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2a[g,i,j] <= Nn[g,j],'P2aX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2a[g,i,j] >= Nx[g,i]+Nn[g,j]-1,'P2aX3-'+str(g)+'-'+str(i)+'-'+str(j))
                
                m.addConstr(P2b[g,i,j] <= Nn[g,i],'P2bX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2b[g,i,j] <= Nx[g,j],'P2bX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2b[g,i,j] >= Nn[g,i]+Nx[g,j]-1,'P2bX3-'+str(g)+'-'+str(i)+'-'+str(j))
                
                m.addConstr(S2a == P2a[g,i,j],'S2a-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(S2b == P2b[g,i,j],'S2b-'+str(g)+'-'+str(i)+'-'+str(j))
                
                m.addConstr(S1 == Dc[g,i,j],'S1-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #P2 = x_i n_j OR n_i x_j
                m.addConstr(P2[g,i,j] <= P2a[g,i,j] + P2b[g,i,j],'P2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2[g,i,j] >= P2a[g,i,j],'P2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(P2[g,i,j] >= P2b[g,i,j],'P2X3-'+str(g)+'-'+str(i)+'-'+str(j))
    
                m.addConstr(T >= P2[g,i,j],'TX1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T <= 1 - Dc[g,i,j] + P2[g,i,j],'TX2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T <= 1,'TUB-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(T >= 0,'TLB-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #e_ij * (1-P1), e_ij*(1-P2)
                m.addConstr(EP1[g,i,j] <= E[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP1[g,i,j] <= 1-Dc[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP1[g,i,j] >= E[g,i,j]-Dc[g,i,j],'EP1-'+str(g)+'-'+str(i)+'-'+str(j))
                
                m.addConstr(EP2[g,i,j] <= E[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP2[g,i,j] <= 1-P2[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP2[g,i,j] >= E[g,i,j]-P2[g,i,j],'EP2-'+str(g)+'-'+str(i)+'-'+str(j))
                
                #T*(1-P1), T*(1-P2)
                m.addConstr(TP1[g,i,j] <= T,'TP1X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP1[g,i,j] <= 1 - Dc[g,i,j],'TP1X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP1[g,i,j] >= T - Dc[g,i,j],'TP1X3-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] <= T,'TP2X1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] <= 1 - P2[g,i,j],'TP2X2-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(TP2[g,i,j] >= T - P2[g,i,j],'TP2X3-'+str(g)+'-'+str(i)+'-'+str(j))
                
                m.addConstr(EP1[g,i,j] <= TP1[g,i,j],'EPTP1-'+str(g)+'-'+str(i)+'-'+str(j))
                m.addConstr(EP2[g,i,j] >= TP2[g,i,j],'EPTP2-'+str(g)+'-'+str(i)+'-'+str(j))
                    
    #print('--------------') 
    
    #####################################################################################################
    ##Initial Guess: from greedy
    ar_ind = 0
    pstart = []
    
    nx_list = get_all_networks(extant, anchor_rem_list)
    for g in range(nG,2,-1):
        a,r = anchor_rem_list[ar_ind]
        curr_nx = nx_list[ar_ind]
    #    print('\ng:',g,'curr_nx:',curr_nx.edges())
    #    print('AR:',a,r)
        
        ND, rND = {}, {} #node_dict
        if g != nG:
            ND[a], ND[r] = 1,2
            rND[1], rND[2] = a,r
            nxt = 3
            for x in curr_nx.nodes():
                if x != a and x != r:
                    ND[x] = nxt
                    rND[nxt] = x
                    nxt += 1     
            Nx[g,1].start = 1
            Nx[g,2].start = 1
            for i in range(3,g+1):
                Nx[g,i].start = 0
        else:
            Nx[g,a].start = 1
            Nx[g,r].start = 1
            extantA, extantR = a,r
            for i in range(1,g+1):
                ND[i], rND[i] = i,i
                if i != a and i != r:
                    Nx[g,i].start = 0
    #    print('ND', ND)
        
        NeighA = set(curr_nx.neighbors(a))
        NeighR = set(curr_nx.neighbors(r))
        CommonN = NeighA & NeighR
        SingleN = (NeighA ^ NeighR) - set([a,r])#(NeighA | NeighR) - CommonN
    #    print('commonN', CommonN)
    #    print('singleN', SingleN)  
        for i in range(1,g+1):
            if rND[i] in CommonN:
                #print('Ny set to 1', i)
                Ny[g,i].start = 1
                Nn[g,i].start = 1
            elif rND[i] in SingleN:
                #print('Nz set to 1', i)
                Nz[g,i].start = 1
                Nn[g,i].start = 1
            else:
                #print('Ny, Nz set to 0', i)
                Ny[g,i].start = 0
                Nz[g,i].start = 0
                Nn[g,i].start = 0
    
        if g != nG:
            edges = curr_nx.edges()
            for i in range(1,g+1):
                for j in range(i+1,g+1):
                    if i != j:
                        #print(i,j,edges)
                        if (rND[i],rND[j]) in edges or (rND[j],rND[i]) in edges:
                            E[g,i,j].start = 1
                            E[g,j,i].start = 1
                        else:
                            E[g,i,j].start = 0
                            E[g,j,i].start = 0
            #print(curr_nx.node, 'futA',futA)               
            for x in curr_nx.nodes():
                if x != futA:
                    #print('adding to pstart', x, ND[x], futND[x])
                    pstart.append((g+1,ND[x],futND[x]))
            if g+1 == nG:
                pstart.extend([(g+1,ND[futA],extantA),(g+1,ND[futA],extantR)])
            else:
                pstart.extend([(g+1,ND[futA],1),(g+1,ND[futA],2)])
            
        ar_ind += 1
        futND = ND
        futrND = rND
        futA = a
    
    a,r = anchor_rem_list[ar_ind]
    ND[a], ND[r] = 1,2  
    pstart.extend([(3,ND[futA],1),(3,ND[futA],2)])
    if ND[futA] == 2:
        pstart.append((3,1,3))
    else:
        pstart.append((3,2,3))
    #print('pstart', pstart,'\n') 
    
    #for g in range(2,nG):
    #    Na[g,1].start = 1
    #    for i in range(2,g+1):
    #        Na[g,i].start = 0
    ##x_ig
    #for g in range(3,nG+1):
    #    Nx[g,1].start = 1
    #    Nx[g,2].start = 1
    #    for i in range(3, g+1):
    #        Nx[g,i].start = 0
    #        
    #Ny[4,1].start = 0
    #Ny[4,2].start = 0
    #Ny[4,3].start = 0
    #Ny[4,4].start = 0
    #Nz[4,1].start = 0
    #Nz[4,2].start = 0
    #Nz[4,3].start = 1
    #Nz[4,4].start = 1
    #Ny[3,1].start = 0
    #Ny[3,2].start = 0
    #Ny[3,3].start = 0
    #Nz[3,1].start = 0
    #Nz[3,2].start = 0
    #Nz[3,3].start = 1
    #                
    #E[3,1,2].start = 1
    #E[3,1,3].start = 1
    #E[3,2,3].start = 0
    #E[3,2,1].start = 1
    #E[3,3,1].start = 1
    #E[3,3,2].start = 0
    ##phantom edges
    #pstart = [(3,1,1),(3,1,2),(3,2,3),(4,1,1),(4,1,2),(4,2,3),(4,3,4)]
    for prevg in range(2, nG):
        g = prevg+1
        for i in range(1, prevg+1):
            for j in range(1,g+1):
                if (g,i,j) in pstart:
                    #print(g,i,j)
                    P[g,i,j].start = 1
                else:
                    P[g,i,j].start = 0
    #####################################################################################################
    # Optimize
    
    m.optimize()
    status = m.status
    
    
    if status == GRB.Status.OPTIMAL or status == GRB.TIME_LIMIT:
    #    print('The optimal objective is %g' % m.objVal)
    #    for v in m.getVars():
    #        if (v.varName.startswith('Edge') or v.varName.startswith('NodeX') or v.varName.startswith('NodeN') or v.varName.startswith('Phantom') )and v.x > 0:
    #            print('%s %g' % (v.varName, v.x))
    ###############################################################################
#        Xnodes = []
#        for g in range(3,nG):
#            prevg = g-1
#            nodelist = []
#            for i in range(1,g+1):
#                namestr = 'NodeX['+str(g)+','+str(i)+']'
#                val = m.getVarByName(namestr).x
#                if val != 0:
#                    #print(namestr, val)
#                    nodelist.append(i)
#            Xnodes.append(nodelist)
#         
#        prevg, g = nG-1, nG
#        nodelist = []
#        for i in range(1,g+1):
#            namestr = 'NodeX['+str(g)+','+str(i)+']'
#            val = m.getVarByName(namestr).x
#            if val != 0:
#                #print(namestr, val)
#                nodelist.append(i)
#        Xnodes.append(nodelist)
#        print('NodeX',Xnodes)

        #print number of solutions
        nSolutions = m.SolCount
        print('#Solutions: %s'% str(nSolutions))        

        # Print objective values of solutions
        for e in range(nSolutions): 
            m.setParam(GRB.Param.SolutionNumber, e) 
            print('%g\n'% m.PoolObjVal)
        
            with open(outputFNAME+'_solution'+str(e)+'.txt','w') as f:
                f.write('Obj: %g\n' % obj.getValue())                
                for v in m.getVars():
                    f.write('%s %g\n' % (v.varName, v.xn))
            
        
            print('Obj: %g' % obj.getValue())
    ###############################################################################
    
    if status == GRB.Status.UNBOUNDED:
        print('The model cannot be solved because it is unbounded')
        
    if status != GRB.Status.INF_OR_UNBD and status != GRB.Status.INFEASIBLE:
        print('Optimization was stopped with status %d' % status)
        
    if status == GRB.Status.INF_OR_UNBD or status == GRB.Status.INFEASIBLE:
        print('Optimization was stopped with status %d' % status)
        m.computeIIS();
        for c in m.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
        m.write('m.mps')
    #####################################################################################################
        
    #    #Relax the constraints to make the model feasible
    #    print('The model is infeasible; relaxing the constraints')
    #    orignumvars = m.NumVars
    #    m.feasRelaxS(0, False, False, True)
    #    m.optimize()
    #    status = m.status
    #    if status in (GRB.Status.INF_OR_UNBD, GRB.Status.INFEASIBLE, GRB.Status.UNBOUNDED):
    #        print('The relaxed model cannot be solved \
    #               because it is infeasible or unbounded')
    #        exit(1)
    #    
    #    if status != GRB.Status.OPTIMAL:
    #        print('Optimization was stopped with status %d' % status)
    #        exit(1)
    #    
    #    print('\nSlack values:')
    #    slacks = m.getVars()[orignumvars:]
    #    for sv in slacks:
    #        if sv.X > 1e-6:
    #            print('%s = %g' % (sv.VarName, sv.X))
    #    
    #    print('\nVariableValues:')        
    #    for v in m.getVars():
    #        print('%s %g' % (v.varName, v.x))
    #    print('Obj: %g' % obj.getValue())

parser = argparse.ArgumentParser(description="ILP to reconstruct maximum likelihood network evolution history using Duplication Mutation with Complementarity (DMC) model. Uses Gurobi solver, see http://www.gurobi.com/documentation/8.0/refman/index.html")

parser.add_argument("-e","--extant",type=str,required=True,help="Extant filename. Format: one edge per line")
parser.add_argument("-g","--greedy",type=str,required=True,help="Filename of Solution from Greedy Approach ReverseDMC")
parser.add_argument("-o","--output",type=str,required=True,help="Output Filename")
parser.add_argument("-c","--q_con",type=float,default=0.7,help="DMC Model Parameter q_con (default: %(default)s)")
parser.add_argument("-m","--q_mod",type=float,default=0.4,help="DMC Model Parameter q_mod (default: %(default)s)")
parser.add_argument("-t","--timelimit",type=int,default=24,help="Time Limit in hours (default: %(default)s)")
parser.add_argument("-n","--numcores",type=int,default=0,help="Number of cores. 0 uses all available cores (default: %(default)s)")
parser.add_argument("-s","--numsolutions",type=int,default=30,help="Number of ILP Solutions (default: %(default)s)")
parser.add_argument("-p","--poolmode",type=int,default=1,help="Pool mode, possible values 0,1,2 (default: %(default)s). 0: finds one optimal solution. 1: find multiple solutions not necessarily the best. 2: find n best multiple solutions.")
parser.add_argument("-f","--focus",type=int,default=1,help="MIP focus, possible values 1,2,3 (default: %(default)s). 1: finds feasible solutions quickly. 2: to prove optimality, if good quality solutions can be found easily. 3: if the best objective bound is moving very slowly or not at all.")
parser.add_argument("-r","--heuristics",type=float,default=0.5,help="Fraction of time spent on MIP heuristics, between 0 and 1 (default: %(default)s). Larger values produce more and better feasible solutions, at a cost of slower progress in the best bound.")

print("Program arguments:")
args,_= parser.parse_known_args()
for arg in vars(args):
    print(arg, getattr(args, arg))


extantFNAME = args.extant 
greedyFNAME = args.greedy 
outputFNAME = args.output 
q_con = args.q_con
q_mod = args.q_mod
timeLimit = args.timelimit
numThreads = args.numcores
numSolutions = args.numsolutions
poolMode = args.poolmode
focus = args.focus
heuristicstime = args.heuristics

#run_ILP(extantFNAME, q_con, q_mod, greedyFNAME, outputFNAME, timeLimit, numThreads)
run_ILP(extantFNAME, q_con, q_mod, greedyFNAME, outputFNAME, timeLimit, numThreads, numSolutions, poolMode, focus, heuristicstime)
"""
if len(sys.argv) != 6:
   print('Usage: $python ILP.py <q_con> <q_mod> <extant filename> <greedy filename> <output filename>\n')
   sys.exit(-1)

numHours = 21
timeLimit = numHours*60.0*60.0 #in seconds
numThreads = 50 #default 0: uses all cores available

q_con = float(sys.argv[1])
q_mod = float(sys.argv[2])
#q_con = 0.7
#q_mod = 0.4
extantFNAME = sys.argv[3] #'extant'
greedyFNAME = sys.argv[4] #'greedy_extant_qmod=0.4_qcon=0.7'
outputFNAME = sys.argv[5] #'ILP_'+extantFNAME+'_qmod='+str(q_mod)+'_qcon='+str(q_con)
run_ILP(extantFNAME, q_con, q_mod, greedyFNAME, outputFNAME, timeLimit, numThreads)
"""

