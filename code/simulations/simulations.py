# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 09:18:26 2017

@author: disvr
"""
import networkx as nx
import random
import pathlib

def forward_DMC(G, q_mod, q_con):
    assert(q_mod < 1 and q_mod > 0)
    assert(q_con < 1 and q_con > 0)
    #select anchor node
    nodes = G.nodes()
    newnode = sorted(nodes)[-1] + 1
    while(True):
        #duplicate node and edges
        H = G.copy()
        H.add_node(newnode)

        anchor = random.choice(nodes)
        #print('chosen anchor', anchor)
        for x in G.neighbors(anchor):
            H.add_edge(newnode, x)
        #remove edges
        for x in H.neighbors(newnode):
            one = random.choice([newnode, anchor])
            if random.random() < q_mod:
                H.remove_edge(one, x)
                #print('deleted', one, x)
        #add edge to anchor
        if random.random() < q_con:
            H.add_edge(anchor, newnode)
            #print('added', anchor, newnode)
        ##############KEEPING ONLY CONNECTED GRAPHS CHANGE LATER#########################
        if (nx.is_connected(H)):
            break
    #print('final edges:', H.edges())
    return H

def simulate_evolution(num_nodes, q_mod, q_con):    
    G1=nx.Graph()
    G1.add_node(1)
    G1.add_node(2)
    G1.add_edge(1,2)
    graph_list = []
    #print('G #1 created')
    #print(G1.edges())
    G = G1
    for i in range(3, num_nodes+1):
        #print('Creating G#', i-1)
        G = forward_DMC(G, q_mod, q_con)
        graph_list.append(G)
        
    return graph_list

def write_graph(filename, graph):
    with open(filename, 'w') as f:
        for x,y in graph.edges():
            f.write(str(x)+' '+str(y)+'\n')
            
            
probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]    
nG = 10

for run in range(1,6):
    folder = 'run'+str(run)
    pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    for q_con in probL:
        for q_mod in probL:
            graphL = simulate_evolution(nG, q_mod, q_con)
            assert(len(graphL) == nG-2)
            for G in graphL:
                n = len(G.nodes())
                name = 'nx_run='+str(run)+'_qmod='+str(q_mod)+'_qcon='+str(q_con)+'_n='+str(n)
                name = folder +'/'+ name
                print(name)
                write_graph(name,G)
            
        