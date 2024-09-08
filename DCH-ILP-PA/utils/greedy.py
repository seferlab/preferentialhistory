import numpy as np
import networkx as nx
import json
import random
from utils.get_anchor_list import *

def dmc_likelihood(G, u, v, qmod, qcon):
    """ Computes the likelihood of u merging with v according to the DMC model. """
    # Ignore edges to each other because those will be taken care of by qcon.
    U = set([n for n in G.neighbors(u)]) - set([v])
    V = set([n for n in G.neighbors(v)]) - set([u])

    Cuv = len(U & V)
    Suv = len(U) + len(V) - 2*Cuv
    gamma = qcon if G.has_edge(u,v) else 1-qcon

    return (1-qmod)**(Cuv) * (qmod/2)**(Suv) * gamma



#==============================================================================
#                           RECONSTRUCTION FUNCTIONS
#==============================================================================
def dmc_delorean(g,qmod,qcon):
    """ Reconstructs the network using the dmc model and delorean algorithm. """
    G = g.copy()
    ARlist = []
    # Make initial pairwise likelihoods.
    L = {}
    for u in G.nodes(): 
        L[u] = {}

    for u in G.nodes():
        for v in G.nodes():
            if u >= v: continue
            L[u][v] = dmc_likelihood(G, u,v,qmod,qcon)


    while (len(G.nodes()) >= 2): # at least two nodes in the graph.

        # Get largest Luv.
        L_list = []
        L_prob = -10000000000

        for u in G:
            for v in G:
                if u >= v: continue

                Luv = L[u][v]
                if Luv > L_prob:
                    L_list = [(u,v)]
                    L_prob = Luv
                elif Luv == L_prob:
                    L_list.append((u,v))

        # Choose random pair; assign random daddy.
        pair = random.choice(L_list)
        (u,v) = (pair[0],pair[1]) if random.random() > 0.5 else (pair[1],pair[0])

        # Nodes whose likelihood values need to be computed.
        altered = (set([n for n in G.neighbors(u)]) | set([n for n in G.neighbors(v)]) | set([u])) - set([v])

        # Prepare to delete v: add new edges in symmetric difference of v to u.
        for neighbor in G.neighbors(v):
            if u == neighbor: continue # Don't add self-edge.
            elif v == neighbor: continue # Don't add, will remove v anyways.
            elif G.has_edge(u,neighbor): continue # Edge already exists.
            else: G.add_edge(u,neighbor)
        G.remove_node(v)

        # append the anchor duplicate pair
        ARlist.append((u,v))

        # Fix up altered Luv values.
        for x in altered:
            for y in G.nodes():
                if x == y: continue
                L[min(x,y)][max(x,y)] = dmc_likelihood(G, x,y,qmod,qcon)

    last_node = [n for n in G.nodes()]
    return ARlist, last_node   

