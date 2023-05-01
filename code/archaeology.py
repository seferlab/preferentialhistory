#!/usr/bin/env python

"""
    Implementation of the Network Archaeology algorithms

    @requires: python 2.5
"""

from __future__ import division
from optparse import OptionParser
from collections import deque
from math import log
import sys,time,logging,random

random.seed(10301949)


#==============================================================================
#                               NODE OBJECT
#==============================================================================
class Node:
    def __init__(self,id):
        self.id = id
        self.nodes = [id]
        self.neighbors = set()

    def __cmp__(self,other):
        if self.id > other.id: return 1
        elif self.id < other.id: return -1
        else: return 0

    def __str__(self):
        s =  "%s --> " %(self.id)
        for index,neighbor in enumerate(self.neighbors):
            s += "(%i) " %(neighbor) # self.weights[index]
        s += "\n"
        return s


#==============================================================================
#                               GRAPH OBJECT
#==============================================================================
class Graph:
    def __init__(self,name):
        self.node_map = {} # Mapping from node id#'s to Node objects.
        self.num_nodes = 0
        self.num_edges = 0
        self.name = name

    def __iter__(self):
        return iter(self.node_map)

    def __str__(self):
        s = "%s: #nodes: %s, #edges: %s\n\n" %(self.name,self.num_nodes,self.num_edges)
        for i in self.node_map.iterkeys():
            s += "%s" %(self.node_map[i])
        return s

    @staticmethod
    def read_graph(filename):
        """ Creates a graph object from an edgelist file. """
        G = Graph("")
        input = open(filename)
        for line in input:
            [u,v] = line.strip().split()
            if u == v: continue
            if not G.has_node(u): G.add_node(u)
            if not G.has_node(v): G.add_node(v)
            G.add_edge(u,v)

        logging.info("Graph: #nodes=%i, #edges=%i" %(G.num_nodes,G.num_edges))
        return G

    def add_node(self,i):
        """ Adds a given node to the graph, as an isolate. """
        self.num_nodes += 1
        self.node_map[i] = Node(i)

    def remove_node(self,i):
        """ Deletes the given node from the graph. """
        I = self.node_map[i]

        # Remove i from the neighbor lists of all of i's neighbors.
        for neighbor in I.neighbors:
            if i != neighbor:
                self.node_map[neighbor].neighbors.remove(i)
            self.num_edges -= 1

        del self.node_map[i]
        self.num_nodes -= 1

    def has_node(self,u):
        """ Returns True if the given node is in the graph; False otherwise. """
        return u in self.node_map

    def add_edge(self,u,v):
        """ Adds edge to the graph. """
        self.node_map[u].neighbors.add(v)
        if u != v: self.node_map[v].neighbors.add(u)
        self.num_edges += 1

    def delete_edge(self,u,v):
        """ Deletes edge from the graph. """
        self.node_map[u].neighbors.remove(v)
        if u != v: self.node_map[v].neighbors.remove(u)
        self.num_edges -= 1

    def has_edge(self,u,v):
        """ Returns True if the edge exists in the graph; False otherwise. """
        return v in self.node_map[u].neighbors

    def nodes(self):
        """ Return all keys in the graph object. """
        return self.node_map.keys()

    def random_node(self):
        """ Returns a random node id. """
        # TODO: something better than populating keys() each time?
        return random.choice(self.node_map.keys())

    def neighbors(self,u):
        """ Returns set of neighbors of u. """
        return self.node_map[u].neighbors

    def dmc_likelihood(self,u,v,qmod,qcon):
        """ Computes the likelihood of u merging with v according to the DMC model. """
        # Ignore edges to each other because those will be taken care of by qcon.
        U = self.node_map[u].neighbors - set([v])
        V = self.node_map[v].neighbors - set([u])

        Cuv = len(U & V)
        Suv = len(U) + len(V) - 2*Cuv
        gamma = qcon if self.has_edge(u,v) else 1-qcon

        return (1-qmod)**(Cuv) * (qmod/2)**(Suv) * gamma



#==============================================================================
#                           RECONSTRUCTION FUNCTIONS
#==============================================================================
def dmc_delorean(G,qmod,qcon):
    """ Reconstructs the network using the dmc model and delorean algorithm. """

    # Make initial pairwise likelihoods.
    L = {}
    for u in G: L[u] = {}

    for u in G:
        for v in G:
            if u >= v: continue
            L[u][v] = G.dmc_likelihood(u,v,qmod,qcon)


    while (G.num_nodes >= 2): # at least two nodes in the graph.

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
        altered = (G.neighbors(u) | G.neighbors(v) | set([u])) - set([v])

        # Prepare to delete v: add new edges in symmetric difference of v to u.
        for neighbor in G.neighbors(v):
            if u == neighbor: continue # Don't add self-edge.
            elif v == neighbor: continue # Don't add, will remove v anyways.
            elif G.has_edge(u,neighbor): continue # Edge already exists.
            else: G.add_edge(u,neighbor)
        G.remove_node(v)

        print "%s\t%s" %(u,v)

        # Fix up altered Luv values.
        for x in altered:
            for y in G:
                if x == y: continue
                L[min(x,y)][max(x,y)] = G.dmc_likelihood(x,y,qmod,qcon)

    last_node = G.nodes()[0]
    print "%s\t%s" %(last_node,last_node)



def ff_delorean(G,p):
    """ Reconstructs the network using the ff model and delorean algorithm. """

    def geometric(p):
        """ Prob wrt. geometric distribution. """
        return int(round(log(random.random()) / log(1-p)))

    NUM_SIMULATIONS = 10

    while (G.num_nodes >= 2):
        L_list = []
        L_overlap = -1
        for u in G:

            u_neighbors = G.neighbors(u)
            for v in u_neighbors: # candidate anchors.
                assert v != u

                # run the ff simulation -- u when anchored from v.
                overlap = 0 # number of times the visited set is exactly the neighbor set.
                for sims in xrange(NUM_SIMULATIONS):
                    queue = deque([v])
                    visited = set([u,v]) # u wasn't in there originally, so can't be visited.

                    while len(queue) > 0:
                        w = queue.popleft()

                        # Flip geometric coin with probability 1-p and choose that many neighbors.
                        x = geometric(1-p)

                        # neighbors, without u.
                        neighbors = G.neighbors(w) - set([u]) # neighbors, without u (not there).
                        spread_to = random.sample(neighbors,x) if len(neighbors) > x else neighbors

                        for candidate in spread_to:
                            if candidate in visited: continue
                            visited.add(candidate)
                            queue.append(candidate)

                    visited.remove(u) # only for simulation. u is not a neighbor of itself.
                    if len(visited ^ u_neighbors) == 0: overlap += 1

                if overlap > L_overlap:
                    L_list = [(u,v)]
                    L_overlap = overlap
                elif overlap == L_overlap:
                    L_list.append((u,v))

        if len(L_list) == 0: # choose pairs randomly.
            u = v = 1
            while u == v:
                u = G.random_node()
                v = G.random_node()
        else:
            pair = random.choice(L_list) # uniformly choose from the best.
            (u,v) = (pair[0],pair[1])

        G.remove_node(u)
        print "%s\t%s" %(v,u)


    last_node = G.nodes()[0]
    print "%s\t%s" %(last_node,last_node)


def pa_delorean(G):
    """ Reconstructs the network using the pa model and delorean algorithm. """

    IMPOSSIBLE = -100000

    while G.num_nodes >= 2:

        # Find the nodes with the lowest degree.
        L_list = []
        L_degree = 10000
        for u in G:
            u_degree = len(G.neighbors(u))
            if u_degree < L_degree:
                L_list = [u]
                L_degree = u_degree
            elif u_degree == L_degree:
                L_list.append(u)


        current_sum = -100
        current_list = []
        for u in L_list:

            # Calculates: \prod_{v \in N(u)} d_v \prod_{v \not\in N(u)} m - d_v
            sum_degrees_of_neighbors = 0
            num_edges = 2*(G.num_edges - len(G.neighbors(u))) # num edges without u!
            for v in G:
                if u == v: continue
                if v in G.neighbors(u): sum_degrees_of_neighbors += log(len(G.neighbors(v))-1) if len(G.neighbors(v)) > 1 else IMPOSSIBLE
                else: sum_degrees_of_neighbors += log(num_edges - len(G.neighbors(v)))

            if sum_degrees_of_neighbors > current_sum:
                current_sum = sum_degrees_of_neighbors
                current_list = [u]
            elif sum_degrees_of_neighbors == current_sum:
                current_list.append(u)

        u = random.choice(current_list) if len(current_list) > 0 else random.choice(G.nodes())
        G.remove_node(u)

        print "%s" %(u)

    last_node = G.nodes()[0]
    print "%s" %(last_node)



#==============================================================================
#                                   MAIN
#==============================================================================
def main():

    # =============== Handle arguments and options ==================
    start = time.time()
    logging.basicConfig(
        level=logging.CRITICAL,
        format='%(levelname)s: %(asctime)s -- %(message)s'
    )

    usage="usage: %prog [options] <network file> <model parameters>"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--model", action="store", dest="model", type="string",default="dmc",help="model to reverse with: dmc, ff, or pa")

    (options, args) = parser.parse_args()
    model = options.model
    # ===============================================================


    # =================== Run NetArch algorithm =====================
    if len(args) == 0:
        logging.critical("No graph specified. Exiting...")
        sys.exit(1)

    G = Graph.read_graph(args[0])

    # DMC
    if model == "dmc" or model == "DMC":
        if len(args) != 3:
            logging.critical("Illegal number of parameters (%i, expected 3): <network> <qmod> <qcon>. Exiting..." %(len(args)))
            sys.exit(1)

        qmod = float(args[1])
        qcon = float(args[2])

        if qmod < 0 or qmod > 1:
            logging.critical("Parameter 'qmod' must be between 0 and 1. Exiting...")
            sys.exit(1)
        elif qcon < 0 or qcon > 1:
            logging.critical("Parameter 'qcon' must be between 0 and 1. Exiting...")
            sys.exit(1)

        logging.info("Running ReverseDMC...")
        print "#ANCHOR\tNODE_REMOVED"
        dmc_delorean(G,qmod,qcon)

    # FF
    elif model == "ff" or model == "FF":
        if len(args) != 2:
            logging.critical("Illegal number of parameters (%i, expected 2): <network> <p>. Exiting..." %(len(args)))
            sys.exit(1)

        p = float(args[1])

        if p <= 0 or p >= 1:
            logging.critical("Parameter 'p' must be between 0 and 1. Exiting...")
            sys.exit(1)

        logging.info("Running ReverseFF...")
        print "#ANCHOR\tNODE_REMOVED"
        ff_delorean(G,p)

    # PA
    elif model == "pa" or model == "PA":
        logging.info("Running ReversePA...")
        print "#NODE_REMOVED"
        pa_delorean(G)

    else:
        logging.critical("Invalid model: %s. Exiting..." %(model))
        sys.exit(1)


    # ========================= Finish ===========================
    logging.info("Time to run: %.2f (mins)" %((time.time()-start) / 60))


if __name__ == "__main__":
    main()
