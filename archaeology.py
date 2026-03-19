#!/usr/bin/env python

"""
    Implementation of the Network Archaeology algorithms

    @requires: python 2.5
"""

from __future__ import division
from optparse import OptionParser
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
        return random.choice(self.node_map.keys())

    def neighbors(self,u):
        """ Returns set of neighbors of u. """
        return self.node_map[u].neighbors



#==============================================================================
#                           RECONSTRUCTION FUNCTIONS
#==============================================================================

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

        if len(current_list) > 0:
            u = random.choice(current_list) 
        else:
            u = random.choice(list(G.nodes()))
        
        G.remove_node(u)

        print ("%s" %(u))

    last_node = list(G.nodes())[0]

    print ("%s" %(last_node))



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

    # PA
    if model == "pa" or model == "PA":
        logging.info("Running ReversePA...")
        print ("#NODE_REMOVED")
        pa_delorean(G)

    else:
        logging.critical("Invalid model: %s. Exiting..." %(model))
        sys.exit(1)


    # ========================= Finish ===========================
    logging.info("Time to run: %.2f (mins)" %((time.time()-start) / 60))


if __name__ == "__main__":
    main()
