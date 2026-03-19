# -*- coding: utf-8 -*-
import networkx as nx
import random
import pathlib

def generate_graph(num_nodes, edge_prob):
    """Each new node connects to existing nodes with a given probability."""
    G = nx.Graph()
    G.add_node(1)
    for newnode in range(2, num_nodes + 1):
        G.add_node(newnode)
        for existing in G.nodes():
            if existing == newnode:
                continue
            if random.random() < edge_prob:
                G.add_edge(newnode, existing)
        # Ensure the new node is not isolated
        if G.degree[newnode] == 0:
            other = random.choice([n for n in G.nodes() if n != newnode])
            G.add_edge(newnode, other)
    return G

def write_graph(filename, graph):
    path = pathlib.Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        for u, v in graph.edges():
            f.write(f"{u} {v}\n")

# ==== CONFIGURATION ====
probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
nrange = range(6, 11)  # n = 6, 7, 8, 9, 10
runrange = range(1, 34)
# ==== MAIN LOOP ====
  
for run in runrange:
    folder = 'run' + str(run)
    pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    for p in probL:
        for n in nrange:
            G = generate_graph(n, p)
            filename = f'{folder}/nx_run={run}_p={p}_n={n}'
            print("Writing:", filename)
            write_graph(filename, G)


for_t2 = [50, 100, 250, 500,1000]
folder = "clusters_for_table_2"
for n in for_t2:
    G = generate_graph(n, 0.5)
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    pathlib.Path(folder_new).mkdir(parents=True, exist_ok=True)
    filename = f'{folder_new}/nx_p=0.5_n={n}'
    print("Writing:", filename)
    write_graph(filename, G)


for_t3 = [12, 20, 50, 100]
folder = "clusters_for_table_3"
for n in for_t3:
    G = generate_graph(n, 0.5)
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    pathlib.Path(folder_new).mkdir(parents=True, exist_ok=True)
    filename = f'{folder_new}/nx_p=0.5_n={n}'
    print("Writing:", filename)
    write_graph(filename, G)



for_t4 = [200]
folder = "clusters_for_table_4"
for n in for_t4:
    G = generate_graph(n, 0.5)
    folder_new = f'{folder}/nx_p=0.5_n={n}'
    pathlib.Path(folder_new).mkdir(parents=True, exist_ok=True)
    filename = f'{folder_new}/nx_p=0.5_n={n}'
    print("Writing:", filename)
    write_graph(filename, G)















