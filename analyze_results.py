import os
import math
import networkx as nx
import numpy as np
import pandas as pd
from grakel import GraphKernel, Graph
from scipy.stats import kendalltau

def read_graph(path):
    G = nx.Graph()
    with open(path) as f:
        for line in f:
            parts = list(map(int, line.strip().split()))
            if len(parts) >= 2:
                G.add_edge(parts[0], parts[1])
            elif len(parts) == 1:
                G.add_node(parts[0])
    return G

def read_ARlist(path):
    anchor_rem = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts or not parts[0].isdigit():
                continue
            parts = list(map(int, parts))
            u = parts[0]
            v = parts[1] if len(parts) > 1 else u
            anchor_rem.append((u, v))
    # print(anchor_rem)
    return anchor_rem


def get_all_networks(extant, anchor_rem):
    net_list = [extant]
    curr = extant.copy()
    for a, r in anchor_rem:
        if r not in curr:
            continue
        Nr = set(curr.neighbors(r))
        Na = set(curr.neighbors(a)) if a in curr else set()
        newN = Nr - Na - {a, r}
        curr.remove_node(r)
        for n in newN:
            curr.add_edge(a, n)
        net_list.append(curr.copy())
        if curr.number_of_nodes() <= 3:
            break
    # print(net_list[0].edges)
    # print(net_list[1].edges)
    # print(len(net_list))
    return net_list

def make_it_seq(seq):
    result = []
    for u, v in seq:
        for node in (u, v):
            if node not in result:
                result.append(node)
    return result

def compute_loglkl(u, v, graph, p):
    if not graph.has_node(u) or not graph.has_node(v):
        return None

    U = set(graph.neighbors(u)) - {v}
    V = set(graph.neighbors(v)) - {u}

    common = len(U & V)
    unique = len(U ^ V)

    if graph.has_edge(u, v):
        gamma = p
    else:
        gamma = 1 - p
    
    loglkl = (common * math.log(p) + unique * math.log(1 - p) + math.log(gamma))
    return loglkl

def compute_total_loglkl(networks, anchor_rem, p):
    total_loglkl = 0.0

    for i, (u, v) in enumerate(anchor_rem):
        if i >= len(networks):
            break
        G = networks[i]
        result = compute_loglkl(u, v, G, p)
        if result is not None:
            total_loglkl += result

    return total_loglkl


def compute_kendalls_tau(true_sequence, pred_sequence):
    true_removed = make_it_seq(true_sequence)
    pred_removed = make_it_seq(pred_sequence)
    # print(true_removed)
    # print(pred_removed)

    common_nodes = list(set(true_removed) & set(pred_removed))

    filtered_true = [v for v in true_removed if v in common_nodes]
    filtered_pred = [v for v in pred_removed if v in common_nodes]

    if len(filtered_true) < 2 or len(filtered_pred) < 2 or len(filtered_true) != len(filtered_pred):
        return 0.0

    tau, _ = kendalltau(filtered_true, filtered_pred)
    # print(tau)
    return tau

def relabel_graph_sequentially(G):
    nodes = sorted(G.nodes())
    mapping = {u:i for i,u in enumerate(nodes)}
    H = nx.relabel_nodes(G, mapping, copy=True)
    rev = {i:u for u,i in mapping.items()}
    return H, rev

def to_grakel_graph(G):
    n = G.number_of_nodes()
    if n == 0:
        return None
    H, _ = relabel_graph_sequentially(G)
    A = nx.to_numpy_array(H, nodelist=sorted(H.nodes()), dtype=np.int8)
    labels = {i: str(i) for i in range(n)}
    return Graph(A, node_labels=labels)

def compute_kernel_similarity(networks, extant, n):
    wl = GraphKernel(
        kernel=[{"name": "weisfeiler_lehman", "n_iter": 5}, {"name": "subtree_wl"}],
        normalize=True
    )
    refG, _ = relabel_graph_sequentially(extant)
    scores = []
    for i in range(1, len(networks)):
        k = n - i
        if k <= 0:
            break
        keep = sorted(refG.nodes())[:k]
        true_sub = refG.subgraph(keep).copy()

        G1 = to_grakel_graph(true_sub)
        G2 = to_grakel_graph(networks[i])
        if G1 is None or G2 is None:
            continue

        K = wl.fit_transform([G1, G2])
        scores.append(float(K[0,1]))

    if not scores:
        return np.nan
    return float(np.mean(scores))



def analyze_all(run_num= list(range(1,34)), nodes=[6,7,8,9,10], probL=None):
    if probL is None:
        probL = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    results = []
    for n in nodes:
        for run in run_num:
            for p in probL:
                folder = f"run{run}"
                base = f"{folder}/nx_run={run}_p={p:.1f}_n={n}"
                extant = read_graph(base)
                if len(extant.nodes()) == 0:
                    # print(base)
                    continue
                
                true_sequence = extant.edges

                configs = {
                    "ReverseLPA": f"{folder}/ILP_greedy_nx_run={run}_p={p:.1f}_n={n}.txt",
                    "ILP-PA1": f"{folder}/ILP_nx_run={run}_p={p:.1f}_n={n}_order0.txt",
                    "ILP-PA2": f"{folder}/ILP_nx_run={run}_p={p:.1f}_n={n}_order1.txt",
                }

                for method, file in configs.items():
                    if not os.path.isfile(file):
                        continue
                    ar = read_ARlist(file)
                    nets = get_all_networks(extant, ar)
                    loglkl = compute_total_loglkl(nets, ar, p)
                    tau = compute_kendalls_tau(true_sequence, ar)
                    # print(file)
                    kernel = compute_kernel_similarity(nets, extant, n)
                    results.append({
                        "n": n, "p": p, "Metric": method,
                        "LogLKL": loglkl, "KendallTau": tau, "KernelSim": kernel
                    })
    df = pd.DataFrame(results)
    # print(results)
    return df
    

def generate_final_summary(df):
    summary = df.groupby("Metric")[["LogLKL", "KendallTau", "KernelSim"]].mean().round(2)
    print(summary)
    df.to_csv("final_summary.csv")
if __name__ == "__main__":
    
    df = analyze_all()
    generate_final_summary(df)
