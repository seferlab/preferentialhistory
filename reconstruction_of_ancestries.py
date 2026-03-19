# -*- coding: utf-8 -*-
import os
import math
import networkx as nx
import numpy as np
import pandas as pd
from grakel import GraphKernel, Graph
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
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

def analyze_all(data, n):
    results = []

    folder = f'data/{data}'
    base = f"{folder}/{data}"
    
    extant = read_graph(base)
    if len(extant.nodes()) == 0:
        return pd.DataFrame()

    true_sequence = extant.edges
    prob_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    
    for p in prob_values:
        best_loglkl = -float('inf')
        best_ar = None
        best_nets = None
        best_order = -1
        
        for order in range(1):
            file_path = f"{folder}/ILP_{data}_p={p}_order{order}.txt"
            
            if not os.path.isfile(file_path):
                continue
                
            ar = read_ARlist(file_path)
            nets = get_all_networks(extant, ar)
            
            loglkl = compute_total_loglkl(nets, ar, p)
            
            if loglkl is not None and loglkl > best_loglkl:
                best_loglkl = loglkl
                best_ar = ar
                best_nets = nets
                best_order = order
                
        if best_ar is not None:
            tau = compute_kendalls_tau(true_sequence, best_ar)
            kernel = compute_kernel_similarity(best_nets, extant, n)
            
            print(f"[{data}] For p={p}, best file was order{best_order} (LogLKL: {best_loglkl:.2f})")
            
            results.append({
                "n": n, 
                "p": p, 
                "Metric": f"p = {p}",
                "LogLKL": best_loglkl, 
                "KendallTau": tau, 
                "KernelSim": kernel
            })
        else:
            print(f"Warning: No valid files found for {data} at p={p}")
            
    df = pd.DataFrame(results)
    return df

def generate_final_summary(dfs_dict):
    """Generate combined summary for multiple datasets"""
    combined_df = pd.concat([df.assign(data=data) for data, df in dfs_dict.items()], ignore_index=True)
    
    pivot_df = combined_df.pivot_table(index='p', columns='data', 
                                        values=['LogLKL', 'KendallTau', 'KernelSim'],
                                        aggfunc='first')
    
    pivot_df.columns = [f"{col[1]}_{col[0]}" for col in pivot_df.columns]
    
    print(pivot_df)
    
    pivot_df.to_csv("data/final_summary_combined.csv")

    plt.figure(figsize=(14, 6))
    plt.axis('off')
    
    header_data = ["p"]
    for data in dfs_dict.keys():
        header_data.extend([f"{data}\nLog-likelihood", f"{data}\nKendall's Tau", f"{data}\nKernel Similarity"])
    
    cell_text = []
    for idx, row in pivot_df.iterrows():
        row_data = [f"p = {idx}"]
        for data in sorted(dfs_dict.keys()):
            row_data.extend([
                f"{row[f'{data}_LogLKL']:.2f}" if pd.notna(row[f'{data}_LogLKL']) else "N/A",
                f"{row[f'{data}_KendallTau']:.2f}" if pd.notna(row[f'{data}_KendallTau']) else "N/A",
                f"{row[f'{data}_KernelSim']:.2f}" if pd.notna(row[f'{data}_KernelSim']) else "N/A"
            ])
        cell_text.append(row_data)
    
    table = plt.table(
        cellText=cell_text,
        colLabels=header_data,
        cellLoc='center',
        loc='center',
        colColours=['#f0f0f0'] * len(header_data)
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    plt.tight_layout()
    plt.savefig("data/table_combined.png", bbox_inches='tight', dpi=600)
    plt.close()

if __name__ == "__main__":
    
    data_list = [['bzip', 13], ['commander', 9]]
    dfs_dict = {}
    
    for data_inner_list in data_list:
        df = analyze_all(data_inner_list[0], data_inner_list[1])
        dfs_dict[data_inner_list[0]] = df
    
    generate_final_summary(dfs_dict)
