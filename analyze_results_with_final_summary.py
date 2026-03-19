import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy

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

def generate_final_summary(df):
    metrics = ["LogLKL", "KendallTau", "KernelSim"]
    summary = df.groupby("Metric")[metrics].mean().round(2)
    print("--- Summary Data ---")
    print(summary)
    
    # 1. Plot as Table
    plt.figure(figsize=(8, 2.5))
    plt.axis('off')
    cell_text = [[f"{row[m]}" for m in metrics] for _, row in summary.iterrows()]

    table = plt.table(
        cellText=cell_text,
        rowLabels=summary.index.tolist(),
        colLabels=["Log-likelihood", "Kendall's Tau", "Kernel Similarity"],
        cellLoc='center', loc='center', colColours=['#f0f0f0'] * 3
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    plt.tight_layout()
    plt.savefig("table-1.png", bbox_inches='tight', dpi=300)
    plt.close()

    # 2. Metric Comparison (ILP vs ReverseLPA)
    df_ilp = df[df["Metric"] == "ILP-PA1"].reset_index(drop=True)
    df_rev = df[df["Metric"] == "ReverseLPA"].reset_index(drop=True)
    
    min_len = min(len(df_ilp), len(df_rev))
    df_merged = pd.concat([
           df_ilp.iloc[:min_len][metrics].add_suffix("_ilp").reset_index(drop=True),
           df_rev.iloc[:min_len][metrics].add_suffix("_rev").reset_index(drop=True)
    ], axis=1)

    categories = ['ILP < G', 'ILP = G', 'ILP > G']
    colors = ['yellow', 'green', 'blue']
    counts = {metric: [0, 0, 0] for metric in metrics}
    tolerance = 1e-4
    
    for metric in metrics:
        for i in range(len(df_merged)):
            diff = df_merged.loc[i, f"{metric}_ilp"] - df_merged.loc[i, f"{metric}_rev"]
            if abs(diff) < tolerance: counts[metric][1] += 1
            elif diff < 0: counts[metric][0] += 1
            else: counts[metric][2] += 1
            

    percent_counts = {m: [round(100 * v / sum(counts[m]), 2) for v in counts[m]] for m in metrics}
    
    
    fig, ax = plt.subplots(figsize=(8, 4))
    y_pos = np.arange(len(metrics))
    left = np.zeros(len(metrics))
    
    for i, (label, color) in enumerate(zip(categories, colors)):
        vals = [percent_counts[metrics[len(metrics) - m - 1]][i] for m in range(len(metrics))]
        if any(vals):
            ax.barh(y_pos, vals, left=left, color=color, height=0.6, label=label)
            left += vals
            
    ax.set_yticks(y_pos)
    ax.set_yticklabels(["Sim", "KTau", "LKL"])
    ax.set_xlim(0, 100)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
    plt.tight_layout()
    plt.savefig("synthetic_result1.png", dpi=600)
    plt.close()


def plot_scatter_graphs(df):
    metrics = ["LogLKL", "KendallTau", "KernelSim"]
    titles = ["Log-likelihood", "Kendall's Tau", "Kernel Similarity"]
    colors = ["blue", "green", "red"]

    df_ilp = df[df["Metric"] == "ILP-PA1"]
    df_rev = df[df["Metric"] == "ReverseLPA"]
    
    on_cols = [c for c in ['n', 'p', 'run'] if c in df.columns]
    
    df_merged = pd.merge(df_ilp, df_rev, on=on_cols, suffixes=('_ilp', '_rev'))

    for i, metric in enumerate(metrics):
        x = df_merged[f"{metric}_rev"]
        y = df_merged[f"{metric}_ilp"]
        
        plt.figure(figsize=(5, 5))
        plt.scatter(x, y, color=colors[i], marker='x', s=1, alpha=0.6)
        
        low = min(x.min(), y.min())
        high = max(x.max(), y.max())
        plt.plot([low, high], [low, high], 'k-', linewidth=1)
        
        plt.xlabel("ReverseLPA")
        plt.ylabel("ILPPA")
        plt.title(titles[i])
        plt.tight_layout()
        plt.savefig(f"synthetic2_{metric}.png", dpi=600)
        plt.close()
  
    
  
def plot_reconstructions():
    folder = "figure-3"
    base_graph_path = f"{folder}/nx_p=0.5_n=6"
    G_base = read_graph(base_graph_path)
    
    configs = {
        "ReverseLPA": f"{folder}/ILP_greedy_nx_p=0.5_n=6.txt",
        "ILP-PA": f"{folder}/ILP_nx_p=0.5_n=6_order5.txt",
        "Truth": list(G_base.edges())
    }
    
    for method, sequence in configs.items():
        if method == "ReverseLPA":
            with open(configs[method]) as f:
                node_sequence = [int(line.strip()) for line in f if line.strip()]
            
            G_recon = nx.Graph()
            steps = []
            
            for node in node_sequence:
                if node not in G_recon:
                    G_recon.add_node(node)
            
            G = nx.Graph()
            for node in node_sequence:
                for neighbor in G_base.neighbors(node):
                    if neighbor in G_recon.nodes and not G.has_edge(node, neighbor):
                        G.add_edge(node, neighbor)
                        steps.append(deepcopy(G))
        
        else:  
            G_recon = nx.Graph()
            steps = []
            
            if method == "ILP-PA":
                with open(sequence) as f:
                    edge_sequence = []
                    for line in f:
                        parts = list(map(int, line.strip().split()))
                        if len(parts) >= 2:
                            edge_sequence.append((parts[0], parts[1]))
                    sequence = edge_sequence
    
            for u, v in sequence:
                if u not in G_recon:
                    G_recon.add_node(u)
                if v not in G_recon:
                    G_recon.add_node(v)
                if not G_recon.has_edge(u, v):
                    G_recon.add_edge(u, v)
                    steps.append(deepcopy(G_recon))
            
        pos = nx.spring_layout(G_base, seed=42)
        if len(steps) > 1:  
            fig, axes = plt.subplots(nrows=len(steps), ncols=1, 
                                   figsize=(20, 13*(len(steps)-1)))
            
            if len(steps)-1 == 1:
                axes = [axes]
           
            
            for idx, ax in enumerate(axes, start=1):  
                # if idx == 1:
                #     ax.set_title(f"{method}", fontsize=10)
                G_step = steps[idx-1]
                nx.draw_networkx_nodes(G_step, pos, node_color='lightgray', node_size=20000, ax=ax)                
                edge_colors = ['black'] * len(G_step.edges())
                nx.draw_networkx_edges(G_step, pos, edge_color=edge_colors, width=0.7, ax=ax)
                nx.draw_networkx_labels(G_step, pos, font_size=60, font_weight="normal", ax=ax)
                ax.axis('off')
            

            plt.tight_layout()
            plt.savefig(f"{folder}/{method.lower()}.png", dpi=300)
            plt.close()
    
    
if __name__ == "__main__":
    df = pd.read_csv("final_summary.csv")
    generate_final_summary(df)
    plot_scatter_graphs(df)
    plot_reconstructions()
