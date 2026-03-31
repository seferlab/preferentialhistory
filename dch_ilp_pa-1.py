import networkx as nx
import igraph as ig
import leidenalg
from pathlib import Path

def read_input(path):
    edges = []
    with open(path) as f:
        for line in f:
            a, b = line.split()
            edges.append((int(a), int(b)))
    return edges

def build_graph(edges, directed=False):
    G = nx.DiGraph() if directed else nx.Graph()
    G.add_edges_from(edges)
    return G

def convert_to_igraph(nx_graph):
    node2idx = {n: i for i, n in enumerate(nx_graph.nodes())}
    idx2node = {i: n for n, i in node2idx.items()}
    g = ig.Graph(directed=nx_graph.is_directed())
    g.add_vertices(len(node2idx))
    g.add_edges([(node2idx[u], node2idx[v]) for u, v in nx_graph.edges()])
    return g, idx2node

def leiden_partitioning(G_t, gamma, seed=None):
    g, idx2node = convert_to_igraph(G_t)
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=gamma,
        seed=seed
    )
    return [[idx2node[i] for i in comm] for comm in part]

def select_roots(subgraphs, G_full):
    roots = {}
    node2cid = {}
    for cid, sg in enumerate(subgraphs):
        for n in sg.nodes():
            node2cid[n] = cid
    for cid, sg in enumerate(subgraphs):
        if sg.number_of_nodes() == 0:
            roots[cid] = -1
            continue
        if sg.number_of_nodes() == 1:
            roots[cid] = next(iter(sg.nodes()))
            continue
        boundary = []
        for u in sg.nodes():
            for v in G_full.neighbors(u):
                if node2cid.get(v, cid) != cid:
                    boundary.append(u)
                    break
        if boundary:
            root = max(boundary, key=lambda x: G_full.degree(x))
        else:
            root = max(sg.nodes(), key=lambda x: sg.degree[x])
        roots[cid] = root
    return roots

if __name__ == "__main__":
    # files = [f"nx_p=0.5_n={n}" for n in [50, 100, 250, 500, 1000]]
    # files = [f"nx_p=0.5_n={n}" for n in [12, 20, 50, 100]]
    files = [f"nx_p=0.5_n={n}" for n in [200]]
    # folder = "clusters_for_table_2"
    # folder = "clusters_for_table_3"
    folder = "clusters_for_table_4"
    # outroot = Path("clusters_for_table_2")     
    # outroot = Path("clusters_for_table_3")   
    outroot = Path("clusters_for_table_4")
    outroot.mkdir(parents=True, exist_ok=True)

    N = 1        
    gamma = 1.15
    p_str = "0.5" 

    for fname in files:
        file_path = f"{folder}/{fname}/{fname}"
        
        G = build_graph(read_input(file_path), directed=False)
        base_dir = outroot / Path(fname).name
        base_dir.mkdir(parents=True, exist_ok=True)

        for run in range(1, N + 1):
            communities = leiden_partitioning(G, gamma=gamma, seed=run)
            subgraphs = [G.subgraph(c).copy() for c in communities]
            roots_orig = select_roots(subgraphs, G)

            run_dir = base_dir / f"{fname}_{run}"
            run_dir.mkdir(parents=True, exist_ok=True)

            for cid, sg in enumerate(subgraphs, start=1):
                nodes = list(sg.nodes())
                local = {u: i + 1 for i, u in enumerate(nodes)}

                root_orig = roots_orig[cid - 1]
                root_local = local.get(root_orig, -1)

               
                (run_dir / f"cluster-{cid}.map").write_text(
                    "".join(f"{local[u]} {u}\n" for u in nodes)
                )
                (run_dir / f"cluster-{cid}.root").write_text(f"{root_local}\n")
                (run_dir / f"cluster-{cid}.root_orig").write_text(f"{root_orig}\n")

                with (run_dir / f"cluster-{cid}.edges").open("w", newline="\n") as fe:
                    # fe.write(f"{Path(fname).name} run {run} cluster {cid}\n")                # 1
                    # fe.write(f"p {p_str}\n")                                                 # 2
                    # fe.write(f"nodes {sg.number_of_nodes()} edges {sg.number_of_edges()}\n") # 3
                    # fe.write("edge\n")                                                       # 4
                    for u, v in sg.edges():
                        fe.write(f"{local[u]} {local[v]}\n")