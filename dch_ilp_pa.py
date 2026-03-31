# -*- coding: utf-8 -*-
import argparse
import glob
import math
import os
import re
from pathlib import Path

import networkx as nx
import igraph as ig
import leidenalg

from ilp_pa import run_ILP


ILP_RE = re.compile(r"ILP_cluster-(\d+)_order0_time_([0-9]+(?:\.[0-9]+)?)\.txt$")
MAP_RE = re.compile(r"cluster-(\d+)\.map$")
ROOT_RE = re.compile(r"cluster-(\d+)\.root$")
ROOT_ORIG_RE = re.compile(r"cluster-(\d+)\.root_orig$")


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
        seed=seed,
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


def read_graph(path):
    G = nx.Graph()
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) == 2:
                G.add_edge(int(s[0]), int(s[1]))
            elif len(s) == 1:
                G.add_node(int(s[0]))
    return G


def read_solution_pairs(path):
    pairs = []
    with open(path) as f:
        for line in f:
            a, b = line.split()
            pairs.append((int(a), int(b)))
    return pairs


def read_removal_list(path):
    if not path.exists():
        return []

    vals = []
    with open(path) as f:
        for line in f:
            vals.append(int(line.split()[0]))
    return vals


def read_int_list(path):
    if not path.exists():
        return []

    vals = []
    with open(path) as f:
        for line in f:
            vals.extend(int(tok) for tok in line.split())
    return vals


def read_single_int(path):
    vals = read_int_list(path)
    return vals[0] if vals else None


def read_cluster_map(path):
    if not path.exists():
        return {}

    m = {}
    with open(path) as f:
        for line in f:
            a, b = line.split()
            m[int(a)] = int(b)
    return m


def compute_loglkl(u, v, graph, p):
    if not graph.has_node(u) or not graph.has_node(v):
        return None

    U = set(graph.neighbors(u)) - {v}
    V = set(graph.neighbors(v)) - {u}

    common = len(U & V)
    unique = len(U ^ V)

    gamma = p if graph.has_edge(u, v) else 1 - p

    return (
        common * math.log(p)
        + unique * math.log(1 - p)
        + math.log(gamma)
    )


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


def get_all_networks(extant, anchor_rem):
    net_list = [extant.copy()]
    curr = extant.copy()

    for a, r in anchor_rem:
        if r not in curr:
            continue

        Nr = set(curr.neighbors(r))
        Na = set(curr.neighbors(a)) if a in curr else set()
        newN = Nr - Na - {a, r}

        curr.remove_node(r)

        if a not in curr:
            curr.add_node(a)

        for n in newN:
            curr.add_edge(a, n)

        net_list.append(curr.copy())

        if curr.number_of_nodes() <= 3:
            break

    return net_list


def root_is_referenced(solution_pairs, local_root):
    if local_root is None:
        return True

    for u, v in solution_pairs:
        if u == local_root or v == local_root:
            return True

    return False


def build_cluster_info(solution_dir):
    cluster_local_to_original = {}
    cluster_root_local = {}
    invalid_clusters = set()

    for mp in sorted(solution_dir.glob("cluster-*.map")):
        mm = MAP_RE.match(mp.name)
        if not mm:
            continue
        cid = int(mm.group(1))
        cluster_local_to_original[cid] = read_cluster_map(mp)

    for rp in sorted(solution_dir.glob("cluster-*.root")):
        mm = ROOT_RE.match(rp.name)
        if not mm:
            continue
        cid = int(mm.group(1))
        cluster_root_local[cid] = read_single_int(rp)

    for rop in sorted(solution_dir.glob("cluster-*.root_orig")):
        mm = ROOT_ORIG_RE.match(rop.name)
        if not mm:
            continue

        cid = int(mm.group(1))
        root_orig = read_single_int(rop)
        local_root = cluster_root_local.get(cid)
        local_to_orig = cluster_local_to_original.get(cid, {})

        if local_root is None or root_orig is None:
            invalid_clusters.add(cid)
            continue

        if local_root not in local_to_orig:
            invalid_clusters.add(cid)
            continue

        if local_to_orig[local_root] != root_orig:
            invalid_clusters.add(cid)

    return cluster_local_to_original, cluster_root_local, invalid_clusters


def translate_pairs_local_to_original(solution_pairs, cluster_id, cluster_local_to_original):
    conv = cluster_local_to_original.get(cluster_id, {})
    out = []

    for u, v in solution_pairs:
        U = conv.get(u)
        V = conv.get(v)

        if U is None or V is None:
            continue
        if U == V:
            continue

        out.append((U, V))

    return out


def choose_anchor_for_removed_node(G, removed_node):
    if removed_node not in G:
        return None

    nbrs = sorted(G.neighbors(removed_node))
    if not nbrs:
        return None

    return nbrs[0]


def greedy_order_to_pairs(extant_graph, removal_order):
    curr = extant_graph.copy()
    pairs = []

    for r in removal_order:
        if r not in curr:
            continue

        a = choose_anchor_for_removed_node(curr, r)
        if a is None:
            continue

        pairs.append((a, r))

        Nr = set(curr.neighbors(r))
        Na = set(curr.neighbors(a)) if a in curr else set()
        newN = Nr - Na - {a, r}

        curr.remove_node(r)

        if a not in curr:
            curr.add_node(a)

        for n in newN:
            curr.add_edge(a, n)

        if curr.number_of_nodes() <= 3:
            break

    return pairs


def find_greedy_file(extant_dir):
    files = sorted(extant_dir.glob("ILP_greedy_*.txt"))
    return files[0] if files else None


def find_solution_dirs(extant_dir, N):
    out = []
    for run in range(1, N + 1):
        cand = extant_dir / f"{extant_dir.name}_{run}"
        if cand.is_dir():
            out.append(cand)
    return out


def get_time_values(solution_dir):
    times = set()

    for pth in solution_dir.glob("ILP_cluster-*_order0_time_*.txt"):
        m = ILP_RE.match(pth.name)
        if m:
            times.add(float(m.group(2)))

    return sorted(times)


def score_cluster_solution_for_time(extant_graph, solution_dir, t, p):
    cluster_local_to_original, cluster_root_local, invalid_clusters = build_cluster_info(solution_dir)

    cluster_files = []
    for pth in solution_dir.glob("ILP_cluster-*_order0_time_*.txt"):
        m = ILP_RE.match(pth.name)
        if not m:
            continue
        cid = int(m.group(1))
        t_file = float(m.group(2))
        if t_file != float(t):
            continue
        cluster_files.append((cid, pth))

    cluster_files.sort(key=lambda x: x[0])

    all_original_pairs = []

    for cid, pth in cluster_files:
        if cid in invalid_clusters:
            continue

        local_pairs = read_solution_pairs(pth)
        if not local_pairs:
            continue

        local_root = cluster_root_local.get(cid)
        if not root_is_referenced(local_pairs, local_root):
            continue

        original_pairs = translate_pairs_local_to_original(
            local_pairs,
            cid,
            cluster_local_to_original,
        )

        if not original_pairs:
            continue

        all_original_pairs.extend(original_pairs)

    if not all_original_pairs:
        return None

    networks = get_all_networks(extant_graph, all_original_pairs)
    return compute_total_loglkl(networks, all_original_pairs, p)


def score_greedy_solution(extant_graph, greedy_file, p):
    removal_order = read_removal_list(greedy_file)
    if not removal_order:
        return None

    pairs = greedy_order_to_pairs(extant_graph, removal_order)
    if not pairs:
        return None

    networks = get_all_networks(extant_graph, pairs)
    return compute_total_loglkl(networks, pairs, p)


def process_extant_dir_all(extant_dir, extant_graph, p, N):
    rows = []

    greedy_file = find_greedy_file(extant_dir)
    solution_dirs = find_solution_dirs(extant_dir, N)

    greedy_score = None
    if greedy_file is not None:
        greedy_score = score_greedy_solution(extant_graph, greedy_file, p)

    if not solution_dirs:
        rows.append([
            extant_dir.name,
            "",
            "",
            "NaN" if greedy_score is None else str(round(greedy_score, 6)),
        ])
        return rows

    for solution_dir in solution_dirs:
        times = get_time_values(solution_dir)

        if not times:
            rows.append([
                extant_dir.name,
                "",
                "",
                "NaN" if greedy_score is None else str(round(greedy_score, 6)),
            ])
            continue

        for t in times:
            cluster_score = score_cluster_solution_for_time(extant_graph, solution_dir, t, p)
            rows.append([
                extant_dir.name,
                str(t),
                "NaN" if cluster_score is None else str(round(cluster_score, 6)),
                "NaN" if greedy_score is None else str(round(greedy_score, 6)),
            ])

    return rows


def process_extant_dir_best(extant_dir, extant_graph, p, N):
    rows = []

    greedy_file = find_greedy_file(extant_dir)
    solution_dirs = find_solution_dirs(extant_dir, N)

    greedy_score = None
    if greedy_file is not None:
        greedy_score = score_greedy_solution(extant_graph, greedy_file, p)

    if not solution_dirs:
        rows.append([
            extant_dir.name,
            "",
            "",
            "NaN" if greedy_score is None else str(round(greedy_score, 6)),
        ])
        return rows

    best_t = None
    best_cluster_score = None

    for solution_dir in solution_dirs:
        times = get_time_values(solution_dir)

        for t in times:
            cluster_score = score_cluster_solution_for_time(extant_graph, solution_dir, t, p)
            if cluster_score is None:
                continue

            if best_cluster_score is None or cluster_score > best_cluster_score:
                best_cluster_score = cluster_score
                best_t = t

    rows.append([
        extant_dir.name,
        "" if best_t is None else str(best_t),
        "NaN" if best_cluster_score is None else str(round(best_cluster_score, 6)),
        "NaN" if greedy_score is None else str(round(greedy_score, 6)),
    ])

    return rows


def write_rows(out_dir, rows, out_name):
    out_path = out_dir / out_name

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("extant_name\ttime_h\tcluster_total_loglkl\tgreedy_loglkl\n")

        if rows is None:
            rows = []

        for row in rows:
            if row is None:
                continue
            if not isinstance(row, (list, tuple)):
                row = [row]
            f.write("\t".join(map(str, row)) + "\n")

    return out_path


def generate_clusters(extant_dir, extant_name, extant_graph, N, gamma):
    for run in range(1, N + 1):
        communities = leiden_partitioning(extant_graph, gamma=gamma, seed=run)
        subgraphs = [extant_graph.subgraph(c).copy() for c in communities]
        roots_orig = select_roots(subgraphs, extant_graph)

        run_dir = extant_dir / f"{extant_name}_{run}"
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
                for u, v in sg.edges():
                    fe.write(f"{local[u]} {local[v]}\n")


def run_ilps(extant_dir, extant_name, N, time_limits, p):
    for run in range(1, N + 1):
        folder_str = f"{extant_dir}/{extant_name}_{run}"

        for t in time_limits:
            for fpath in sorted(glob.glob(os.path.join(folder_str, "cluster-*.edges"))):
                short_name = os.path.splitext(os.path.basename(fpath))[0]
                extantFNAME = fpath
                outputFNAME = folder_str + "/ILP_" + short_name

                run_ILP(
                    extantFNAME,
                    outputFNAME,
                    p,
                    timeLimit=t,
                    numSolutions=10,
                    poolMode=1,
                    focus=1,
                    numThreads=0,
                    memusage=10,
                )


def run_dch_ilp_pa(root_dir, extant_name, time_limits, N, gamma, p):
    extant_dir = Path(root_dir) / extant_name
    extant_file = extant_dir / extant_name

    if not extant_dir.is_dir():
        raise FileNotFoundError(f"Missing directory: {extant_dir}")
    if not extant_file.is_file():
        raise FileNotFoundError(f"Missing extant file: {extant_file}")

    extant_graph = build_graph(read_input(extant_file), directed=False)

    generate_clusters(extant_dir, extant_name, extant_graph, N, gamma)
    run_ilps(extant_dir, extant_name, N, time_limits, p)

    rows_all = process_extant_dir_all(extant_dir, extant_graph, p, N)
    out_all = write_rows(extant_dir, rows_all, "results_by_timelimit.tsv")

    rows_best = process_extant_dir_best(extant_dir, extant_graph, p, N)
    out_best = write_rows(extant_dir, rows_best, "results_best.tsv")

    print("Finished:", extant_dir)
    print("All results:", out_all)
    print("Best results:", out_best)



def main():
    
    parser = argparse.ArgumentParser(description="DCH-ILP for PA")
  
    parser.add_argument("-r", "--root", required=True, help="Parent directory containing the extant folder.")
    parser.add_argument("-e","--extant", required=True, help="Extant filename. Format: one edge per line")
    parser.add_argument("-t", "--timelimit", type=float, nargs="+", default=[4], help="Time limits in hours.")
    parser.add_argument("-n","--N", type=int, default=10, help="Number of runs / partitions to generate.")
    parser.add_argument("-g", "--gamma", type=float, required=True, help="Leiden resolution parameter.")
    parser.add_argument("-p","--p_value", type=float, default=0.5, help="p value for the model ILP-PA")
    
    args,_= parser.parse_known_args()
    
    extantFNAME = args.extant 
    root_dir = args.root
    p_value = args.p_value
    timeLimit = args.timelimit * 3600
    N = args.N
    gamma = args.gamma
    run_dch_ilp_pa(
        root_dir,
        extantFNAME,
        timeLimit,
        N,
        gamma,
        p_value,
    )


if __name__ == "__main__":
    main()