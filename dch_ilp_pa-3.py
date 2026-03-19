import re
import math
from pathlib import Path
import networkx as nx

P = 0.5

ROOT_DIRS = [
    Path("clusters_for_table_2"),
    Path("clusters_for_table_3"),
    Path("clusters_for_table_4"),
]

SOLUTION_SUFFIX = "_1"

ILP_RE = re.compile(r"ILP_cluster-(\d+)_order0_time_(\d+)\.txt$")
MAP_RE = re.compile(r"cluster-(\d+)\.map$")
ROOT_RE = re.compile(r"cluster-(\d+)\.root$")
ROOT_ORIG_RE = re.compile(r"cluster-(\d+)\.root_orig$")


def read_graph(path):
    G = nx.Graph()
    with open(path) as f:
        for line in f:
            s = line.strip().split()
            if not s:
                continue
            if s[0].startswith("#"):
                continue

            vals = []
            ok = True
            for x in s:
                if x.lstrip("-").isdigit():
                    vals.append(int(x))
                else:
                    ok = False
                    break

            if not ok or not vals:
                continue

            if len(vals) >= 2:
                G.add_edge(vals[0], vals[1])
            elif len(vals) == 1:
                G.add_node(vals[0])

    return G


def read_solution_pairs(path):
    pairs = []
    with open(path) as f:
        for line in f:
            s = line.strip().split()
            if not s:
                continue
            if s[0].startswith("#"):
                continue

            vals = []
            for x in s:
                if x.lstrip("-").isdigit():
                    vals.append(int(x))

            if len(vals) >= 2:
                pairs.append((vals[0], vals[1]))

    return pairs


def read_removal_list(path):
    vals = []
    if not path.exists():
        return vals

    with open(path) as f:
        for line in f:
            s = line.strip().split()
            if not s:
                continue
            if s[0].startswith("#"):
                continue

            first_num = None
            for tok in s:
                if tok.lstrip("-").isdigit():
                    first_num = int(tok)
                    break

            if first_num is not None:
                vals.append(first_num)

    return vals


def read_int_list(path):
    vals = []
    if not path.exists():
        return vals

    with open(path) as f:
        for line in f:
            s = line.strip().split()
            for tok in s:
                if tok.lstrip("-").isdigit():
                    vals.append(int(tok))

    return vals


def read_single_int(path):
    vals = read_int_list(path)
    if vals:
        return vals[0]
    return None


def read_cluster_map(path):
    m = {}
    if not path.exists():
        return m

    with open(path) as f:
        for line in f:
            s = line.strip().split()
            if len(s) < 2:
                continue
            if not s[0].lstrip("-").isdigit():
                continue
            if not s[1].lstrip("-").isdigit():
                continue
            m[int(s[0])] = int(s[1])

    return m


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

    loglkl = (
        common * math.log(p) +
        unique * math.log(1 - p) +
        math.log(gamma)
    )
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


def find_extant_file(extant_dir):
    for cand in extant_dir.iterdir():
        if cand.is_file() and not cand.name.startswith("ILP_greedy_"):
            return cand
    return None


def find_greedy_file(extant_dir):
    files = sorted(extant_dir.glob("ILP_greedy_*.txt"))
    if files:
        return files[0]
    return None


def find_solution_dir(extant_dir):
    for cand in extant_dir.iterdir():
        if cand.is_dir() and cand.name.endswith(SOLUTION_SUFFIX):
            return cand
    return None


def get_time_values(solution_dir):
    times = set()

    for pth in solution_dir.glob("ILP_cluster-*_order0_time_*.txt"):
        m = ILP_RE.match(pth.name)
        if m:
            times.add(int(m.group(2)))

    return sorted(times)


def score_cluster_solution_for_time(extant_graph, solution_dir, t):
    cluster_local_to_original, cluster_root_local, invalid_clusters = build_cluster_info(solution_dir)

    total = 0.0
    seen = 0

    cluster_files = []
    for pth in solution_dir.glob(f"ILP_cluster-*_order0_time_{t}.txt"):
        m = ILP_RE.match(pth.name)
        if not m:
            continue
        cid = int(m.group(1))
        cluster_files.append((cid, pth))

    cluster_files.sort(key=lambda x: x[0])

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
            cluster_local_to_original
        )

        if not original_pairs:
            continue

        networks = get_all_networks(extant_graph, original_pairs)
        sc = compute_total_loglkl(networks, original_pairs, P)

        total += sc
        seen += 1

    if seen == 0:
        return None

    return total


def score_greedy_solution(extant_graph, greedy_file):
    removal_order = read_removal_list(greedy_file)
    if not removal_order:
        return None

    pairs = greedy_order_to_pairs(extant_graph, removal_order)
    if not pairs:
        return None

    networks = get_all_networks(extant_graph, pairs)
    return compute_total_loglkl(networks, pairs, P)


def process_root_dir(root_dir):
    rows = []

    for extant_dir in sorted(root_dir.iterdir()):
        if not extant_dir.is_dir():
            continue

        extant_file = find_extant_file(extant_dir)
        if extant_file is None:
            continue

        greedy_file = find_greedy_file(extant_dir)
        solution_dir = find_solution_dir(extant_dir)

        extant_graph = read_graph(extant_file)

        greedy_score = None
        if greedy_file is not None:
            greedy_score = score_greedy_solution(extant_graph, greedy_file)

        if solution_dir is None:
            rows.append([
                extant_dir.name,
                "",
                "",
                "NaN" if greedy_score is None else str(round(greedy_score, 6))
            ])
            continue

        times = get_time_values(solution_dir)

        if not times:
            rows.append([
                extant_dir.name,
                "",
                "",
                "NaN" if greedy_score is None else str(round(greedy_score, 6))
            ])
            continue

        for t in times:
            cluster_score = score_cluster_solution_for_time(extant_graph, solution_dir, t)

            rows.append([
                extant_dir.name,
                str(t),
                "NaN" if cluster_score is None else str(round(cluster_score, 6)),
                "NaN" if greedy_score is None else str(round(greedy_score, 6))
            ])

    return rows


def write_rows(root_dir, rows):
    out_path = root_dir / "results_by_timelimit.tsv"

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("extant_name\ttime_h\tcluster_total_loglkl\tgreedy_loglkl\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    return out_path


if __name__ == "__main__":
    for root_dir in ROOT_DIRS:
        if not root_dir.exists():
            print("Skipping missing folder:", root_dir)
            continue

        rows = process_root_dir(root_dir)
        out_path = write_rows(root_dir, rows)

        print("Finished:", root_dir)
        print("Output:", out_path)
