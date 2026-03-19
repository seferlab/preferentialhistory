import pandas as pd
import re

def read_mapping(path):
    df = pd.read_csv(path, sep=r"\s+|\t", engine="python", header=None, names=["Gene","no"])
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df["no"]   = pd.to_numeric(df["no"], errors="coerce").astype("Int64")
    return df.dropna(subset=["no"]).astype({"no": int})

def read_orthologs(path):
    df = pd.read_csv(path, sep=r"\s+|\t", engine="python", header=None, names=["Gene","Orthologs"])
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df["Orthologs"] = pd.to_numeric(df["Orthologs"], errors="coerce").fillna(0).astype(int)
    return df

def read_reverse_nodes(path, mapping):
    df = pd.read_csv(path, sep=r"\s+|\t", engine="python", header=None, names=["no"])
    df["no"] = df["no"].astype(str).str.extract(r"(\d+)").astype(float).astype("Int64")
    df = df.dropna(subset=["no"]).astype({"no": int})
    df = df.merge(mapping, on="no", how="left")
    df["step"] = range(1, len(df)+1)
    return df[["Gene","step"]]

def edge_based_time_steps(path, mapping):
    step_counter = 1
    node_steps = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            a = re.findall(r"\d+", line)
            if len(a) < 2:
                continue
            u, v = int(a[0]), int(a[1])

            if u not in node_steps:
                node_steps[u] = step_counter
                step_counter += 1

            if v not in node_steps:
                node_steps[v] = step_counter
                step_counter += 1

    df = pd.DataFrame(list(node_steps.items()), columns=["no", "step"])
    df = df.merge(mapping, on="no", how="left")[["Gene","step"]]
    return df

def make_table_txt(reverse_file, ilp1_file, ilp2_file, mapping_file, ortholog_file, output_txt):
    mapping = read_mapping(mapping_file)
    orth = read_orthologs(ortholog_file).sort_values("Orthologs", ascending=False).reset_index(drop=True)
    orth["TimeStep"] = range(1, len(orth)+1)

    rev_df  = read_reverse_nodes(reverse_file, mapping)
    ilp1_df = edge_based_time_steps(ilp1_file, mapping)
    ilp2_df = edge_based_time_steps(ilp2_file, mapping)

    def step_to_text(df):
        g = df.groupby("step")["Gene"].apply(list).to_dict()
        return {k: ", ".join(sorted(v, key=str)) for k,v in g.items()}

    rev_step = step_to_text(rev_df)
    i1_step  = step_to_text(ilp1_df)
    i2_step  = step_to_text(ilp2_df)

    def lookup(step, d): return d.get(step, "")

    table = orth.copy()
    table["ReverseLPA"] = table["TimeStep"].map(lambda s: lookup(s, rev_step))
    table["ILP1"]       = table["TimeStep"].map(lambda s: lookup(s, i1_step))
    table["ILP2"]       = table["TimeStep"].map(lambda s: lookup(s, i2_step))

    table = table[["Gene","Orthologs","TimeStep","ReverseLPA","ILP1","ILP2"]]

    with open(output_txt, "w", encoding="utf-8") as f:
        header = f"{'Gene':<10} {'Orthologs':<10} {'TimeStep':<10} {'ReverseLPA':<25} {'ILP1':<25} {'ILP2':<25}\n"
        f.write(header)
        f.write("=" * len(header) + "\n")
        for _, row in table.iterrows():
            if pd.isna(row["Gene"]) or str(row["Gene"]).strip().lower() in ["nan", "none", "gene", "0"]:
                continue
            f.write(f"{row['Gene']:<10} {row['Orthologs']:<10} {row['TimeStep']:<10} {row['ReverseLPA']:<25} {row['ILP1']:<25} {row['ILP2']:<25}\n")
    print(table.head(5))

if __name__ == "__main__":
    for file in ['bzip', 'commander']:
        folder = f"data/{file}/"
        if file == 'bzip':
            make_table_txt(
                reverse_file   = folder + f"ILP_greedy_{file}.txt",
                ilp1_file      = folder + f"ILP_{file}_order0.txt",
                ilp2_file      = folder + f"ILP_{file}_order1.txt",
                mapping_file   = folder + f"{file}_mapping.txt",
                ortholog_file  = folder + f"orthologs_{file}.txt",
                output_txt     = folder + f"comparison_table_{file}.txt"
            )
        else:
            make_table_txt(
                reverse_file   = folder + f"ILP_greedy_{file}.txt",
                ilp1_file      = folder + f"ILP_{file}_order0.txt",
                ilp2_file      = folder + f"ILP_{file}_order3.txt",
                mapping_file   = folder + f"{file}_mapping.txt",
                ortholog_file  = folder + f"orthologs_{file}.txt",
                output_txt     = folder + f"comparison_table_{file}.txt"
            )
            
