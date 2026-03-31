# ILP-PA: Network Archaeology with ILP and Greedy Algorithms

Reconstructs how networks evolved by finding the sequence of node additions under the Preferential Attachment model.

---

## Core Code Files

### `ilp_pa.py`
Integer Linear Programming solver using Gurobi. Reconstructs maximum likelihood network evolution within the Preferential Attachment model.

**Usage**:
```bash
python ilp_pa.py -e EXTANT -o OUTPUT [-p P_VALUE] [-t TIMELIMIT] [-s NUMSOLUTIONS] 
                 [-m POOLMODE] [-f FOCUS] [-n NUMCORES] [-u MEMUSAGE]
```

**Arguments**:
- `-e EXTANT`: Network file (edge list format)
- `-o OUTPUT`: Output file prefix
- `-p P_VALUE`: PA parameter p (default: 0.5)
- `-t TIMELIMIT`: Time limit in hours (default: 2)
- `-n NUMCORES`: Cores to use, 0=all (default: 0)
- `-s NUMSOLUTIONS`: Number of solutions (default: 30)
- `-m POOLMODE`: Pool mode 0/1/2 (default: 1)
- `-f FOCUS`: MIP focus 1/2/3 (default: 1)
- `-u MEMUSAGE`: Memory limit in GB (default: 10)

### `ilp_pa_1.py`
Simplified ILP solver without time limits. Runs until optimal solution found.

### `archaeology.py`
Greedy heuristic solver for network reconstruction. Fast approximation algorithm.

### `sim.py`
Generates synthetic random networks using preferential attachment model.

---

# DCH-ILP-PA

Runs the DCH-ILP-PA algortihm for a single known extant graph.

## Required arguments
- `-r, --root`  
  Parent directory containing the extant folder.

- `-e, --extant`  
  Extant filename. The extant file is expected at:

  ```text
  <root>/<extant>/<extant>

---

## Analysis Scripts

### `run_ilp_sim.py`
Runs ILP-PA on multiple networks. Time limit: 2 hours per problem, 10 solutions per network.

### `run_greedy.py`
Runs greedy algorithm on networks for baseline comparison.

### `analyze_results.py`
Analyzes solutions. Computes: Kernel Similarity, Kendall's Tau, Log-likelihood.

### `analyze_results_with_final_summary.py`
Extended analysis with aggregated statistics and rankings.

### `run_data.py`
Runs ILP solver on multiple bZip and commander networks. No time limit, 10 solutions per network.

### `table_5_7.py` & `table_6_8.py`
Generate comparison tables for the networks in the research paper.

### `dch_ilp_pa.py`
DCH-ILP-PA algortihm on extant files.
Step 1: Community detection using Leiden partitioning.
Step 2: ILP solving on subgraph of networks.
Step 3: Log-likelihood computation and analysis.

### `run_dch_ilp_pa.py`
Runs DCH-ILP-PA on mutilple networks.

### `reconstruction_of_ancestries.py`
Generate reconstruction of ancestries table for the networks in the research paper.

---

## Requirements

```bash
pip install networkx numpy pandas scipy matplotlib grakel
pip install gurobipy  # Requires Gurobi license (free academic available)
```

## Quick Usage

**Run ILP with time limit**:
```bash
python ilp_pa.py -e input_graph -o output -p 0.5
```

**Run single ILP problem**:
```bash
python ilp_pa_1.py -e input_graph -o output 0.5
```

**Generate networks**:
```bash
python sim.py
```

**Run ILP on batch**:
```bash
python run_ilp_sim.py
```

**Run greedy**:
```bash
python run_greedy.py
```

**Run ILP for the networks**:
```bash
python run_data.py
```

**Analyze results**:
```bash
python analyze_results.py
```

**Generate paper tables**:
```bash
python table_5_7.py
python table_6_8.py
python reconstruction_of_ancestries.py
```

**Run DCH-ILP-PA for a single extant**:
```bash
python dch_ilp_pa.py -r clusters_for_table_4 -e nx_p=0.5_n=200 -g 1.15
```

**Run DCH-ILP-PA on batch**:
```bash
python run_dch_ilp_pa.py 
```