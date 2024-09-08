# DCH-ILP-PA (Divide  and  Conquer  Heuristic  using the ILP)

## Description
**DCH-ILP-PA** is a parallel heuristic designed specifically for the large PPI network reconstruction.

## Dependency
```
Python >= 3.6

gurobi >= 8.0

igraph >= 0.7.1

networkx >= 2.4

leidenalg >= 0.7.0
```

## Usage
Running on cluster is strongly recommended as the ILP algorithm is computationally expensive. Recommended running on a machine with at least `32GB` of RAM.

* `ILP.py` and `run_ILP_for_sim.py` are the key part of the DCH-ILP-PA algorithm. One can run the pre-clustered sample graph (with `n_nodes = 50`, `q_con = 0.4`, `q_mod = 0.7`) using

```{bash}
python3 run_ILP_for_sim.py
```

* `benchmark_lkl.ipynb` is the jupyter-notebook that benchmark the result(based on likelihood value) obtained from **DCH-ILP** algorithm and **ReverseDMC(Greedy)** algorithm. One can run `benchmark_lkl.ipynb` directly with pre-runed **DCH-ILP** solution.

* `prerun_partitioning.ipynb` is the jupyter-notebook that pre-clusters the PPI network before running **DCH-ILP**. Current clustering result of the sample network is stored in `./Example_data/nG50/`, sub-graph of size no larger than $11$ is recommended.

* `simulation.ipynb` simulate the PPI network given network size `n_nodes`, DMC parameters `q_con` and `q_mod`.

* `./Example_data/nG50/` stores the sample network, clustering result and the pre-runed solution of the network.






