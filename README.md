# WSCAN
Instruction for experiment

WSCAN is implemented in C++ and compiled via the provided `Makefile`.

After cloning the repository, simply run:
```
make
```

This will generate an executable file named:
```
wscan
```

You can then execute it with the desired parameters.

---

## Parameters for Experiment

- eps : Similarity threshold parameter.

- mu : Minimum number of similar neighbors required.

- gamma : Balance parameter (used by the `WSCAN++` family; ignored by the others).

- similarity : Name of the similarity function to use. Available options:
  - `SCAN`         — structural similarity (unweighted)
  - `WSCAN`        — weighted structural similarity
  - `WSCAN++`      — proposed weighted similarity (min aggregation)
  - `WSCAN++_max`  — WSCAN++ variant with max aggregation
  - `WSCAN++_avg`  — WSCAN++ variant with average aggregation

- network : Path of the dataset folder (relative to `dataset/`). The folder has to contain a file named `network.dat`.

---

### Edge Perturbation Parameters

- edge_p : Percentage of edges to perturb.

- delta_p : Percentage for delta (delta is the parameter for perturbing edge weight).

- weight_method : Method to compute edge weight.
  Options:
  - "avg"
  - "max"

- seed : RNG seed for edge perturbation (default `42`). Must be `>= 1`.
  Only affects runs where perturbation is active (`edge_p > 0` and `delta_p > 0`).

---

### Parallel Execution

- parallel : Enable parallel mode.
  (Boolean flag — no value required. Just include `--parallel`.)

- threads : Number of threads to use in parallel mode. `0` uses the OpenMP default.

---

### Synthetic Mode

- synthetic : Enable synthetic dataset mode.
  (Boolean flag — no value required. Just include `--synthetic`.)

---

### Output

- output_file : Name (without extension) of the CSV file where results will be stored.
  Results are written to `output/<output_file>.csv` (the `output/` directory must exist).
  Each run appends one row; the header is written once.

  Columns:
  ```
  dataset, similarity, eps, mu, gamma, edge_p, delta_p, weight_method,
  is_parallel, num_threads, similarity_time, runtime_sec,
  ARI, NMI, modularity, conductance,
  num_clusters, num_hubs, num_outliers,
  clustered_node_ratio, hub_ratio, outlier_ratio
  ```
  In synthetic mode, three columns `mut, muw, idx` are inserted right after `dataset`.
  `ARI`/`NMI` are written as `NULL` when no ground truth is available, and all
  quality metrics are `NULL` when no cluster is found.

---

## Examples

### Basic Usage
```
./wscan --network karate --similarity WSCAN++ --eps 0.2 --mu 8 --gamma 0.7 --output_file result
```

---

### With Edge Perturbation
```
./wscan --network sociopattern_primaryschool \
        --similarity WSCAN++ \
        --eps 0.2 --mu 8 --gamma 0.7 \
        --edge_p 0.1 --delta_p 0.05 --weight_method max --seed 1 --output_file result
```

---

### Parallel Execution
```
./wscan --network sociopattern_primaryschool \
        --similarity WSCAN++ \
        --eps 0.2 --mu 8 --gamma 0.7 \
        --parallel --threads 4 --output_file result
```

---

### Synthetic Mode with Custom Output File
```
./wscan --network synthetic_dataset \
        --eps 0.2 --mu 8 --gamma 0.7 \
        --synthetic \
        --output_file synthetic_result
```

---

## Notes

1. Ground-truth Detection
   The program automatically checks whether the dataset has ground-truth labels by verifying the existence of a `labels.dat` file in the dataset folder (if synthetic, `community.dat`). No additional argument is required. Datasets without a label file are evaluated with internal metrics only (ARI/NMI are `NULL`).

2. Node Indexing
   Node indices must start from 1.
   Although non-consecutive numbering does not cause runtime errors, memory allocation is based on the maximum node ID.
   Therefore, non-sequential numbering may result in unnecessary memory usage.

3. Undirected Graph Assumption
   The dataset is currently processed in undirected mode.
   Therefore, duplicate edges with reversed node order must not appear in the input.
   For example, the following two edges must not coexist:

   1 2 5
   2 1 5

   Only one of them should be included in the dataset.
