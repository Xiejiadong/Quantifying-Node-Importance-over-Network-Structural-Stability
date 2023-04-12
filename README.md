# Quantifying Node Importance over Network Structural Stability

This is the source code of our paper "Quantifying Node Importance over Network Structural Stability", Fan Zhang, Qingyuan Linghu, Jiadong Xie, Kai Wang, Xuemin Lin, Wenjie Zhang, accepted by KDD 2023.

# Experimental Environment

The experiments are performed on a CentOS Linux serve (Release 7.5.1804) with Quad-Core Intel Xeon CPU (E5-2640 v4 @ 2.20GHz) and 128G memory. All the algorithms are implemented in C++. The source code is compiled by GCC (7.3.0) under the -O3 optimization. We adopt OpenMP to utilize the multithreads of the machine.

# Compile

```shell
g++ -o main -std=c++11 -fopenmp DataGraph.cpp Subgraph.cpp Combination.cpp Partition.cpp Maintenance.cpp Experiment.cpp Master.cpp Monitoring.cpp main.cpp -O3
```
# Run

## DataGraph

Save datasets in /DataGraph, e.g., sample_graph.txt

## Run code

```shell
./main
```
## Output

The experiment results will be saved in `./experiments`.

### Effectiveness

Experiments of effectiveness in the paper.

```shell
dataset: sample_graph.txt
vertices = 58228 ; edges = 214078
Loaded dataset successfully!
index building: xxxxxs
Output file: sample_effectiveness.txt
```

### Static

Evaluation of static algorithm.

```shell
dataset: sample_graph.txt
vertices = 58228 ; edges = 214078
index building: xxxxxs
Output file: sample_static.txt
Time cost: xxxxxs
```

### Maintenance

Evaluation of maintenance algorithm.

```shell
dataset: sample_graph.txt
vertices = 58228 ; edges = 214078
index building: xxxxxs
Output file: sample_maintenance-2.txt
Number of threads: 2
Time cost: xxxxxs
```

