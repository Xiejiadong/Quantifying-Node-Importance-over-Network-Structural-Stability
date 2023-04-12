#pragma once

#include "DataGraph.h"
#include "Partition.h"

using namespace std;

class Maintenance
{
public:

	DataGraph* datagraph;

	Partition* partition;

	Maintenance(DataGraph* datagraph, Partition* partition);

	void insertEdge(double src, double dst, unordered_set<Subgraph*>* deleted_nodes, int num_threads);

	void removeEdge(double src, double dst, unordered_set<Subgraph*>* deleted_nodes, int num_threads);

	void insertEdge(double src, double dst, int num_threads);

	void removeEdge(double src, double dst, int num_threads);
};