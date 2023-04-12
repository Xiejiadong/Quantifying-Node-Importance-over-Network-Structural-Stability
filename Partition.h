#pragma once

#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <omp.h>

#include "DataGraph.h"
#include "Subgraph.h"

#define INITIAL -1

using namespace std;

bool compareDegree(const double& a, const double& b);

class Partition
{
public:

	DataGraph* datagraph;
	int num_nodes;

	unordered_set<Subgraph*> CC_Partitions;
	vector<Subgraph*> self_node; // the k-shell component of itself
	vector<Subgraph*> new_nodes;
	vector<unordered_map<Subgraph*, int> > Followers;
	vector<unordered_map<Subgraph*, int> > Collapsers;

	/*
	* vec_seq: vertext ID
	* vec_ele: core number of the vertex
	*/
	static vector<int> coreness; // VLDB2015's d

	/*
	* vec_seq: vertext ID
	* pair_1: k-shell number
	* pair_2: layer number of its shell
	*/
	vector<pair<int, int> > shell_layer;

	/*
	* vec_seq: ascending total order of (#shell, #layer)
	* vec_ele: vertex ID
	*/
	vector<double> total_order; // VLDB2015's D

	/*
	* vec_seq: vertex ID
	* vec_ele: vector position of the vertex in total_order
	*/
	vector<double> order_pointer; // VLDB2015's p

	/*
	* vec_seq: k-shell number
	* vec_ele: the start position of vertices having this #shell in total_order
	*/
	vector<double> shell_tag; // VLDB2015's b

	unordered_set<double> anchor_verts; // the set of chosen anchored vertices

	/*
	* vec_seq: vertex ID
	* vec_ele: if the vertex is allocated in DAG
	*/
	vector<bool> allocated;

	int max_degree; // auxiliary variable help adapt to VLDB2015's algorithm

	Partition(DataGraph* datagraph);

	void SL_Decomposition();

	void P_Decomposition();

	void P_Maintenance(double anchor); // only for followers reuse

	void shellConnect(double src, Subgraph* shell_node);

	void shellConnect(double src, Subgraph* shell_node, unordered_set<Subgraph*>* deleted_nodes);
};
