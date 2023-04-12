#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <queue>
#include <unordered_set>

#define INITIAL -1
#define SURVIVED 0
#define UNEXPLORED 1
#define DISCARDED 2

using namespace std;

struct Node
{
	int layer;
	int id;

	Node(int l, int i) : layer(l), id(i) { }

	friend bool operator < (const struct Node& n1, const struct Node& n2)
	{
		if (n1.layer > n2.layer)
			return true;
		else if (n1.layer == n2.layer && n1.id > n2.id)
			return true;
		else
			return false;
	}
};

class Subgraph
{
public:

	int k; // the coreness of this subgraph vertices except for layer 0
	int num_nodes;
	// vertex id -> the vertex sequence in the following vectors
	// including the vertices <= k
	unordered_map<double, int> vertIndex;

	vector<int> shell, layer;

	vector<vector<int> > adj, adj1, adj2; // all, lower or equal layer, higher layer
	
	vector<int> highKSupport;

	unordered_map<double, vector<int> > lowKNeighbors; // key: vertex id -> values: vertices' sequences

	vector<int> degreeBound, vertexStatus;

	void layerDecomposition();

	vector<double> findFollowers_D(double x_);

	int findFollowers_H(double x_);

	vector<double> collectFollowers(double x_);

	void shrink(int x, int u, unordered_set<int>* Survivals,vector<int>& degreeBound,vector<int>& vertexStatus);

	int findCollapsers_in(double x_);

	int findCollapsers_out(double x_);

	vector<double> collectCollapsers_in(double x_);

	vector<double> collectCollapsers_out(double x_);

	vector<double> collectFollowers_inE(double x_, double y_);

	vector<double> collectFollowers_outE(double x_);

	vector<double> collectCollapsers_inE(double x_, double y_);

	vector<double> collectCollapsers_outE(double x_);

	void shrink_E(int u);

};