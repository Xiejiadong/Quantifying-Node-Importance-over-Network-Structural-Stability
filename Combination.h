#pragma once

#include <omp.h>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <fstream>

#include "DataGraph.h"

#define INITIAL -1

using namespace std;

bool compareDegree1(const int& a, const int& b);

bool comparePotential(const int& a, const int& b);

class Combination
{
public:

	DataGraph* datagraph;
	
	static vector<int> anch_degrees;
	static vector<double> vert_potentials;

	vector<int> coreness_init;

	vector<vector<int> > comb_verts; // each vector item is a comination of sampled vertices
	vector<vector<int> > v_comb_seqs; // stores the combinations' seqs in comb_verts for each vertex
	vector<vector<double> > v_comb_shares; // the share a vertex in a combination, following v_comb_seqs

	vector<int> v_samp_num;
	vector<double> v_share_avg;
	vector<double> v_share_var;

	Combination(DataGraph* datagraph);

	vector<int> coreDecomposition(vector<int>& anchors);

	void sampleRand(int samp_num, double br);

	void updStatistics(double r1, double r2, double threshold, double k_, double b_);

	int topbTest(double br);

	void trackSampling(string dataset, double br, int track_freq, double threshold, 
		double r1, double r2, double k_, double b_);

	void assignAnchors();
};