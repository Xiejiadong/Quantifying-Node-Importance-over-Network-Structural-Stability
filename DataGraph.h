#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <iostream>

using namespace std;

class DataGraph
{
public:

	/*
	* dimension 1: vertex ID := iterator sequence
	* dimension 2: the list of adjacent verticess
	*/
	vector<vector<double> > AdjList;
	unordered_map<double, double> id2seq; // map the vertex from ordinary dataset id to AdjList sequence
	unordered_map<double, double> seq2id; // map the sequence back to the id for presentation convenience

	DataGraph(string fileName,int &N,int &M);
};
