#pragma once

#include <omp.h>
#include <fstream>
#include <map>
#include <set>
#include <random>
#include "Partition.h"
#include "Maintenance.h"
#include "DataGraph.h"

bool comparePowerScale(const double& a, const double& b);
bool compareRemoveTime(const int& a, const int& b);
bool compareInsertTime(const int& a, const int& b);
bool compareRemoveUpdate(const int& a, const int& b);
bool compareInsertUpdate(const int& a, const int& b);

class Experiment
{
public:

	Partition* partition;

	Maintenance* maintenance;

	DataGraph* datagraph;

	static vector<double> power_scale; // collapser_power / anchor_power;

	static vector<double> removing_time;
	static vector<double> inserting_time;

	static vector<double> changed_percent_rmv;
	static vector<double> changed_percent_ist;

	Experiment(Partition* partition, Maintenance* maintenance, DataGraph* datagraph);

	void effectiveness();

	void scalability();

	void streaming();

	void shellsize();

};