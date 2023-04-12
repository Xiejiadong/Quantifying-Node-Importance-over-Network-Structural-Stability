#include <iostream> 
#include <ctime>
#include <string>
#include <set>

#include "omp.h"
#include "DataGraph.h"
#include "Subgraph.h"
#include "Partition.h"
#include "Master.h"
#include "Maintenance.h"
#include "Combination.h"
#include "Experiment.h"
#include "Monitoring.h"

using namespace std;

int main(int argc, char* argv[])
{
	cout << "dataset: ";
	string dataset;
	getline(cin, dataset);
	int n,m;
	DataGraph datagraph(dataset,n,m);
	cout << "Loaded dataset successfully!" << endl;
	Partition partition(&datagraph);
	Maintenance maintenance(&datagraph, &partition);
	Master master(&partition, &maintenance);
	Experiment experiment(&partition, &maintenance, &datagraph);
	Monitoring monitoring(&partition);
	double time1 = omp_get_wtime();
	partition.P_Decomposition();
	double time2 = omp_get_wtime();
	cout << "index building: " << (time2 - time1) << "s" << endl;


	/* Please run different experiments by commenting/uncommenting the following lines */

	/*
	* Run different parts of current experiments by different functions in the class 'Experiment'
	*/
	experiment.effectiveness();
	experiment.scalability();
	experiment.streaming();
	return 0;
}
