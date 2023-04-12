#include "Monitoring.h"

Monitoring::Monitoring(Partition* partition)
{
	this->partition = partition;
}

void Monitoring::monitor()
{
	cout << "Output file: ";
	string output;
	cin >> output;
	output = "experiments/" + output;
	ofstream out(output);

	vector<unordered_map<Subgraph*, vector<double> > > Followers;
	vector<unordered_map<Subgraph*, vector<double> > > Collapsers;

	vector<double> empty_vector;
	unordered_map<Subgraph*, vector<double> > empty_map;
	for (double i = 0; i < partition->Followers.size(); i++)
	{
		Followers.push_back(empty_map);
		for (unordered_map<Subgraph*, int>::iterator it = partition->Followers[i].begin(); it != partition->Followers[i].end(); it++)
			Followers[i][it->first] = empty_vector;
	}

	for (double i = 0; i < partition->Collapsers.size(); i++)
	{
		Collapsers.push_back(empty_map);
		for (unordered_map<Subgraph*, int>::iterator it = partition->Collapsers[i].begin(); it != partition->Collapsers[i].end(); it++)
			Collapsers[i][it->first] = empty_vector;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
		{
			double v = it->first;
			int seqv = it->second;

			vector<double> followers = s->collectFollowers(v);
			Followers[v][s] = followers;

			if (s->shell[seqv] == s->k)
			{
				vector<double> collapsers = s->collectCollapsers_in(v);
				Collapsers[v][s] = collapsers;
			}
		}

		for (unordered_map<double, vector<int> >::iterator it = s->lowKNeighbors.begin(); it != s->lowKNeighbors.end(); it++)
		{
			double v = it->first;

			vector<double> collapsers = s->collectCollapsers_out(v);
			Collapsers[v][s] = collapsers;
		}
	}

	for (double i = 0; i < partition->datagraph->AdjList.size(); i++)
	{
		int u = (int)i;
		out << "[" << u << "] k" << partition->coreness[i] << endl;
		for (unordered_map<Subgraph*, vector<double> >::iterator it = Collapsers[i].begin(); it != Collapsers[i].end(); it++)
		{
			if ((it->second).size() == 0) continue;

			Subgraph* s = it->first;
			out << "C - " << s->k << ": ";
			for (int j = 0; j < (it->second).size(); j++)
				out << (int)(it->second)[j] << "; ";
			out << endl;
		}
		for (unordered_map<Subgraph*, vector<double> >::iterator it = Followers[i].begin(); it != Followers[i].end(); it++)
		{
			if ((it->second).size() == 0) continue;

			Subgraph* s = it->first;
			out << "F - " << s->k << ": ";
			for (int j = 0; j < (it->second).size(); j++)
				out << (int)(it->second)[j] << "; ";
			out << endl;
		}
		out << endl;
	}
}