#include "Master.h"

Master::Master(Partition* partition, Maintenance* maintenance)
{
	this->partition = partition;
	this->maintenance = maintenance;
}

void Master::Anchoring(int b)
{
	ofstream out("output.txt");
	
	for (int i = 0; i < b; i++)
	{
		double time1 = omp_get_wtime();

		double amount = 0;

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < partition->new_nodes.size(); i++)
		{
			Subgraph* s = partition->new_nodes[i];

			for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			{
				double v = it->first;

				int followers = s->findFollowers_H(v);
				partition->Followers[v][s] = followers;
			}

			amount += s->vertIndex.size();
		}

		double time2 = omp_get_wtime();

		int lambda = 0;
		double anchor = INITIAL;
		for (int i = 0; i < partition->datagraph->AdjList.size(); i++)
		{
			double v = i;

			int f = 0;
			for (unordered_map<Subgraph*, int>::iterator it = partition->Followers[v].begin(); it != partition->Followers[v].end(); it++)
				f += it->second;

			if (f > lambda)
			{
				lambda = f;
				anchor = v;
			}
		}
		double time3 = omp_get_wtime();
		partition->P_Maintenance(anchor);

		double time4 = omp_get_wtime();
		cout << i + 1 << ": " << anchor << "'s followers = "<< lambda << "; amount = " << amount << "; t1~t2: " << (time2 - time1) << "; t2~t3: " << (time3 - time2) <<  "; t3~t4: " << (time4 - time3) << endl;
		out << i + 1 << ": " << anchor << "'s followers = " << lambda << "; amount = " << amount << "; t1~t2: " << (time2 - time1) << "; t2~t3: " << (time3 - time2) << "; t3~t4: " << (time4 - time3) << endl;
	}
}

void Master::dynamicMaintain()
{
	double time1 = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
		{
			double v = it->first;
			int seqv = it->second;

			int followers = s->findFollowers_H(v);
			partition->Followers[v][s] = followers;

			if (s->shell[seqv] == s->k)
			{
				int collapsers = s->findCollapsers_in(v);
				partition->Collapsers[v][s] = collapsers;
			}
		}

		for (unordered_map<double, vector<int> >::iterator it = s->lowKNeighbors.begin(); it != s->lowKNeighbors.end(); it++)
		{
			double v = it->first;

			int collapsers = s->findCollapsers_out(v);
			partition->Collapsers[v][s] = collapsers;
		}
	}
	double time2 = omp_get_wtime();

	cout << "Finished followers and collapsers computation; " << "Time cost: " << (time2 - time1) << "s." << endl;

	/*for (int i = 0; i < partition->Collapsers.size(); i += 100)
	{
		int u = i;
		int collapsers = 0;
		for (unordered_map<Subgraph*, int>::iterator it = partition->Collapsers[u].begin(); it != partition->Collapsers[u].end(); it++)
			collapsers += it->second;

		cout << u << "-" << partition->datagraph->seq2id[u] << ":" << collapsers << endl;
	}*/

	double time3 = omp_get_wtime();

	for (int i = 0; i < 100; i++)
	{
		random_device rd1;
		uniform_int_distribution<int> dist1(0, partition->datagraph->AdjList.size() - 1);
		double _src = (double)dist1(rd1);
		double src = partition->datagraph->seq2id[_src];

		random_device rd2;
		uniform_int_distribution<int> dist2(0, partition->datagraph->AdjList[_src].size() - 1);
		double _dst = partition->datagraph->AdjList[_src][dist2(rd2)];
		double dst = partition->datagraph->seq2id[_dst];

		cout << src << " " << dst << endl;

		//vector<double> changed_verts;
		//maintenance->removeEdge(src, dst, &changed_verts);
		maintenance->removeEdge(src, dst,1);
	}
	
	for (int i = 0; i < 100; i++)
	{
		random_device rd1, rd2;
		uniform_int_distribution<int> dist(0, partition->datagraph->AdjList.size() - 1);
		double _src = (double)dist(rd1);
		double _dst = (double)dist(rd2);
		double src = partition->datagraph->seq2id[_src];
		double dst = partition->datagraph->seq2id[_dst];

		cout << src << " " << dst << endl;
		
		//vector<double> changed_verts;
		//maintenance->insertEdge(src, dst, &changed_verts);
		maintenance->insertEdge(src, dst,1);
	}

	double time4 = omp_get_wtime();

	cout << "Maintaining time: " << (time4 - time3) << "s. " << endl;
}
