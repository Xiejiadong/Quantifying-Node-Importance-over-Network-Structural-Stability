#include "Experiment.h"

vector<double> Experiment::power_scale;
vector<double> Experiment::removing_time;
vector<double> Experiment::inserting_time;
vector<double> Experiment::changed_percent_rmv;
vector<double> Experiment::changed_percent_ist;

bool comparePowerScale(const double& a, const double& b)
{
	return Experiment::power_scale[a] < Experiment::power_scale[b];
}

bool compareRemoveTime(const int& a, const int& b)
{
	return Experiment::removing_time[a] < Experiment::removing_time[b];
}

bool compareInsertTime(const int& a, const int& b)
{
	return Experiment::inserting_time[a] < Experiment::inserting_time[b];
}

bool compareRemoveUpdate(const int& a, const int& b)
{
	return Experiment::changed_percent_rmv[a] < Experiment::changed_percent_rmv[b];
}

bool compareInsertUpdate(const int& a, const int& b)
{
	return Experiment::changed_percent_ist[a] < Experiment::changed_percent_ist[b];
}

Experiment::Experiment(Partition* partition, Maintenance* maintenance, DataGraph* datagraph)
{
	this->partition = partition;
	this->maintenance = maintenance;
	this->datagraph = datagraph;
}

void Experiment::effectiveness()
{

	cout << "Output file: ";
	string output;
	cin >> output;
	output = "experiments/" + output;
	ofstream out(output);

	double time1 = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
		{
			double v = it->first;
			int seqv = it->second;

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
		}
	}
	double time3 = omp_get_wtime();

	out << "Collapse time: " << (time2 - time1) << endl;
	out << "Anchor time: " << (time3 - time2) << endl;

	map<int, pair<int, double> > degree_dist_collapsers, coreness_dist_collapsers, degree_dist_followers, coreness_dist_followers;
	
	vector<double> power_order;
	vector<double> collapser_power, anchor_power;
	for (double i = 0; i < partition->datagraph->AdjList.size(); i++)
	{
		double u = i;
		power_order.push_back(u);

		int degree_u = partition->datagraph->AdjList[u].size();
		int coreness_u = partition->coreness[u];

		double collapsers_u = 0;
		for (unordered_map<Subgraph*, int>::iterator it = partition->Collapsers[u].begin(); it != partition->Collapsers[u].end(); it++)
			collapsers_u += it->second;
		collapser_power.push_back(collapsers_u);

		double followers_u = 0;
		for (unordered_map<Subgraph*, int>::iterator it = partition->Followers[u].begin(); it != partition->Followers[u].end(); it++)
			followers_u += it->second;
		anchor_power.push_back(followers_u);

		if (followers_u == 0)
			if (collapsers_u == 0)
				power_scale.push_back(100);
			else
				power_scale.push_back(999999999);
		else
			power_scale.push_back(collapsers_u / followers_u * 100);

		if (degree_dist_collapsers.find(degree_u) == degree_dist_collapsers.end())
		{
			degree_dist_collapsers[degree_u] = make_pair(0, 0);
			degree_dist_followers[degree_u] = make_pair(0, 0);
		}

		int collapsers_num = degree_dist_collapsers[degree_u].first;
		double collapsers_avg = degree_dist_collapsers[degree_u].second;

		collapsers_avg = collapsers_avg * collapsers_num + collapsers_u;
		collapsers_avg /= (++collapsers_num);

		degree_dist_collapsers[degree_u].first = collapsers_num;
		degree_dist_collapsers[degree_u].second = collapsers_avg;

		int followers_num = degree_dist_followers[degree_u].first;
		double followers_avg = degree_dist_followers[degree_u].second;

		followers_avg = followers_avg * followers_num + followers_u;
		followers_avg /= (++followers_num);

		degree_dist_followers[degree_u].first = followers_num;
		degree_dist_followers[degree_u].second = followers_avg;

		if (coreness_dist_collapsers.find(coreness_u) == coreness_dist_collapsers.end())
		{
			coreness_dist_collapsers[coreness_u] = make_pair(0, 0);
			coreness_dist_followers[coreness_u] = make_pair(0, 0);
		}

		collapsers_num = coreness_dist_collapsers[coreness_u].first;
		collapsers_avg = coreness_dist_collapsers[coreness_u].second;

		collapsers_avg = collapsers_avg * collapsers_num + collapsers_u;
		collapsers_avg /= (++collapsers_num);

		coreness_dist_collapsers[coreness_u].first = collapsers_num;
		coreness_dist_collapsers[coreness_u].second = collapsers_avg;

		followers_num = coreness_dist_followers[coreness_u].first;
		followers_avg = coreness_dist_followers[coreness_u].second;

		followers_avg = followers_avg * followers_num + followers_u;
		followers_avg /= (++followers_num);

		coreness_dist_followers[coreness_u].first = followers_num;
		coreness_dist_followers[coreness_u].second = followers_avg;
	}

	out << endl << "Degree_distribution_Collapsers: " << endl;
	for (map<int, pair<int, double> >::iterator it = degree_dist_collapsers.begin(); it != degree_dist_collapsers.end(); it++)
	{
		int degree = it->first;
		int vert_num = (it->second).first;
		double collapsers_avg = (it->second).second;

		out << degree << "\t" << vert_num << "\t" << collapsers_avg << endl;
	}

	out << endl << "Degree_distribution_Followers: " << endl;
	for (map<int, pair<int, double> >::iterator it = degree_dist_followers.begin(); it != degree_dist_followers.end(); it++)
	{
		int degree = it->first;
		int vert_num = (it->second).first;
		double followers_avg = (it->second).second;

		out << degree << "\t" << vert_num << "\t" << followers_avg << endl;
	}

	out << endl << "Coreness_distribution_Collapsers: " << endl;
	for (map<int, pair<int, double> >::iterator it = coreness_dist_collapsers.begin(); it != coreness_dist_collapsers.end(); it++)
	{
		int coreness = it->first;
		int vert_num = (it->second).first;
		double collapsers_avg = (it->second).second;

		out << coreness << "\t" << vert_num << "\t" << collapsers_avg << endl;
	}

	out << endl << "Coreness_distribution_Followers: " << endl;
	for (map<int, pair<int, double> >::iterator it = coreness_dist_followers.begin(); it != coreness_dist_followers.end(); it++)
	{
		int coreness = it->first;
		int vert_num = (it->second).first;
		double followers_avg = (it->second).second;

		out << coreness << "\t" << vert_num << "\t" << followers_avg << endl;
	}

	sort(power_order.begin(), power_order.end(), comparePowerScale);

	out << endl << "Two-power Scale: " << endl;
	for (double i = 0; i < power_order.size(); i++)
	{
		double u = power_order[i];
		out << collapser_power[u] << "\t" << anchor_power[u] << "\t" << power_scale[u] << endl;
	}
}

void Experiment::scalability()
{
	cout << "Output file: ";
	string output;
	cin >> output;
	output = "experiments/" + output;
	ofstream out(output);
	double time1 = omp_get_wtime();
	
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(4)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			if (s->shell[seqv] == s->k)
			{
				int collapsers = s->findCollapsers_in(v);
				partition->Collapsers[v][s] = collapsers;
			}
		}

		tmp.clear();
		for (unordered_map<double, vector<int> >::iterator it = s->lowKNeighbors.begin(); it != s->lowKNeighbors.end(); it++)
			tmp.push_back({it->first,-1});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(4)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;

			int collapsers = s->findCollapsers_out(v);
			partition->Collapsers[v][s] = collapsers;
		}
	}


	double time2 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(4)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			int followers = s->findFollowers_H(v);
			partition->Followers[v][s] = followers;
		}
	}

	double time3 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(2)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			if (s->shell[seqv] == s->k)
			{
				int collapsers = s->findCollapsers_in(v);
				partition->Collapsers[v][s] = collapsers;
			}
		}

		tmp.clear();
		for (unordered_map<double, vector<int> >::iterator it = s->lowKNeighbors.begin(); it != s->lowKNeighbors.end(); it++)
			tmp.push_back({it->first,-1});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(2)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;

			int collapsers = s->findCollapsers_out(v);
			partition->Collapsers[v][s] = collapsers;
		}
	}
	double time4 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(2)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			int followers = s->findFollowers_H(v);
			partition->Followers[v][s] = followers;
		}
	}

	double time5 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
		{
			double v = it->first;
			int seqv = it->second;

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
	double time6 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
		{
			double v = it->first;
			int seqv = it->second;

			int followers = s->findFollowers_H(v);
			partition->Followers[v][s] = followers;
		}
	}
	double time7 = omp_get_wtime();
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(8)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			if (s->shell[seqv] == s->k)
			{
				int collapsers = s->findCollapsers_in(v);
				partition->Collapsers[v][s] = collapsers;
			}
		}

		tmp.clear();
		for (unordered_map<double, vector<int> >::iterator it = s->lowKNeighbors.begin(); it != s->lowKNeighbors.end(); it++)
			tmp.push_back({it->first,-1});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(8)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;

			int collapsers = s->findCollapsers_out(v);
			partition->Collapsers[v][s] = collapsers;
		}
	}
	double time8 = omp_get_wtime();
	
	for (int i = 0; i < partition->new_nodes.size(); i++)
	{
		Subgraph* s = partition->new_nodes[i];

		vector<pair<double,int>> tmp;
		int tmp_len;
		for (unordered_map<double, int>::iterator it = s->vertIndex.begin(); it != s->vertIndex.end(); it++)
			tmp.push_back({it->first,it->second});
		tmp_len=tmp.size();
		#pragma omp parallel for num_threads(8)
		for (int i=0;i<tmp_len;i++)
		{
			double v = tmp[i].first;
			int seqv = tmp[i].second;

			int followers = s->findFollowers_H(v);
			partition->Followers[v][s] = followers;
		}
	}
	double time9 = omp_get_wtime();

	out << " - 8 threads - " << endl;
	out << "Collapse time: " << (time8 - time7) << endl;
	out << "Anchor time: " << (time9 - time8) << endl;
	out << endl;

	out << " - 4 threads - " << endl;
	out << "Collapse time: " << (time2 - time1) << endl;
	out << "Anchor time: " << (time3 - time2) << endl;
	out << endl;

	out << " - 2 threads - " << endl;
	out << "Collapse time: " << (time4 - time3) << endl;
	out << "Anchor time: " << (time5 - time4) << endl;
	out << endl;

	out << " - 1 thread - " << endl;
	out << "Collapse time: " << (time6 - time5) << endl;
	out << "Anchor time: " << (time7 - time6) << endl;
	out << endl;
}

void Experiment::streaming()
{
	cout << "Output file: ";
	string output;
	cin >> output;
	output = "experiments/" + output;
	ofstream out(output);
	cout << "Number of threads: ";
	int num_threads;
	cin >> num_threads;

	
	double time1 = omp_get_wtime();
	#pragma omp parallel for num_threads(num_threads)
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

	out << "Finished followers and collapsers computation; " << "Time cost: " << (time2 - time1) << "s." << endl << endl;

	
	double avg_time_rmv = 0;
	double avg_update_rmv = 0;
	vector<pair<double,double>> remove_edges;

	mt19937 rd1(2022), rd2(2023);
	for (int i = 0; i < 100; i++)
	{
		uniform_int_distribution<int> dist1(0, partition->datagraph->AdjList.size() - 1);
		double _src = (double)dist1(rd1);
		double src = partition->datagraph->seq2id[_src];

		uniform_int_distribution<int> dist2(0, partition->datagraph->AdjList[_src].size() - 1);
		double _dst = partition->datagraph->AdjList[_src][dist2(rd2)];
		double dst = partition->datagraph->seq2id[_dst];

		remove_edges.push_back({src,dst});
		double time3 = omp_get_wtime();

		unordered_set<Subgraph*> deleted_nodes;
		maintenance->removeEdge(src, dst, &deleted_nodes, num_threads);

		double time4 = omp_get_wtime();

		removing_time.push_back(time4 - time3);
		avg_time_rmv += (time4 - time3);

		unordered_set<double> changed_verts_set;

		for (unordered_set<Subgraph*>::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); it++)
		{
			Subgraph* ss = *it;
			for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				changed_verts_set.insert(itt->first);

			for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				changed_verts_set.insert(itt->first);
		}

		for (int i = 0; i < partition->new_nodes.size(); i++)
		{
			Subgraph* ss = partition->new_nodes[i];
			for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				changed_verts_set.insert(it->first);

			for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				changed_verts_set.insert(it->first);
		}

		int change_vert_num = changed_verts_set.size() > 0 ? changed_verts_set.size() : 2;
		changed_percent_rmv.push_back(change_vert_num);
		avg_update_rmv += change_vert_num;
	}
	avg_time_rmv /= 100;
	avg_update_rmv /= 100;

	double avg_time_ist = 0;
	double avg_update_ist = 0;
	for (auto i:remove_edges)
	{
		uniform_int_distribution<int> dist(0, partition->datagraph->AdjList.size() - 1);
		double _src = (double)dist(rd1);
		double _dst = (double)dist(rd2);
		double src = i.first;
		double dst = i.second;
		double time3 = omp_get_wtime();

		unordered_set<Subgraph*> deleted_nodes;
		maintenance->insertEdge(src, dst, &deleted_nodes, num_threads);

		double time4 = omp_get_wtime();

		inserting_time.push_back(time4 - time3);
		avg_time_ist += (time4 - time3);

		unordered_set<double> changed_verts_set;

		for (unordered_set<Subgraph*>::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); it++)
		{
			Subgraph* ss = *it;
			for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				changed_verts_set.insert(itt->first);

			for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				changed_verts_set.insert(itt->first);
		}

		for (int i = 0; i < partition->new_nodes.size(); i++)
		{
			Subgraph* ss = partition->new_nodes[i];
			for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				changed_verts_set.insert(it->first);

			for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				changed_verts_set.insert(it->first);
		}

		int change_vert_num = changed_verts_set.size() > 0 ? changed_verts_set.size() : 2;
		changed_percent_ist.push_back(change_vert_num);
		avg_update_ist += change_vert_num;
	}
	avg_time_ist /= 100;
	avg_update_ist /= 100;

	vector<int> rmv_time_order;
	vector<int> ist_time_order;
	vector<int> rmv_update_order;
	vector<int> ist_update_order;
	for (int i = 0; i < 100; i++)
	{
		rmv_time_order.push_back(i);
		ist_time_order.push_back(i);
		rmv_update_order.push_back(i);
		ist_update_order.push_back(i);
	}
	
	sort(rmv_time_order.begin(), rmv_time_order.end(), compareRemoveTime);
	sort(ist_time_order.begin(), ist_time_order.end(), compareInsertTime);
	sort(rmv_update_order.begin(), rmv_update_order.end(), compareRemoveUpdate);
	sort(ist_update_order.begin(), ist_update_order.end(), compareInsertUpdate);
	
	int min = 0;
	int mid = 50;
	int max = 99;

	double min_time_rmv = removing_time[rmv_time_order[min]];
	double mid_time_rmv = removing_time[rmv_time_order[mid]];
	double max_time_rmv = removing_time[rmv_time_order[max]];
	
	double min_time_ist = inserting_time[ist_time_order[min]];
	double mid_time_ist = inserting_time[ist_time_order[mid]];
	double max_time_ist = inserting_time[ist_time_order[max]];

	double min_update_rmv = changed_percent_rmv[rmv_update_order[min]];
	double mid_update_rmv = changed_percent_rmv[rmv_update_order[mid]];
	double max_update_rmv = changed_percent_rmv[rmv_update_order[max]];

	double min_update_ist = changed_percent_ist[ist_update_order[min]];
	double mid_update_ist = changed_percent_ist[ist_update_order[mid]];
	double max_update_ist = changed_percent_ist[ist_update_order[max]];
	
	out << "- Removing time -" << endl;
	out << "Avg: " << avg_time_rmv << endl;
	out << "Min: " << min_time_rmv << endl;
	out << "Mid: " << mid_time_rmv << endl;
	out << "Max: " << max_time_rmv << endl;
	out << endl;

	out << "- Inserting time -" << endl;
	out << "Avg: " << avg_time_ist << endl;
	out << "Min: " << min_time_ist << endl;
	out << "Mid: " << mid_time_ist << endl;
	out << "Max: " << max_time_ist << endl;
	out << endl;

	out << "- Removing update -" << endl;
	out << "Avg: " << avg_update_rmv << endl;
	out << "Min: " << min_update_rmv << endl;
	out << "Mid: " << mid_update_rmv << endl;
	out << "Max: " << max_update_rmv << endl;
	out << endl;

	out << "- Inserting update -" << endl;
	out << "Avg: " << avg_update_ist << endl;
	out << "Min: " << min_update_ist << endl;
	out << "Mid: " << mid_update_ist << endl;
	out << "Max: " << max_update_ist << endl;
	out << endl;
	
	out << "Removing time (1 - 100):" << endl;
	for (int i = 0; i < 100; i++)
		out << removing_time[i] << endl;
	out << endl;

	out << "Inserting time (1 - 100):" << endl;
	for (int i = 0; i < 100; i++)
		out << inserting_time[i] << endl;
	out << endl;

	out << "Updated vertices percentage (removal; 1 - 100):" << endl;
	for (int i = 0; i < 100; i++)
		out << changed_percent_rmv[i] << endl;
	out << endl;

	out << "Updated vertices percentage (insertion; 1 - 100):" << endl;
	for (int i = 0; i < 100; i++)
		out << changed_percent_ist[i] << endl;
	out << endl;
}
