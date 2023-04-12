#include "Maintenance.h"

Maintenance::Maintenance(DataGraph* datagraph, Partition* partition)
{
	this->datagraph = datagraph;
	this->partition = partition;
}

void Maintenance::insertEdge(double src, double dst, unordered_set<Subgraph*>* deleted_nodes, int num_thr)
{
	if (src != dst && datagraph->id2seq.find(src) != datagraph->id2seq.end() && datagraph->id2seq.find(dst) != datagraph->id2seq.end())
	{
		double u = datagraph->id2seq[src];
		double v = datagraph->id2seq[dst];

		int flag1 = 0;
		for (double i = 0; i < datagraph->AdjList[u].size(); i++)
		{
			if (datagraph->AdjList[u][i] == v)
			{
				flag1 = 1;
				break;
			}
		}
		if (flag1 == 0) datagraph->AdjList[u].push_back(v);

		int flag2 = 0;
		for (double i = 0; i < datagraph->AdjList[v].size(); i++)
		{
			if (datagraph->AdjList[v][i] == u)
			{
				flag2 = 1;
				break;
			}
		}
		if (flag2 == 0) datagraph->AdjList[v].push_back(u);

		if (flag1 == 0 && flag2 == 0)
		{
			vector<double> followers;
			Subgraph* s;

			if (partition->coreness[u] < partition->coreness[v])
			{
				followers = partition->self_node[u]->collectFollowers_outE(u);
				s = partition->self_node[u];
				
				Subgraph* s_ = partition->self_node[v];

				if (s_->vertIndex.find(u) == s_->vertIndex.end())
				{
					int tmp = s_->vertIndex.size();
					s_->vertIndex[u] = tmp;
					s_->shell.push_back(partition->coreness[u]);
					s_->layer.push_back(INITIAL);
					s_->highKSupport.push_back(0);
					vector<int> Neigh;
					s_->adj1.push_back(Neigh);
					s_->adj2.push_back(Neigh);

				}

				s_->adj1[s_->vertIndex[v]].push_back(s_->vertIndex[u]);
				s_->adj2[s_->vertIndex[u]].push_back(s_->vertIndex[v]);
				partition->Followers[u][s_] = s_->findFollowers_H(u);
				// adding (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "+<: " << followers.size() << endl;
			}
			else if (partition->coreness[u] == partition->coreness[v])
			{
				followers = partition->self_node[u]->collectFollowers_inE(u, v);
				s = partition->self_node[u];
				// cout << "+=: " << followers.size() << endl;
			}
			else // partition->coreness[u] > partition->coreness[v]
			{
				followers = partition->self_node[v]->collectFollowers_outE(v);
				s = partition->self_node[v];
				
				Subgraph* s_ = partition->self_node[u];

				if (s_->vertIndex.find(v) == s_->vertIndex.end())
				{
					int tmp = s_->vertIndex.size();
					s_->vertIndex[v] = tmp;
					s_->shell.push_back(partition->coreness[v]);
					s_->layer.push_back(INITIAL);
					s_->highKSupport.push_back(0);

					vector<int> Neigh;
					s_->adj1.push_back(Neigh);
					s_->adj2.push_back(Neigh);

				}
				s_->adj1[s_->vertIndex[u]].push_back(s_->vertIndex[v]);
				s_->adj2[s_->vertIndex[v]].push_back(s_->vertIndex[u]);

				partition->Followers[v][s_] = s_->findFollowers_H(v);
				// adding (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "+>: " << followers.size() << endl;
			}

			for (int i = 0; i < followers.size(); i++)
			{
				double w = followers[i];
				partition->coreness[w]++;
				partition->shell_layer[w].first++;
				partition->shell_layer[w].second = INITIAL;
			}

			partition->new_nodes.clear();
			deleted_nodes->insert(s);

			for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
			{
				double w = itt->first;
				int seqw = itt->second;

				if (s->shell[seqw] == s->k)
				{
					int flag = 0;
					for (int i = 0; i < partition->new_nodes.size(); i++)
					{
						if (partition->self_node[w] == partition->new_nodes[i])
						{
							flag = 1;
							break;
						}
					}

					if (flag == 1) continue;

					Subgraph* shell_node = new Subgraph();
					partition->new_nodes.push_back(shell_node);

					partition->Followers[w][shell_node] = INITIAL;
					partition->Collapsers[w][shell_node] = INITIAL;
					partition->self_node[w] = shell_node;

					shell_node->k = partition->coreness[w];
					shell_node->shell.push_back(partition->coreness[w]);
					shell_node->layer.push_back(INITIAL);
					vector<int> Neigh;
					shell_node->adj.push_back(Neigh);
					shell_node->highKSupport.push_back(0);
					shell_node->vertIndex[w] = shell_node->highKSupport.size() - 1;
					partition->shellConnect(w, shell_node, deleted_nodes);
				}
			}

			for (unordered_set<Subgraph*>::iterator it = deleted_nodes->begin(); it != deleted_nodes->end(); it++)
			{
				Subgraph* ss = *it;
				partition->CC_Partitions.erase(ss);
				for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				{
					double w = itt->first;
					partition->Followers[w].erase(ss);
				}

				for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				{
					double w = itt->first;
					partition->Collapsers[w].erase(ss);
				}
			}

			// #pragma omp parallel for schedule(dynamic)
			vector<double> thred_time;
			for (int i=0;i<num_thr;i++)
				thred_time.push_back(0);
			// #pragma omp parallel for num_threads(num_thr)
			for (int i = 0; i < partition->new_nodes.size(); i++)
			{
				Subgraph* ss = partition->new_nodes[i];
				ss->layerDecomposition();
				partition->CC_Partitions.insert(ss);

				vector<pair<double,int>> tmp;
				int tmp_len;

				for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
					tmp.push_back({it->first,it->second});
				for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
					tmp.push_back({it->first,-1});
				tmp_len=tmp.size();


				
				#pragma omp parallel for num_threads(num_thr)
				for (int i=0;i<tmp_len;i++)
				// for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				{
					double beginTime=omp_get_wtime(),id=omp_get_thread_num();
					if (tmp[i].second==-1)
					{
						double w=tmp[i].first;
						int collapsers = ss->findCollapsers_out(w);
						partition->Collapsers[w][ss] = collapsers;
					}
					else{
						double w = tmp[i].first;
						int seqw = tmp[i].second;

						int followers = ss->findFollowers_H(w);
						partition->Followers[w][ss] = followers;

						if (ss->shell[seqw] == ss->k)
						{
							int collapsers = ss->findCollapsers_in(w);
							partition->Collapsers[w][ss] = collapsers;
						}
					}
					double endTime=omp_get_wtime();
					thred_time[id]+=endTime-beginTime;
					
				}

				// vector<double> tmp;
				// int tmp_len;
				// for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				// 	tmp.push_back(it->first);
				// tmp_len=tmp.size();
				// #pragma omp parallel for num_threads(num_threads)
				// for (int i=0;i<tmp_len;i++)
				// // for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				// {
				// 	// double w = it->first;
				// 	double w=tmp[i];

				// 	int collapsers = ss->findCollapsers_out(w);
				// 	partition->Collapsers[w][ss] = collapsers;
				// }
			}
			for (int i=0;i<num_thr;i++)
				// cout<<i<<" "<<thred_time[i]<<endl;
		}
		else
		{
			cout << "ERROR: no need to insert already existed edges!" << endl;
		}
	}
	else
	{
		cout << "ERROR: inserting edges between illegal vertices!" << endl;
	}
}

void Maintenance::removeEdge(double src, double dst, unordered_set<Subgraph*>* deleted_nodes, int num_thr)
{
	if (src != dst && datagraph->id2seq.find(src) != datagraph->id2seq.end() && datagraph->id2seq.find(dst) != datagraph->id2seq.end())
	{
		double u = datagraph->id2seq[src];
		double v = datagraph->id2seq[dst];

		if (datagraph->AdjList[u].size() == 1 || datagraph->AdjList[v].size() == 1)
		{
			cout << "ERROR: removing edges between illegal vertices!" << endl;
			return;
		}

		int flag1 = 0;
		for (double i = 0; i < datagraph->AdjList[u].size(); i++)
		{
			if (datagraph->AdjList[u][i] == v)
			{
				flag1 = 1;
				datagraph->AdjList[u].erase(datagraph->AdjList[u].begin() + i);
				break;
			}
		}

		int flag2 = 0;
		for (double i = 0; i < datagraph->AdjList[v].size(); i++)
		{
			if (datagraph->AdjList[v][i] == u)
			{
				flag2 = 1;
				datagraph->AdjList[v].erase(datagraph->AdjList[v].begin() + i);
				break;
			}
		}

		if (flag1 == 1 && flag2 == 1)
		{
			vector<double> collapsers;
			Subgraph* s;

			if (partition->coreness[u] < partition->coreness[v])
			{
				// if (src == 109 && dst == 4103) cout << u << "$" << v << endl;
				collapsers = partition->self_node[u]->collectCollapsers_outE(u);
				s = partition->self_node[u];
				
				Subgraph* s_ = partition->self_node[v];
					
				int u_seq = s_->vertIndex[u];
				//s_->vertIndex.erase(u);
				for (double i = 0; i < s_->adj2[u_seq].size(); i++)
				{
					int neigh = s_->adj2[u_seq][i];
					for (int j = 0; j < s_->adj1[neigh].size(); j++)
					{
						if (s_->adj1[neigh][j] == u_seq)
						{
							s_->adj1[neigh][j] = s_->adj1[neigh][s_->adj1[neigh].size() - 1];
							s_->adj1[neigh].pop_back();
								
							break;
						}
					}

					if (s_->adj.size() == 0)
						continue;

					for (int j = 0; j < s_->adj[neigh].size(); j++)
					{
						if (s_->adj[neigh][j] == u_seq)
						{
							s_->adj[neigh][j] = s_->adj[neigh][s_->adj[neigh].size() - 1];
							s_->adj[neigh].pop_back();
								
							break;
						}
					}
				}
				s_->adj2[u_seq].clear();
				if (s_->adj.size() != 0) s_->adj[u_seq].size();
				
				partition->Followers[u].erase(s_);
				// removing (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "-<: " << collapsers.size() << endl;
			}
			else if (partition->coreness[u] == partition->coreness[v])
			{
				collapsers = partition->self_node[u]->collectCollapsers_inE(u, v);
				s = partition->self_node[u];
				// cout << "-=: " << collapsers.size() << endl;
			}
			else // partition->coreness[u] > partition->coreness[v]
			{
				collapsers = partition->self_node[v]->collectCollapsers_outE(v);
				s = partition->self_node[v];

				Subgraph* s_ = partition->self_node[u];

				int v_seq = s_->vertIndex[v];
				//s_->vertIndex.erase(v);
				for (double i = 0; i < s_->adj2[v_seq].size(); i++)
				{
					int neigh = s_->adj2[v_seq][i];
					for (int j = 0; j < s_->adj1[neigh].size(); j++)
					{
						if (s_->adj1[neigh][j] == v_seq)
						{
							s_->adj1[neigh][j] = s_->adj1[neigh][s_->adj1[neigh].size() - 1];
							s_->adj1[neigh].pop_back();

							break;
						}
					}

					if (s_->adj.size() == 0)
						continue;

					for (int j = 0; j < s_->adj[neigh].size(); j++)
					{
						if (s_->adj[neigh][j] == v_seq)
						{
							s_->adj[neigh][j] = s_->adj[neigh][s_->adj[neigh].size() - 1];
							s_->adj[neigh].pop_back();

							break;
						}
					}
				}
				s_->adj2[v_seq].clear();
				if (s_->adj.size() != 0) s_->adj[v_seq].size();

				partition->Followers[v].erase(s_);
				// removing (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "->: " << collapsers.size() << endl;
			}
			
			for (int i = 0; i < collapsers.size(); i++)
			{
				double w = collapsers[i];
				partition->coreness[w]--;
				partition->shell_layer[w].first--;
				partition->shell_layer[w].second = INITIAL;
			}
			
			partition->new_nodes.clear();
			deleted_nodes->insert(s);
			
			for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
			{
				double w = itt->first;
				int seqw = itt->second;

				if (s->shell[seqw] == s->k)
				{
					int flag = 0;
					for (int i = 0; i < partition->new_nodes.size(); i++)
					{
						if (partition->self_node[w] == partition->new_nodes[i])
						{
							flag = 1;
							break;
						}
					}

					if (flag == 1) continue;

					Subgraph* shell_node = new Subgraph();
					partition->new_nodes.push_back(shell_node);

					partition->Followers[w][shell_node] = INITIAL;
					partition->Collapsers[w][shell_node] = INITIAL;
					partition->self_node[w] = shell_node;

					shell_node->k = partition->coreness[w];
					shell_node->shell.push_back(partition->coreness[w]);
					shell_node->layer.push_back(INITIAL);
					vector<int> Neigh;
					shell_node->adj.push_back(Neigh);
					shell_node->highKSupport.push_back(0);
					shell_node->vertIndex[w] = shell_node->highKSupport.size() - 1;
					partition->shellConnect(w, shell_node, deleted_nodes);
				}
			}
			
			for (unordered_set<Subgraph*>::iterator it = deleted_nodes->begin(); it != deleted_nodes->end(); it++)
			{
				Subgraph* ss = *it;
				partition->CC_Partitions.erase(ss);
				for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				{
					double w = itt->first;
					partition->Followers[w].erase(ss);
				}

				for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				{
					double w = itt->first;
					partition->Collapsers[w].erase(ss);
				}
			}
			
			// #pragma omp parallel for num_threads(num_threads)
			// #pragma omp parallel for num_threads(num_thr)
			for (int i = 0; i < partition->new_nodes.size(); i++)
			{
				Subgraph* ss = partition->new_nodes[i];
				ss->layerDecomposition();
				partition->CC_Partitions.insert(ss);
				
				vector<pair<double,int>> tmp;
				int tmp_len;

				for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
					tmp.push_back({it->first,it->second});
				for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
					tmp.push_back({it->first,-1});
				tmp_len=tmp.size();

				
				#pragma omp parallel for num_threads(num_thr)
				for (int i=0;i<tmp_len;i++)
				// for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				{
					if(tmp[i].second==-1)
					{
						double w = tmp[i].first;
						int collapsers = ss->findCollapsers_out(w);
						partition->Collapsers[w][ss] = collapsers;
					}
					else{
						double w = tmp[i].first;
						int seqw = tmp[i].second;

						int followers = ss->findFollowers_H(w);
						partition->Followers[w][ss] = followers;

						if (ss->shell[seqw] == ss->k)
						{
							int collapsers = ss->findCollapsers_in(w);
							partition->Collapsers[w][ss] = collapsers;
						}
					}
				}

				// vector<double> tmp;
				// int tmp_len;
				// for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				// 	tmp.push_back(it->first);
				// tmp_len=tmp.size();
				// #pragma omp parallel for num_threads(num_threads)
				// for (int i=0;i<tmp_len;i++)
				// // for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				// {
				// 	double w = tmp[i];

				// 	int collapsers = ss->findCollapsers_out(w);
				// 	partition->Collapsers[w][ss] = collapsers;
				// }
			}
		}
		else
		{
			cout << "ERROR: removing nonexistent edges!" << endl;
		}
	}
	else
	{
		cout << "ERROR: removing edges between illegal vertices!" << endl;
	}
}

void Maintenance::insertEdge(double src, double dst, int num_thr)
{
	if (src != dst && datagraph->id2seq.find(src) != datagraph->id2seq.end() && datagraph->id2seq.find(dst) != datagraph->id2seq.end())
	{
		double u = datagraph->id2seq[src];
		double v = datagraph->id2seq[dst];

		int flag1 = 0;
		for (double i = 0; i < datagraph->AdjList[u].size(); i++)
		{
			if (datagraph->AdjList[u][i] == v)
			{
				flag1 = 1;
				break;
			}
		}
		if (flag1 == 0) datagraph->AdjList[u].push_back(v);

		int flag2 = 0;
		for (double i = 0; i < datagraph->AdjList[v].size(); i++)
		{
			if (datagraph->AdjList[v][i] == u)
			{
				flag2 = 1;
				break;
			}
		}
		if (flag2 == 0) datagraph->AdjList[v].push_back(u);

		if (flag1 == 0 && flag2 == 0)
		{
			vector<double> followers;
			Subgraph* s;

			if (partition->coreness[u] < partition->coreness[v])
			{
				followers = partition->self_node[u]->collectFollowers_outE(u);
				s = partition->self_node[u];

				Subgraph* s_ = partition->self_node[v];
				if (s_->vertIndex.find(u) == s_->vertIndex.end())
				{
					int tmp = s_->vertIndex.size();
					s_->vertIndex[u] = tmp;
					s_->shell.push_back(partition->coreness[u]);
					s_->layer.push_back(INITIAL);
					s_->highKSupport.push_back(0);

					vector<int> Neigh;
					s_->adj1.push_back(Neigh);
					s_->adj2.push_back(Neigh);
				}
				s_->adj1[s_->vertIndex[v]].push_back(s_->vertIndex[u]);
				s_->adj2[s_->vertIndex[u]].push_back(s_->vertIndex[v]);

				partition->Followers[u][s_] = s_->findFollowers_H(u);
				// adding (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "+<: " << followers.size() << endl;
			}
			else if (partition->coreness[u] == partition->coreness[v])
			{
				followers = partition->self_node[u]->collectFollowers_inE(u, v);
				s = partition->self_node[u];
				// cout << "+=: " << followers.size() << endl;
			}
			else // partition->coreness[u] > partition->coreness[v]
			{
				followers = partition->self_node[v]->collectFollowers_outE(v);
				s = partition->self_node[v];

				Subgraph* s_ = partition->self_node[u];
				if (s_->vertIndex.find(v) == s_->vertIndex.end())
				{
					int tmp = s_->vertIndex.size();
					s_->vertIndex[v] = tmp;
					s_->shell.push_back(partition->coreness[v]);
					s_->layer.push_back(INITIAL);
					s_->highKSupport.push_back(0);

					vector<int> Neigh;
					s_->adj1.push_back(Neigh);
					s_->adj2.push_back(Neigh);
				}

				s_->adj1[s_->vertIndex[u]].push_back(s_->vertIndex[v]);
				s_->adj2[s_->vertIndex[v]].push_back(s_->vertIndex[u]);

				partition->Followers[v][s_] = s_->findFollowers_H(v);
				// adding (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "+>: " << followers.size() << endl;
			}

			for (int i = 0; i < followers.size(); i++)
			{
				double w = followers[i];
				partition->coreness[w]++;
				partition->shell_layer[w].first++;
				partition->shell_layer[w].second = INITIAL;
			}

			partition->new_nodes.clear();
			unordered_set<Subgraph*> deleted_nodes;
			deleted_nodes.insert(s);

			for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
			{
				double w = itt->first;
				int seqw = itt->second;

				if (s->shell[seqw] == s->k)
				{
					int flag = 0;
					for (int i = 0; i < partition->new_nodes.size(); i++)
					{
						if (partition->self_node[w] == partition->new_nodes[i])
						{
							flag = 1;
							break;
						}
					}

					if (flag == 1) continue;

					Subgraph* shell_node = new Subgraph();
					partition->new_nodes.push_back(shell_node);

					partition->Followers[w][shell_node] = INITIAL;
					partition->Collapsers[w][shell_node] = INITIAL;
					partition->self_node[w] = shell_node;

					shell_node->k = partition->coreness[w];
					shell_node->shell.push_back(partition->coreness[w]);
					shell_node->layer.push_back(INITIAL);
					vector<int> Neigh;
					shell_node->adj.push_back(Neigh);
					shell_node->highKSupport.push_back(0);
					shell_node->vertIndex[w] = shell_node->highKSupport.size() - 1;
					partition->shellConnect(w, shell_node, &deleted_nodes);
				}
			}

			for (unordered_set<Subgraph*>::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); it++)
			{
				Subgraph* ss = *it;
				partition->CC_Partitions.erase(ss);
				for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				{
					double w = itt->first;
					partition->Followers[w].erase(ss);
				}

				for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				{
					double w = itt->first;
					partition->Collapsers[w].erase(ss);
				}
			}

			// #pragma omp parallel for schedule(dynamic)
			#pragma omp parallel for schedule(dynamic, num_thr)
			for (int i = 0; i < partition->new_nodes.size(); i++)
			{
				Subgraph* ss = partition->new_nodes[i];
				ss->layerDecomposition();
				partition->CC_Partitions.insert(ss);
				// #pragma omp parallel for num_threads(num_threads)
				for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				{
					double w = it->first;
					int seqw = it->second;

					int followers = ss->findFollowers_H(w);
					partition->Followers[w][ss] = followers;

					if (ss->shell[seqw] == ss->k)
					{
						int collapsers = ss->findCollapsers_in(w);
						partition->Collapsers[w][ss] = collapsers;
					}
				}
				// #pragma omp parallel for num_threads(num_threads)
				for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				{
					double w = it->first;

					int collapsers = ss->findCollapsers_out(w);
					partition->Collapsers[w][ss] = collapsers;
				}
			}
		}
		else
		{
			cout << "ERROR: no need to insert already existed edges!" << endl;
		}
	}
	else
	{
		cout << "ERROR: inserting edges between illegal vertices!" << endl;
	}
}

void Maintenance::removeEdge(double src, double dst, int num_thr)
{
	if (src != dst && datagraph->id2seq.find(src) != datagraph->id2seq.end() && datagraph->id2seq.find(dst) != datagraph->id2seq.end())
	{
		double u = datagraph->id2seq[src];
		double v = datagraph->id2seq[dst];

		if (datagraph->AdjList[u].size() == 1 || datagraph->AdjList[v].size() == 1)
		{
			cout << "ERROR: removing edges between illegal vertices!" << endl;
			return;
		}

		int flag1 = 0;
		for (double i = 0; i < datagraph->AdjList[u].size(); i++)
		{
			if (datagraph->AdjList[u][i] == v)
			{
				flag1 = 1;
				datagraph->AdjList[u].erase(datagraph->AdjList[u].begin() + i);
				break;
			}
		}

		int flag2 = 0;
		for (double i = 0; i < datagraph->AdjList[v].size(); i++)
		{
			if (datagraph->AdjList[v][i] == u)
			{
				flag2 = 1;
				datagraph->AdjList[v].erase(datagraph->AdjList[v].begin() + i);
				break;
			}
		}

		if (flag1 == 1 && flag2 == 1)
		{
			vector<double> collapsers;
			Subgraph* s;

			if (partition->coreness[u] < partition->coreness[v])
			{
				// if (src == 109 && dst == 4103) cout << u << "$" << v << endl;
				collapsers = partition->self_node[u]->collectCollapsers_outE(u);
				s = partition->self_node[u];

				Subgraph* s_ = partition->self_node[v];

				int u_seq = s_->vertIndex[u];
				//s_->vertIndex.erase(u);
				for (double i = 0; i < s_->adj2[u_seq].size(); i++)
				{
					int neigh = s_->adj2[u_seq][i];
					for (int j = 0; j < s_->adj1[neigh].size(); j++)
					{
						if (s_->adj1[neigh][j] == u_seq)
						{
							s_->adj1[neigh][j] = s_->adj1[neigh][s_->adj1[neigh].size() - 1];
							s_->adj1[neigh].pop_back();

							break;
						}
					}

					if (s_->adj.size() == 0)
						continue;

					for (int j = 0; j < s_->adj[neigh].size(); j++)
					{
						if (s_->adj[neigh][j] == u_seq)
						{
							s_->adj[neigh][j] = s_->adj[neigh][s_->adj[neigh].size() - 1];
							s_->adj[neigh].pop_back();

							break;
						}
					}
				}
				s_->adj2[u_seq].clear();
				if (s_->adj.size() != 0) s_->adj[u_seq].size();

				partition->Followers[u].erase(s_);
				// removing (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "-<: " << collapsers.size() << endl;
			}
			else if (partition->coreness[u] == partition->coreness[v])
			{
				collapsers = partition->self_node[u]->collectCollapsers_inE(u, v);
				s = partition->self_node[u];
				// cout << "-=: " << collapsers.size() << endl;
			}
			else // partition->coreness[u] > partition->coreness[v]
			{
				collapsers = partition->self_node[v]->collectCollapsers_outE(v);
				s = partition->self_node[v];

				Subgraph* s_ = partition->self_node[u];

				int v_seq = s_->vertIndex[v];
				//s_->vertIndex.erase(v);
				for (double i = 0; i < s_->adj2[v_seq].size(); i++)
				{
					int neigh = s_->adj2[v_seq][i];
					for (int j = 0; j < s_->adj1[neigh].size(); j++)
					{
						if (s_->adj1[neigh][j] == v_seq)
						{
							s_->adj1[neigh][j] = s_->adj1[neigh][s_->adj1[neigh].size() - 1];
							s_->adj1[neigh].pop_back();

							break;
						}
					}

					if (s_->adj.size() == 0)
						continue;

					for (int j = 0; j < s_->adj[neigh].size(); j++)
					{
						if (s_->adj[neigh][j] == v_seq)
						{
							s_->adj[neigh][j] = s_->adj[neigh][s_->adj[neigh].size() - 1];
							s_->adj[neigh].pop_back();

							break;
						}
					}
				}
				s_->adj2[v_seq].clear();
				if (s_->adj.size() != 0) s_->adj[v_seq].size();

				partition->Followers[v].erase(s_);
				// removing (u, v) has no effect to s_ in terms of collapsers
				// And s will be replaced anyway
				// cout << "->: " << collapsers.size() << endl;
			}

			for (int i = 0; i < collapsers.size(); i++)
			{
				double w = collapsers[i];
				partition->coreness[w]--;
				partition->shell_layer[w].first--;
				partition->shell_layer[w].second = INITIAL;
			}

			partition->new_nodes.clear();
			unordered_set<Subgraph*> deleted_nodes;
			deleted_nodes.insert(s);

			for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
			{
				double w = itt->first;
				int seqw = itt->second;

				if (s->shell[seqw] == s->k)
				{
					int flag = 0;
					for (int i = 0; i < partition->new_nodes.size(); i++)
					{
						if (partition->self_node[w] == partition->new_nodes[i])
						{
							flag = 1;
							break;
						}
					}

					if (flag == 1) continue;

					Subgraph* shell_node = new Subgraph();
					partition->new_nodes.push_back(shell_node);

					partition->Followers[w][shell_node] = INITIAL;
					partition->Collapsers[w][shell_node] = INITIAL;
					partition->self_node[w] = shell_node;

					shell_node->k = partition->coreness[w];
					shell_node->shell.push_back(partition->coreness[w]);
					shell_node->layer.push_back(INITIAL);
					vector<int> Neigh;
					shell_node->adj.push_back(Neigh);
					shell_node->highKSupport.push_back(0);
					shell_node->vertIndex[w] = shell_node->highKSupport.size() - 1;
					partition->shellConnect(w, shell_node, &deleted_nodes);
				}
			}

			for (unordered_set<Subgraph*>::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); it++)
			{
				Subgraph* ss = *it;
				partition->CC_Partitions.erase(ss);
				for (unordered_map<double, int>::iterator itt = ss->vertIndex.begin(); itt != ss->vertIndex.end(); itt++)
				{
					double w = itt->first;
					partition->Followers[w].erase(ss);
				}

				for (unordered_map<double, vector<int> >::iterator itt = ss->lowKNeighbors.begin(); itt != ss->lowKNeighbors.end(); itt++)
				{
					double w = itt->first;
					partition->Collapsers[w].erase(ss);
				}
			}

			// #pragma omp parallel for schedule(dynamic)
			#pragma omp parallel for schedule(dynamic, num_thr)
			for (int i = 0; i < partition->new_nodes.size(); i++)
			{
				Subgraph* ss = partition->new_nodes[i];
				ss->layerDecomposition();
				partition->CC_Partitions.insert(ss);
				// #pragma omp parallel for num_threads(num_threads)
				for (unordered_map<double, int>::iterator it = ss->vertIndex.begin(); it != ss->vertIndex.end(); it++)
				{
					double w = it->first;
					int seqw = it->second;

					int followers = ss->findFollowers_H(w);
					partition->Followers[w][ss] = followers;

					if (ss->shell[seqw] == ss->k)
					{
						int collapsers = ss->findCollapsers_in(w);
						partition->Collapsers[w][ss] = collapsers;
					}
				}
				// #pragma omp parallel for num_threads(num_threads)
				for (unordered_map<double, vector<int> >::iterator it = ss->lowKNeighbors.begin(); it != ss->lowKNeighbors.end(); it++)
				{
					double w = it->first;

					int collapsers = ss->findCollapsers_out(w);
					partition->Collapsers[w][ss] = collapsers;
				}
			}
		}
		else
		{
			cout << "ERROR: removing nonexistent edges!" << endl;
		}
	}
	else
	{
		cout << "ERROR: removing edges between illegal vertices!" << endl;
	}
}