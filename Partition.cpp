#include "Partition.h"

bool compareDegree(const double& a, const double& b)
{
	return Partition::coreness[a] < Partition::coreness[b];
}

vector<int> Partition::coreness;

Partition::Partition(DataGraph* datagraph)
{
	this->datagraph = datagraph;
	max_degree = 0;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
	{
		allocated.push_back(false);
		self_node.push_back(NULL);
		unordered_map<Subgraph*, int> f, c;
		Followers.push_back(f);
		Collapsers.push_back(c);
		if (datagraph->AdjList[i].size() > max_degree)
			max_degree = datagraph->AdjList[i].size();
	}
}

void Partition::SL_Decomposition()
{
	coreness.clear(); //d
	shell_tag.clear(); //b
	total_order.clear(); //D
	order_pointer.clear(); //p

	shell_layer.clear();

	for (double i = 0; i < datagraph->AdjList.size(); i++)
	{
		total_order.push_back(i);
		int degree = datagraph->AdjList[i].size();
		coreness.push_back(degree);
		order_pointer.push_back(INITIAL);
		shell_layer.push_back(make_pair(INITIAL, INITIAL));
	}

	for (unordered_set<double>::iterator it = anchor_verts.begin(); it != anchor_verts.end(); it++)
	{
		int inf_degree = max_degree + datagraph->AdjList[*it].size();
		coreness[*it] = inf_degree;
	}

	sort(total_order.begin(), total_order.end(), compareDegree);
	int bi = 0;
	for (double i = 0; i < total_order.size(); i++)
	{
		double id = total_order[i];
		int degree = coreness[id];
		order_pointer[id] = i;

		while (degree >= bi)
		{
			shell_tag.push_back(i);
			bi++;
		}
	}
	
	for (double i = 0; i < total_order.size(); i++)
	{
		double v = total_order[i];

		if (shell_layer[v].first < coreness[v])
		{
			shell_layer[v].first = coreness[v];
			shell_layer[v].second = 1;
		}
		else if (shell_layer[v].first == coreness[v])
		{
			shell_layer[v].second = shell_layer[v].second + 1;
		}
		else;
		
		for (int j = 0; j < datagraph->AdjList[v].size(); j++)
		{
			double u = datagraph->AdjList[v][j];
			if (coreness[u] > coreness[v])
			{
				int du = coreness[u]; double pu = order_pointer[u];
				double pw = shell_tag[du]; double w = total_order[pw];

				if (u != w)
				{
					total_order[pu] = w; total_order[pw] = u;
					order_pointer[u] = pw; order_pointer[w] = pu;
				}

				shell_tag[du]++; coreness[u]--;

				shell_layer[u].first = shell_layer[v].first;
				shell_layer[u].second = shell_layer[v].second;
			}
		}
	}
}

void Partition::P_Decomposition()
{
	SL_Decomposition();
	int num=0,Max=0,Min=1000000000;
	double average=0;
	for (double i = 0; i < total_order.size(); i++)
	{
		double u = total_order[i];

		if (allocated[u] == false && anchor_verts.find(u) == anchor_verts.end())
		{
			num++;
			Subgraph* shell_node = new Subgraph();
			CC_Partitions.insert(shell_node);

			allocated[u] = true;
			Followers[u][shell_node] = INITIAL;
			Collapsers[u][shell_node] = INITIAL;
			self_node[u] = shell_node;

			shell_node->k = coreness[u];
			shell_node->shell.push_back(shell_layer[u].first);
			shell_node->layer.push_back(shell_layer[u].second); 
			vector<int> lowNeigh, highNeigh;
			shell_node->adj1.push_back(lowNeigh);
			shell_node->adj2.push_back(highNeigh);
			shell_node->highKSupport.push_back(0);
			shell_node->vertIndex[u] = shell_node->highKSupport.size() - 1;
			shell_node->num_nodes=1;
			shellConnect(u, shell_node);
		}
	}
	
	for (unordered_set<Subgraph*>::iterator it = CC_Partitions.begin(); it != CC_Partitions.end(); it++)
		new_nodes.push_back(*it);

}

void Partition::P_Maintenance(double anchor)
{
	vector<double> followers;
	for (unordered_map<Subgraph*, int>::iterator it = Followers[anchor].begin(); it != Followers[anchor].end(); it++)
	{
		vector<double> ff = it->first->collectFollowers(anchor);
		for (int i = 0; i < ff.size(); i++)
			followers.push_back(ff[i]);
	}
	
	anchor_verts.insert(anchor);
	coreness[anchor] = max_degree;
	shell_layer[anchor].first = max_degree;
	shell_layer[anchor].second = INITIAL;

	for (int i = 0; i < followers.size(); i++)
	{
		double u = followers[i];
		coreness[u]++;
		shell_layer[u].first++;
		shell_layer[u].second = INITIAL;
	}
	
	unordered_set<Subgraph*> deleted_nodes;
	new_nodes.clear();

	for (unordered_map<Subgraph*, int>::iterator it = Followers[anchor].begin(); it != Followers[anchor].end(); it++)
	{
		Subgraph* s = it->first;
		deleted_nodes.insert(s);
		for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
		{
			double v = itt->first;
			int seqv = itt->second;

			if (v == anchor) continue;

			if (s->shell[seqv] == s->k)
			{
				int flag = 0;
				for (int i = 0; i < new_nodes.size(); i++)
				{
					if (self_node[v] == new_nodes[i])
					{
						flag = 1;
						break;
					}
				}
				
				if (flag == 1) continue;

				Subgraph* shell_node = new Subgraph();
				new_nodes.push_back(shell_node);

				Followers[v][shell_node] = INITIAL;
				self_node[v] = shell_node;

				shell_node->k = coreness[v];
				shell_node->shell.push_back(coreness[v]);
				shell_node->layer.push_back(INITIAL);
				vector<int> Neigh;
				shell_node->adj.push_back(Neigh);
				shell_node->highKSupport.push_back(0);
				shell_node->vertIndex[v] = shell_node->highKSupport.size() - 1;
				shellConnect(v, shell_node, &deleted_nodes);
			}	
		}
	}

	for (unordered_set<Subgraph*>::iterator it = deleted_nodes.begin(); it != deleted_nodes.end(); it++)
	{
		Subgraph* s = *it;
		CC_Partitions.erase(s);
		for (unordered_map<double, int>::iterator itt = s->vertIndex.begin(); itt != s->vertIndex.end(); itt++)
		{
			double v = itt->first;
			Followers[v].erase(s);
		}		
	}

	for (int i = 0; i < new_nodes.size(); i++)
	{
		new_nodes[i]->layerDecomposition();
		CC_Partitions.insert(new_nodes[i]);
	}
}

void Partition::shellConnect(double u, Subgraph* shell_node)
{
	int sequ = shell_node->vertIndex[u];

	for (int i = 0; i < datagraph->AdjList[u].size(); i++)
	{
		double v = datagraph->AdjList[u][i];

		if (coreness[v] < shell_node->k)
		{
			if (shell_node->vertIndex.find(v) == shell_node->vertIndex.end())
			{
				shell_node->shell.push_back(shell_layer[v].first);
				shell_node->layer.push_back(shell_layer[v].second);
				vector<int> lowNeigh, highNeigh;
				shell_node->adj1.push_back(lowNeigh);
				shell_node->adj2.push_back(highNeigh);
				shell_node->highKSupport.push_back(0);
				shell_node->vertIndex[v] = shell_node->highKSupport.size() - 1;

				Followers[v][shell_node] = INITIAL;
			}

			int seqv = shell_node->vertIndex[v];
			shell_node->adj1[sequ].push_back(seqv);
			shell_node->adj2[seqv].push_back(sequ);
		}
		else if (coreness[v] == shell_node->k)
		{
			if (allocated[v] == true)
			{
				int seqv = shell_node->vertIndex[v];
				if (shell_layer[v].second < shell_layer[u].second)
					shell_node->adj1[sequ].push_back(seqv);
				else if (shell_layer[v].second == shell_layer[u].second)
					shell_node->adj1[sequ].push_back(seqv);
				else //shell_layer[v].second > shell_layer[u].second
					shell_node->adj2[sequ].push_back(seqv);
			}
			else //allocated[v] == false
			{				
				shell_node->shell.push_back(shell_layer[v].first);
				shell_node->layer.push_back(shell_layer[v].second);
				vector<int> lowNeigh, highNeigh;
				shell_node->adj1.push_back(lowNeigh);
				shell_node->adj2.push_back(highNeigh);
				shell_node->highKSupport.push_back(0);
				shell_node->vertIndex[v] = shell_node->highKSupport.size() - 1;

				int seqv = shell_node->vertIndex[v];
				if (shell_layer[v].second < shell_layer[u].second)
					shell_node->adj1[sequ].push_back(seqv);
				else if (shell_layer[v].second == shell_layer[u].second)
					shell_node->adj1[sequ].push_back(seqv);
				else //shell_layer[v].second > shell_layer[u].second
					shell_node->adj2[sequ].push_back(seqv);
				
				allocated[v] = true;
				self_node[v] = shell_node;
				Followers[v][shell_node] = INITIAL;
				Collapsers[v][shell_node] = INITIAL;
				shell_node->num_nodes++;
				shellConnect(v, shell_node);
			}
		}
		else if (coreness[v] > shell_node->k)
		{
			shell_node->highKSupport[sequ]++;

			if (shell_node->lowKNeighbors.find(v) == shell_node->lowKNeighbors.end())
			{
				vector<int> neighs;
				shell_node->lowKNeighbors[v] = neighs;

				Collapsers[v][shell_node] = INITIAL;
			}

			shell_node->lowKNeighbors[v].push_back(sequ);
		}
		else;
	}
}

void Partition::shellConnect(double u, Subgraph* shell_node, unordered_set<Subgraph*>* deleted_nodes)
{
	int sequ = shell_node->vertIndex[u];

	for (int i = 0; i < datagraph->AdjList[u].size(); i++)
	{
		double v = datagraph->AdjList[u][i];

		if (coreness[v] < shell_node->k)
		{
			if (shell_node->vertIndex.find(v) == shell_node->vertIndex.end())
			{
				shell_node->shell.push_back(coreness[v]);
				shell_node->layer.push_back(INITIAL);
				vector<int> Neigh;
				shell_node->adj.push_back(Neigh);
				shell_node->highKSupport.push_back(0);
				shell_node->vertIndex[v] = shell_node->highKSupport.size() - 1;

				Followers[v][shell_node] = INITIAL;
			}

			int seqv = shell_node->vertIndex[v];
			shell_node->adj[sequ].push_back(seqv);
			shell_node->adj[seqv].push_back(sequ);
		}
		else if (coreness[v] == shell_node->k)
		{
			if (self_node[v] != shell_node)
			{
				deleted_nodes->insert(self_node[v]);

				shell_node->shell.push_back(coreness[v]);
				shell_node->layer.push_back(INITIAL);
				vector<int> Neigh;
				shell_node->adj.push_back(Neigh);
				shell_node->highKSupport.push_back(0);
				shell_node->vertIndex[v] = shell_node->highKSupport.size() - 1;

				int seqv = shell_node->vertIndex[v];
				shell_node->adj[sequ].push_back(seqv);

				self_node[v] = shell_node;
				Followers[v][shell_node] = INITIAL;
				Collapsers[v][shell_node] = INITIAL;

				shellConnect(v, shell_node, deleted_nodes);
			}
			else
			{
				int seqv = shell_node->vertIndex[v];
				shell_node->adj[sequ].push_back(seqv);
			}
		}
		else if (coreness[v] > shell_node->k)
		{
			shell_node->highKSupport[sequ]++;

			if (shell_node->lowKNeighbors.find(v) == shell_node->lowKNeighbors.end())
			{
				vector<int> neighs;
				shell_node->lowKNeighbors[v] = neighs;

				Collapsers[v][shell_node] = INITIAL;
			}

			shell_node->lowKNeighbors[v].push_back(sequ);
		}
		else;
	}
}
