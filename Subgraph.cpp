#include "Subgraph.h"

void Subgraph::layerDecomposition()
{
	vector<int> d, b, D, p;

	for (int i = 0; i < adj.size(); i++)
	{
		D.push_back(i);
		d.push_back(adj[i].size() + highKSupport[i]);
		p.push_back(INITIAL);
	}
	
	for (int i = 1; i < D.size(); i++)
	{
		int u = D[i];

		int j = 0;
		while (j < i && d[u] > d[D[j]])
			j++;

		for (int m = i; m > j; m--)
			D[m] = D[m - 1];

		D[j] = u;
	}
	
	int bi = 0;
	for (int i = 0; i < D.size(); i++)
	{
		int id = D[i];
		int degree = d[id];
		p[id] = i;

		while (degree >= bi)
		{
			b.push_back(i);
			bi++;
		}
	}
	
	for (int i = 0; i < D.size(); i++)
	{
		int v = D[i];

		if (shell[v] == k)
			layer[v] += 1;

		for (int j = 0; j < adj[v].size(); j++)
		{
			double u = adj[v][j];
			if (d[u] > d[v])
			{
				int du = d[u]; int pu = p[u];
				int pw = b[du]; int w = D[pw];

				if (u != w)
				{
					D[pu] = w; D[pw] = u;
					p[u] = pw; p[w] = pu;
				}

				b[du]++; d[u]--;

				if (shell[v] == k)
					layer[u] = layer[v];
			}
		}
	}
	
	for (int i = 0; i < adj.size(); i++)
	{
		int u = i;
		vector<int> empty;
		adj1.push_back(empty);
		adj2.push_back(empty);
		for (int j = 0; j < adj[u].size(); j++)
		{
			int v = adj[u][j];

			if (shell[v] < shell[u])
				adj1[u].push_back(v);
			else if (shell[u] == shell[v] && layer[v] <= layer[u])
				adj1[u].push_back(v);
			else if (shell[u] == shell[v] && layer[v] > layer[u])
				adj2[u].push_back(v);
			else //shell[v] > shell[u]
				adj2[u].push_back(v);
		}
	}
}

vector<double> Subgraph::findFollowers_D(double x_)
{
	vector<int> d, b, D, p;

	for (int i = 0; i < adj1.size(); i++)
	{
		D.push_back(i);
		d.push_back(adj1[i].size() + adj2[i].size() + highKSupport[i]);
		p.push_back(INITIAL);
	}

	int x = vertIndex[x_];
	d[x] = vertIndex.size() + k + 1;

	for (int i = 0; i < D.size(); i++)
	{
		if (i == 0)
			continue;

		int u = D[i];

		int j = 0;
		while (d[u] > d[D[j]] && j < i)
			j++;

		for (int m = i; m > j; m--)
			D[m] = D[m - 1];

		D[j] = u;
	}

	int bi = 0;
	for (int i = 0; i < D.size(); i++)
	{
		int id = D[i];
		int degree = d[id];
		p[id] = i;

		while (degree >= bi)
		{
			b.push_back(i);
			bi++;
		}
	}

	for (int i = 0; i < D.size(); i++)
	{
		int v = D[i];

		for (int j = 0; j < adj1[v].size(); j++)
		{
			double u = adj1[v][j];
			if (d[u] > d[v])
			{
				int du = d[u]; int pu = p[u];
				int pw = b[du]; int w = D[pw];

				if (u != w)
				{
					D[pu] = w; D[pw] = u;
					p[u] = pw; p[w] = pu;
				}

				b[du]++; d[u]--;
			}
		}

		for (int j = 0; j < adj2[v].size(); j++)
		{
			double u = adj2[v][j];
			if (d[u] > d[v])
			{
				int du = d[u]; int pu = p[u];
				int pw = b[du]; int w = D[pw];

				if (u != w)
				{
					D[pu] = w; D[pw] = u;
					p[u] = pw; p[w] = pu;
				}

				b[du]++; d[u]--;
			}
		}
	}

	vector<double> followers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (d[seqv] == k + 1 && v != x_)
			followers.push_back(v);
	}

	return followers;
}

int Subgraph::findFollowers_H(double x_)
{
	vector<int> degreeBound, vertexStatus;
	degreeBound.clear();
	vertexStatus.clear();
	for (int i = 0; i < vertIndex.size(); i++)
	{
		degreeBound.push_back(INITIAL);
		vertexStatus.push_back(UNEXPLORED);
	}
	
	priority_queue<Node> H;
	unordered_set<int> H_set;
	unordered_set<int> Survivals;

	int x = vertIndex[x_];
	degreeBound[x] = vertIndex.size() + k + 1;
	vertexStatus[x] = SURVIVED;
	H.push(Node(0, x));
	H_set.insert(x);

	while (H.size() != 0)
	{
		Node u = H.top();
		H.pop();

		if (u.id != x && vertexStatus[u.id] != DISCARDED)
		{
			int d_plus_u = 0;
			for (int i = 0; i < adj1[u.id].size(); i++)
			{
				int v = adj1[u.id][i];
				if (vertexStatus[v] != DISCARDED && H_set.find(v) != H_set.end())
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (vertexStatus[v] != DISCARDED)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u.id];

			degreeBound[u.id] = d_plus_u;
		}

		if (degreeBound[u.id] >= k + 1)
		{
			vertexStatus[u.id] = SURVIVED;
			Survivals.insert(u.id);
			
			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (H_set.find(v) == H_set.end())
				{
					H.push(Node(layer[v], v));
					H_set.insert(v);
				}
			}
		}
		else
		{
			vertexStatus[u.id] = DISCARDED;
			Survivals.erase(u.id);
			shrink(x, u.id, &Survivals,degreeBound,vertexStatus);
		}
	}

	int followers = Survivals.size();

	return --followers;
}

vector<double> Subgraph::collectFollowers(double x_)
{
	degreeBound.clear();
	vertexStatus.clear();
	for (int i = 0; i < vertIndex.size(); i++)
	{
		degreeBound.push_back(INITIAL);
		vertexStatus.push_back(UNEXPLORED);
	}

	priority_queue<Node> H;
	unordered_set<int> H_set;
	unordered_set<int> Survivals;

	int x = vertIndex[x_];
	degreeBound[x] = vertIndex.size() + k + 1;
	vertexStatus[x] = SURVIVED;
	H.push(Node(0, x));
	H_set.insert(x);

	while (H.size() != 0)
	{
		Node u = H.top();
		H.pop();

		if (u.id != x && vertexStatus[u.id] != DISCARDED)
		{
			int d_plus_u = 0;
			for (int i = 0; i < adj1[u.id].size(); i++)
			{
				int v = adj1[u.id][i];
				if (vertexStatus[v] != DISCARDED && H_set.find(v) != H_set.end())
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (vertexStatus[v] != DISCARDED)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u.id];

			degreeBound[u.id] = d_plus_u;
		}

		if (degreeBound[u.id] >= k + 1)
		{
			vertexStatus[u.id] = SURVIVED;
			Survivals.insert(u.id);

			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (H_set.find(v) == H_set.end())
				{
					H.push(Node(layer[v], v));
					H_set.insert(v);
				}
			}
		}
		else
		{
			vertexStatus[u.id] = DISCARDED;
			Survivals.erase(u.id);
			shrink(x, u.id, &Survivals,degreeBound,vertexStatus);
		}
	}

	vector<double> followers;

	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seq = it->second;

		if (vertexStatus[seq] == SURVIVED && v != x_)
			followers.push_back(v);
	}

	return followers;
}

void Subgraph::shrink(int x, int u, unordered_set<int>* Survivals,vector<int>& degreeBound,vector<int>& vertexStatus)
{
	vector<int> shrink_neighbours;

	for (int i = 0; i < adj1[u].size(); i++)
	{
		int v = adj1[u][i];
		if (v != x && vertexStatus[v] == SURVIVED)
		{
			degreeBound[v]--;
			if (degreeBound[v] < k + 1)
			{
				vertexStatus[v] = DISCARDED;
				Survivals->erase(v);
				shrink_neighbours.push_back(v);
			}
		}
	}
	for (int i = 0; i < adj2[u].size(); i++)
	{
		int v = adj2[u][i];
		if (v != x && vertexStatus[v] == SURVIVED)
		{
			degreeBound[v]--;
			if (degreeBound[v] < k + 1)
			{
				vertexStatus[v] = DISCARDED;
				Survivals->erase(v);
				shrink_neighbours.push_back(v);
			}
		}
	}

	for (int i = 0; i < shrink_neighbours.size(); i++)
	{
		int neigh = shrink_neighbours[i];
		shrink(x, neigh, Survivals,degreeBound,vertexStatus);
	}
}

int Subgraph::findCollapsers_in(double x_)
{
	int x = vertIndex[x_];

	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	explored_queue.push(x);
	explored_set.insert(x);

	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		if (u != x)
		{
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && collapsed[v] == false)
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (collapsed[v] == false)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u];
		}

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}

	int collapsers = 0;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true)
			collapsers += 1;
	}
	collapsers -= 1; // x itself
	
	return collapsers;
}

int Subgraph::findCollapsers_out(double x_)
{
	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	for (int i = 0; i < lowKNeighbors[x_].size(); i++)
	{
		int u = lowKNeighbors[x_][i];
		explored_queue.push(u);
		explored_set.insert(u);

		highKSupport[u]--;
	}

	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		for (int i = 0; i < adj1[u].size(); i++)
		{
			int v = adj1[u][i];
			if (shell[v] == k && collapsed[v] == false)
				d_plus_u++;
		}
		for (int i = 0; i < adj2[u].size(); i++)
		{
			int v = adj2[u][i];
			if (collapsed[v] == false)
				d_plus_u++;
		}
		d_plus_u += highKSupport[u];

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}

	for (int i = 0; i < lowKNeighbors[x_].size(); i++)
	{
		highKSupport[lowKNeighbors[x_][i]]++;
	}
	
	int collapsers = 0;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true)
			collapsers += 1;
	}
	
	return collapsers;
}

vector<double> Subgraph::collectCollapsers_in(double x_)
{
	int x = vertIndex[x_];

	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	explored_queue.push(x);
	explored_set.insert(x);

	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		if (u != x)
		{
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && collapsed[v] == false)
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (collapsed[v] == false)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u];
		}

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}

	vector<double> collapsers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true && v != x_)
			collapsers.push_back(v);
	}

	return collapsers;
}

vector<double> Subgraph::collectCollapsers_out(double x_)
{
	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	for (int i = 0; i < lowKNeighbors[x_].size(); i++)
	{
		int u = lowKNeighbors[x_][i];
		explored_queue.push(u);
		explored_set.insert(u);

		highKSupport[u]--;
	}

	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		for (int i = 0; i < adj1[u].size(); i++)
		{
			int v = adj1[u][i];
			if (shell[v] == k && collapsed[v] == false)
				d_plus_u++;
		}
		for (int i = 0; i < adj2[u].size(); i++)
		{
			int v = adj2[u][i];
			if (collapsed[v] == false)
				d_plus_u++;
		}
		d_plus_u += highKSupport[u];

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}

	for (int i = 0; i < lowKNeighbors[x_].size(); i++)
	{
		highKSupport[lowKNeighbors[x_][i]]++;
	}

	vector<double> collapsers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true && v != x_)
			collapsers.push_back(v);
	}

	return collapsers;
}

vector<double> Subgraph::collectFollowers_inE(double x_, double y_)
{
	int x = vertIndex[x_];
	int y = vertIndex[y_];

	if (layer[x] == layer[y])
	{
		int d_plus_x = adj2[x].size() + highKSupport[x] + 1;
		int d_plus_y = adj2[y].size() + highKSupport[y] + 1;

		if (d_plus_x >= k + 1 && d_plus_y >= k + 1)
		{
			adj1[x].push_back(y);
			adj1[y].push_back(x);

			degreeBound.clear();
			vertexStatus.clear();
			for (int i = 0; i < vertIndex.size(); i++)
			{
				degreeBound.push_back(INITIAL);
				vertexStatus.push_back(UNEXPLORED);
			}

			priority_queue<Node> H;
			unordered_set<int> H_set;

			degreeBound[x] = d_plus_x;
			vertexStatus[x] = SURVIVED;
			H.push(Node(layer[x], x));
			H_set.insert(x);

			degreeBound[y] = d_plus_y;
			vertexStatus[y] = SURVIVED;
			H.push(Node(layer[y], y));
			H_set.insert(y);

			while (H.size() != 0)
			{
				Node u = H.top();
				H.pop();

				if (u.id != x && u.id != y && vertexStatus[u.id] != DISCARDED)
				{
					int d_plus_u = 0;
					for (int i = 0; i < adj1[u.id].size(); i++)
					{
						int v = adj1[u.id][i];
						if (vertexStatus[v] != DISCARDED && H_set.find(v) != H_set.end())
							d_plus_u++;
					}
					for (int i = 0; i < adj2[u.id].size(); i++)
					{
						int v = adj2[u.id][i];
						if (vertexStatus[v] != DISCARDED)
							d_plus_u++;
					}
					d_plus_u += highKSupport[u.id];

					degreeBound[u.id] = d_plus_u;
				}

				if (degreeBound[u.id] >= k + 1)
				{
					vertexStatus[u.id] = SURVIVED;

					for (int i = 0; i < adj2[u.id].size(); i++)
					{
						int v = adj2[u.id][i];
						if (H_set.find(v) == H_set.end())
						{
							H.push(Node(layer[v], v));
							H_set.insert(v);
						}
					}
				}
				else
				{
					vertexStatus[u.id] = DISCARDED;;
					shrink_E(u.id);
				}
			}

			adj1[x].pop_back();
			adj1[y].pop_back();
		}
		else
		{
			return vector<double>();
		}
	}
	else
	{
		if (layer[y] < layer[x])
		{
			int z = x;
			x = y;
			y = z;
		}

		int d_plus_x = adj2[x].size() + highKSupport[x] + 1;

		if (d_plus_x >= k + 1)
		{
			adj2[x].push_back(y);
			adj1[y].push_back(x);

			degreeBound.clear();
			vertexStatus.clear();
			for (int i = 0; i < vertIndex.size(); i++)
			{
				degreeBound.push_back(INITIAL);
				vertexStatus.push_back(UNEXPLORED);
			}

			priority_queue<Node> H;
			unordered_set<int> H_set;

			degreeBound[x] = d_plus_x;
			vertexStatus[x] = SURVIVED;
			H.push(Node(layer[x], x));
			H_set.insert(x);

			while (H.size() != 0)
			{
				Node u = H.top();
				H.pop();

				if (u.id != x && vertexStatus[u.id] != DISCARDED)
				{
					int d_plus_u = 0;
					for (int i = 0; i < adj1[u.id].size(); i++)
					{
						int v = adj1[u.id][i];
						if (vertexStatus[v] != DISCARDED && H_set.find(v) != H_set.end())
							d_plus_u++;
					}
					for (int i = 0; i < adj2[u.id].size(); i++)
					{
						int v = adj2[u.id][i];
						if (vertexStatus[v] != DISCARDED)
							d_plus_u++;
					}
					d_plus_u += highKSupport[u.id];

					degreeBound[u.id] = d_plus_u;
				}

				if (degreeBound[u.id] >= k + 1)
				{
					vertexStatus[u.id] = SURVIVED;

					for (int i = 0; i < adj2[u.id].size(); i++)
					{
						int v = adj2[u.id][i];
						if (H_set.find(v) == H_set.end())
						{
							H.push(Node(layer[v], v));
							H_set.insert(v);
						}
					}
				}
				else
				{
					vertexStatus[u.id] = DISCARDED;
					shrink_E(u.id);
				}
			}

			adj2[x].pop_back();
			adj1[y].pop_back();
		}
		else
		{
			return vector<double>();
		}
	}

	vector<double> followers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (vertexStatus[seqv] == SURVIVED)
			followers.push_back(v);
	}

	return followers;
}

vector<double> Subgraph::collectFollowers_outE(double x_)
{
	int x = vertIndex[x_];
	int d_plus_x = adj2[x].size() + highKSupport[x] + 1;

	if (d_plus_x < k + 1)
		return vector<double>();

	degreeBound.clear();
	vertexStatus.clear();
	for (int i = 0; i < vertIndex.size(); i++)
	{
		degreeBound.push_back(INITIAL);
		vertexStatus.push_back(UNEXPLORED);
	}

	priority_queue<Node> H;
	unordered_set<int> H_set;

	degreeBound[x] = d_plus_x;
	vertexStatus[x] = SURVIVED;
	H.push(Node(layer[x], x));
	H_set.insert(x);
	
	while (H.size() != 0)
	{
		Node u = H.top();
		H.pop();

		if (u.id != x && vertexStatus[u.id] != DISCARDED)
		{
			int d_plus_u = 0;
			for (int i = 0; i < adj1[u.id].size(); i++)
			{
				int v = adj1[u.id][i];
				if (vertexStatus[v] != DISCARDED && H_set.find(v) != H_set.end())
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (vertexStatus[v] != DISCARDED)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u.id];

			degreeBound[u.id] = d_plus_u;
		}

		if (degreeBound[u.id] >= k + 1)
		{
			vertexStatus[u.id] = SURVIVED;

			for (int i = 0; i < adj2[u.id].size(); i++)
			{
				int v = adj2[u.id][i];
				if (H_set.find(v) == H_set.end())
				{
					H.push(Node(layer[v], v));
					H_set.insert(v);
				}
			}
		}
		else
		{
			vertexStatus[u.id] = DISCARDED;
			shrink_E(u.id);
		}
	}
	
	vector<double> followers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;
		
		if (vertexStatus[seqv] == SURVIVED)
			followers.push_back(v);
	}

	return followers;
}

vector<double> Subgraph::collectCollapsers_inE(double x_, double y_)
{
	int x = vertIndex[x_];
	int d_plus_x = adj2[x].size() + highKSupport[x] - 1;
	for (int i = 0; i < adj1[x].size(); i++)
		if (shell[adj1[x][i]] == k)
			d_plus_x++;

	int y = vertIndex[y_];
	int d_plus_y = adj2[y].size() + highKSupport[y] - 1;
	for (int i = 0; i < adj1[y].size(); i++)
		if (shell[adj1[y][i]] == k)
			d_plus_y++;

	if (d_plus_x >= k && d_plus_y >= k)
		return vector<double>();

	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	if (d_plus_x >= k && d_plus_y < k)
	{
		int z = x;
		x = y;
		y = z;
	}

	explored_queue.push(x);
	explored_set.insert(x);

	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		if (u != x)
		{
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && collapsed[v] == false)
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (collapsed[v] == false)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u];
		}

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}

	vector<double> collapsers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true)
			collapsers.push_back(v);
	}

	return collapsers;
}

vector<double> Subgraph::collectCollapsers_outE(double x_)
{
	int x = vertIndex[x_];
	int d_plus_x = adj2[x].size() + highKSupport[x] - 1;
	for (int i = 0; i < adj1[x].size(); i++)
		if (shell[adj1[x][i]] == k)
			d_plus_x++;

	if (d_plus_x >= k)
		return vector<double>();

	vector<bool> collapsed;
	collapsed.resize(vertIndex.size(), false);

	queue<int> explored_queue;
	unordered_set<int> explored_set;

	explored_queue.push(x);
	explored_set.insert(x);
	
	while (explored_queue.size() != 0)
	{
		int u = explored_queue.front();
		explored_queue.pop();
		explored_set.erase(u);

		int d_plus_u = 0;

		if (u != x)
		{
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && collapsed[v] == false)
					d_plus_u++;
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (collapsed[v] == false)
					d_plus_u++;
			}
			d_plus_u += highKSupport[u];
		}

		if (d_plus_u < k)
		{
			collapsed[u] = true;
			for (int i = 0; i < adj1[u].size(); i++)
			{
				int v = adj1[u][i];
				if (shell[v] == k && explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
			for (int i = 0; i < adj2[u].size(); i++)
			{
				int v = adj2[u][i];
				if (explored_set.find(v) == explored_set.end() && collapsed[v] == false)
				{
					explored_queue.push(v);
					explored_set.insert(v);
				}
			}
		}
	}
	
	vector<double> collapsers;
	for (unordered_map<double, int>::iterator it = vertIndex.begin(); it != vertIndex.end(); it++)
	{
		double v = it->first;
		int seqv = it->second;

		if (collapsed[seqv] == true)
			collapsers.push_back(v);
	}

	return collapsers;
}

void Subgraph::shrink_E(int u)
{
	vector<int> shrink_neighbours;

	for (int i = 0; i < adj1[u].size(); i++)
	{
		int v = adj1[u][i];
		if (vertexStatus[v] == SURVIVED)
		{
			degreeBound[v]--;
			if (degreeBound[v] < k + 1)
			{
				vertexStatus[v] = DISCARDED;
				shrink_neighbours.push_back(v);
			}
		}
	}
	for (int i = 0; i < adj2[u].size(); i++)
	{
		int v = adj2[u][i];
		if (vertexStatus[v] == SURVIVED)
		{
			degreeBound[v]--;
			if (degreeBound[v] < k + 1)
			{
				vertexStatus[v] = DISCARDED;
				shrink_neighbours.push_back(v);
			}
		}
	}

	for (int i = 0; i < shrink_neighbours.size(); i++)
	{
		int neigh = shrink_neighbours[i];
		shrink_E(neigh);
	}
}