#include "Combination.h"

bool compareDegree1(const int& a, const int& b)
{
	return Combination::anch_degrees[a] < Combination::anch_degrees[b];
}

bool comparePotential(const int& a, const int& b)
{
	return Combination::vert_potentials[a] > Combination::vert_potentials[b];
}

vector<int> Combination::anch_degrees;
vector<double> Combination::vert_potentials;

Combination::Combination(DataGraph* datagraph)
{
	this->datagraph = datagraph;
}

vector<int> Combination::coreDecomposition(vector<int>& anchors)
{
	vector<int>d, b, D, p;

	for (int i = 0; i < datagraph->AdjList.size(); i++)
	{
		D.push_back(i);
		int degree = datagraph->AdjList[i].size();
		d.push_back(degree);
		p.push_back(INITIAL);
	}

	anch_degrees = d;

	for (int i = 0; i < anchors.size(); i++)
	{
		/* set the degrees of anchor vertices as infinity */
		int inf_degree = datagraph->AdjList.size() * 2;
		anch_degrees[anchors[i]] = inf_degree;
		d[anchors[i]] = inf_degree;
	}

	/* Initialization of d, b, D and p */
	sort(D.begin(), D.end(), compareDegree1);
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

	/* Bottom-Up k-core decomposition */
	for (int i = 0; i < D.size(); i++)
	{
		int v = D[i];
		for (int j = 0; j < datagraph->AdjList[v].size(); j++)
		{
			int u = datagraph->AdjList[v][j];
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

	return d;
}

void Combination::sampleRand(int samp_num, double br)
{
	vector<vector<int> >().swap(comb_verts);
	vector<vector<int> >().swap(v_comb_seqs);
	vector<vector<double> >().swap(v_comb_shares);

	v_comb_seqs.resize(datagraph->AdjList.size(), vector<int>());
	v_comb_shares.resize(datagraph->AdjList.size(), vector<double>());

	int rand_range = 1 / br;
	for (int i = 0; i < samp_num; i++)
	{
		vector<int> anchor_set;
		for (int j = 0; j < datagraph->AdjList.size(); j++)
		{
			random_device rd;
			uniform_int_distribution<int> dist(0, rand_range - 1);
			int seed = dist(rd);
			if (seed == 0)
			{
				int anc_v = j;
				anchor_set.push_back(anc_v);
				v_comb_seqs[anc_v].push_back(comb_verts.size());
				v_comb_shares[anc_v].push_back(0.0);
			}	
		}
		comb_verts.push_back(anchor_set);
	}
}

void Combination::updStatistics(double r1, double r2, double threshold, double k_, double b_)
{
	vector<double> comb_gain_share;
	for (int i = 0; i < comb_verts.size(); i++)
	{
		vector<int>& anchors = comb_verts[i];
		vector<int> coreness_anc = coreDecomposition(anchors);

		double coreness_gain = 0;
		for (int j = 0; j < coreness_init.size(); j++)
			if (coreness_anc[j] < datagraph->AdjList.size())
				coreness_gain += ((double)coreness_anc[j] - (double)coreness_init[j]);

		cout << "random: " << coreness_gain << endl;
		coreness_gain /= anchors.size();
		comb_gain_share.push_back(coreness_gain);
	}
	
	for (int i = 0; i < v_comb_seqs.size(); i++)
		for (int j = 0; j < v_comb_seqs[i].size(); j++)
			v_comb_shares[i][j] = comb_gain_share[v_comb_seqs[i][j]];

	double avg_samp_num = 0;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
	{
		int v = i;

		double avg_v = v_samp_num[v] * v_share_avg[v];
		for (int j = 0; j < v_comb_shares[v].size(); j++)
			avg_v += v_comb_shares[v][j];
		if (v_samp_num[v] + v_comb_shares[v].size() != 0)
			avg_v = avg_v / (v_samp_num[v] + (double)v_comb_shares[v].size());
		else
			avg_v = 0;
		
		double var_v = v_samp_num[v] * (v_share_var[v] + (v_share_avg[v] - avg_v) * (v_share_avg[v] - avg_v));
		for (int j = 0; j < v_comb_shares[v].size(); j++)
			var_v += (v_comb_shares[v][j] - avg_v) * (v_comb_shares[v][j] - avg_v);
		if (v_samp_num[v] + v_comb_shares[v].size() != 0)
			var_v = var_v / (v_samp_num[v] + (double)v_comb_shares[v].size());
		else
			var_v = 0;

		v_samp_num[v] += v_comb_shares[v].size();
		v_share_avg[v] = avg_v;
		v_share_var[v] = var_v;

		avg_samp_num += v_samp_num[v];
	}
	avg_samp_num /= datagraph->AdjList.size();

	double samp_num_thresh = avg_samp_num * threshold;
	vector<double>().swap(vert_potentials);
	vert_potentials.resize(datagraph->AdjList.size(), 0.0);
	int effec_num = 0;
	double effec_avg = 0;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
	{
		int v = i;

		if (v_samp_num[v] < samp_num_thresh)
			continue;

		double v_potential = k_ * (r1 * v_share_avg[v] - r2 * v_share_var[v]) + b_;
		vert_potentials[v] = v_potential;

		effec_num += 1;
		effec_avg += v_potential;
	}

	effec_avg /= effec_num;
	double effec_std = 0;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
	{
		int v = i;

		if (v_samp_num[v] < samp_num_thresh)
			continue;

		double v_gap_avg = vert_potentials[v] - effec_avg;
		if (v_gap_avg >= 0)
			effec_std += v_gap_avg;
		else
			effec_std -= v_gap_avg;
	}
	effec_std /= effec_num;

	cout << "Effective avg: " << effec_avg << "; " << "Effective std: " << effec_std << endl;
}

int Combination::topbTest(double br)
{
	vector<int> potential_order;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
		potential_order.push_back(i);
	
	sort(potential_order.begin(), potential_order.end(), comparePotential);

	vector<int> anchors;
	int b = datagraph->AdjList.size() * br;
	for (int i = 0; i < b; i++)
		anchors.push_back(potential_order[i]);

	vector<int> coreness_anc = coreDecomposition(anchors);

	int coreness_gain = 0;
	for (int i = 0; i < coreness_init.size(); i++)
		if (coreness_anc[i] < datagraph->AdjList.size())
			coreness_gain += (coreness_anc[i] - coreness_init[i]);

	return coreness_gain;
}

void Combination::trackSampling(string dataset, double br, int track_freq, double threshold, 
	double r1, double r2, double k_, double b_)
{
	ofstream out("tracks/track_" + dataset);
	out << "Budget ratio: " << br << endl;
	out << "Tracking frequency: " << track_freq << endl;
	out << "Effective threshold ratio: " << threshold << endl;
	out << "r1(avg): " << r1 << endl;
	out << "r2(var): " << r2 << endl;
	out << "k(scal): " << k_ << endl;
	out << "b(scal): " << b_ << endl;
	out << endl;

	vector<int> none_anchors;
	coreness_init = coreDecomposition(none_anchors);

	v_samp_num.resize(datagraph->AdjList.size(), 0);
	v_share_avg.resize(datagraph->AdjList.size(), 0.0);
	v_share_var.resize(datagraph->AdjList.size(), 0.0);
	int track_round = 0;

	double time1 = omp_get_wtime();
	while (1)
	{
		sampleRand(track_freq, br);
		updStatistics(r1, r2, threshold, k_, b_);
		int coreness_gain = topbTest(br);
		double time2 = omp_get_wtime();

		track_round += track_freq;
		double time_cost = time2 - time1;
		out << "Round " << track_round << ": coreness gain = " << coreness_gain << "; running time = "<< time_cost << endl;
		cout << "Round " << track_round << ": coreness gain = " << coreness_gain << "; running time = " << time_cost << endl;

		ofstream out_train("train_data/potential_" + to_string(track_round) + "_" + dataset);
		for (int i = 0; i < datagraph->AdjList.size(); i++)
			out_train << i << "\t" << vert_potentials[i] << endl;
	}
}

void Combination::assignAnchors()
{
	/*cout << "potential: ";
	string potential;
	getline(cin, potential);

	potential = "potentials/" + potential;
	ifstream in(potential.c_str());

	string line;
	while (getline(in, line))
	{
		string ptl = line.substr(line.find("\t") + 1, line.find("\n") - line.find("\t") - 1);
		vert_potentials.push_back(stof(ptl));
	}

	vector<int> potential_order;
	for (int i = 0; i < datagraph->AdjList.size(); i++)
		potential_order.push_back(i);

	sort(potential_order.begin(), potential_order.end(), comparePotential);

	vector<int> anchors;
	int b = datagraph->AdjList.size() * 0.05;
	for (int i = 0; i < b; i++)
		anchors.push_back(potential_order[i]);*/
	

	cout << "gac anchors: ";
	string gac_anchors;
	getline(cin, gac_anchors);

	gac_anchors = "gac_anchors/" + gac_anchors;
	ifstream in(gac_anchors.c_str());

	vector<int> anchors_full;
	string line;
	while (getline(in, line))
	{
		if (line.substr(0, 1) == ":") break;
		string anchor_s = line.substr(line.find("-") + 2, line.find("'") - line.find("-") - 2);
		anchors_full.push_back(stoi(anchor_s));
	}

	vector<int> anchors;
	for (int i = 0; i < anchors_full.size() * 0.33; i++)
		anchors.push_back(anchors_full[i]);

	//////////////////////////////////////////////////////////////
	vector<int> none_anchors;
	coreness_init = coreDecomposition(none_anchors);

	vector<int> coreness_anc = coreDecomposition(anchors);

	int coreness_gain = 0;
	for (int i = 0; i < coreness_init.size(); i++)
		if (coreness_anc[i] < datagraph->AdjList.size())
			coreness_gain += (coreness_anc[i] - coreness_init[i]);

	cout << "coreness gain: " << coreness_gain << endl;
}
