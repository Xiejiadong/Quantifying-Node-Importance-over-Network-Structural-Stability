#include "DataGraph.h"

DataGraph::DataGraph(string fileName,int &N,int &M)
{
	vector<vector<double> >().swap(AdjList);
	unordered_map<double, double>().swap(id2seq);
	unordered_map<double, double>().swap(seq2id);

	fileName = "DataGraph/" + fileName;
	ifstream in(fileName.c_str());

	if (!in)
	{
		cout << "Fail to read " << fileName << "." << endl;
		return;
	}

	N = 0, M = 0;

	string line;
	
	while (getline(in, line))
	{
		string src_s = line.substr(0, line.find("\t"));
		string dst_s = line.substr(line.find("\t") + 1, line.find("\n") - line.find("\t") - 1);
		double src = stof(src_s);
		double dst = stof(dst_s);
		if (src != dst)
		{
			// cout<<src<<" "<<dst<<endl;
			unordered_map<double, double>::iterator it;
			int seq1, seq2;

			it = id2seq.find(src);
			if (it == id2seq.end())
			{
				seq1 = AdjList.size();
				id2seq[src] = seq1;
				seq2id[seq1] = src;
				vector<double> newVertex;
				AdjList.push_back(newVertex);
			}
			else
			{
				seq1 = it->second;
			}

			it = id2seq.find(dst);
			if (it == id2seq.end())
			{
				seq2 = AdjList.size();
				id2seq[dst] = seq2;
				seq2id[seq2] = dst;
				vector<double> newVertex;
				AdjList.push_back(newVertex);
			}
			else
			{
				seq2 = it->second;
			}

			int flag;

			flag = 0;
			for (double i = 0; i < AdjList[seq1].size(); i++)
			{
				if (AdjList[seq1][i] == seq2)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0) AdjList[seq1].push_back(seq2);

			flag = 0;
			for (double i = 0; i < AdjList[seq2].size(); i++)
			{
				if (AdjList[seq2][i] == seq1)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0) AdjList[seq2].push_back(seq1), M++;
			// cout<<M<<endl;
		}

	}
	N = id2seq.size();
	cout<< "vertices = " << N << " ; edges = " << M << endl;
}
