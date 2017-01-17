#include <vector>
#include <string>
#include <sstream>

using namespace std;

struct SequenceVertex
{
	string seq;
	string seqId, seqString;
	int begin;
	
	SequenceVertex(string seqId, string seq, int begin)
	{
		this->seq = seq;
		this->seqId = seqId;
		this->begin = begin;
	}

	SequenceVertex()
	{
	}
};

struct Sequnece
{
	string seqID, seq, scoreString;
	vector<int> score;

	void makeScoreTable()
	{
		stringstream stream(scoreString);
		while(stream)
		{
			int buffer;
			stream >> buffer;
			score.push_back(buffer);
		}
	}
};

struct Connection
{
	int vertexId, score;
	Connection(){};

	Connection(int vertex_id, int score)
		: vertexId(vertex_id),
		  score(score)
	{
	}
	bool operator ==(const Connection &v) const
	{
		if (this->vertexId == v.vertexId) return true;
		return false;
	}
};

struct SequenceFragment
{
	string seqId;
	string seq;
	int begining, end, id;
	SequenceFragment(string seq_id, int begining, int end, int id,string seq)
		: seqId(seq_id),
		  begining(begining),
		  end(end),
		  id(id),
		  seq(seq)
	{
	}
};

struct Motive
{
	string consensus;
	vector<SequenceFragment> fragments;

	Motive()
	{
	}
};
