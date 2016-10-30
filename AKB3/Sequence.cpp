#include <vector>
#include <string>
#include <sstream>

using namespace std;

struct SequenceFragment {
	vector<char> seq;
	vector<int>score;
	string seqId;
	int begin;

	SequenceFragment(string seqId, vector<char> seq, vector<int> score,int begin) {
		this->seq = seq;
		this->score = score;
		this->seqId = seqId;
		this->begin = begin;
	}
};

struct Sequnece {
	string seqID, seq, scoreString;
	vector<int> score;

	void makeScoreTable() {
		stringstream stream(scoreString);
		while (stream) {
			int buffer;
			stream >> buffer;
			score.push_back(buffer);
		}
	}
};