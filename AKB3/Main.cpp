#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include "Sequence.cpp"
using namespace std;

void loadFile(string);
void print(vector<vector<int>>);
void printSeq();
void printGraph();
void makeVertexes();
vector<int> findPivots(vector<int>, vector<int>);
vector<int> findCommonNodes(vector<int>, int);
void findCliques(vector<int>, vector<int>, vector<int>);
void dotPlot();

void findCliques(vector<int> potentialClique, vector<int> nodes, vector<int> skiped_nodes);

vector<Sequnece> seqences;
vector<SequenceFragment*> vertexes;
vector<vector<int>> connections,cliques;
int window, minScore;
random_device rd;
mt19937 mt(rd());

int main() {
	cout << "podaj nazwe pliku do otwarcia bez rozszerzenia" << endl;
	string fileName;
	cin >> fileName;
	cout << endl;

	loadFile(fileName);
	//for debuging
	printSeq();
	cout << "wczytano pliki" << endl;
	cout << "podaj rozmiar okna(4-7)" << endl;
	cin >> window;
	cout << endl << "podaj minimalna tolerowana ocene nukleotydu" << endl;
	cin >> minScore;

	makeVertexes();
	cout << "utworzylo wierzcho³chiki, zaczynam tworzyc polaczenia" << endl;
	dotPlot();
	cout << "polaczenia utworzne, zaczynam szukac max motywu" << endl;
	//for debuging
	printGraph();
	
	vector<int> emptyPotentialClique, emptySkipedNodes, vertexList;
	for (int i = 0; i < connections.size(); i++)
		vertexList.push_back(i);
	findCliques(emptyPotentialClique, vertexList, emptySkipedNodes);
	
	//for debuging
	print(cliques);
	cout << "zakonczono prace programu, wcisnij dowolny klawisz by wyjsc...." << endl;
	getchar();
	getchar();

	return 0;
}

void printSeq() {
	cout << endl;
	cout << "sekwencje i ich identyfikatory" << endl;
	for (auto seq : seqences) {
		cout << seq.seqID << endl;
		cout << seq.seq << endl;
		cout << seq.scoreString << endl;
	}
	cout << endl;
}
void printGraph() {
	for (int i = 0; i < connections.size(); i++) {
		if (connections[i].size() != 0) {
			cout << i << ": ";
			for (int v : connections[i]) {
				cout << v << " ";
			}
			cout << endl;
		}
	}
}
void loadFile(string fileName) {
	fstream seqFile, scoreFile;
	vector<string> seqIds;
	string line, scoreline;

	seqFile.open(fileName+".fasta", ios::in);
	scoreFile.open(fileName+".qual", ios::in);

	while (!seqFile.eof()) {
		getline(seqFile, line, '\n');
		if (line[0] == '>') {
			seqIds.push_back(line);
		}
	}
	cout << "lista identyfikatorow sekwencji w pliku: " << endl;
	for (int i = 0; i < seqIds.size();i++) {
		cout << i<<". " << seqIds[i] << endl;
	}
	
	//rewinding file to begining
	seqFile.clear();
	seqFile.seekg(0);

	getline(seqFile, line, '\n');
	while (!seqFile.eof()) {
		
		if (line[0] == '>') {
			Sequnece seq;
			seq.seqID = line;
			getline(seqFile, line, '\n');
			while (line[0] != '>' && !seqFile.eof()) {
				line.erase(remove(line.begin(), line.end(), '\n'), line.end());
				seq.seq += line;
				getline(seqFile, line, '\n');
			}
			while (!scoreFile.eof()) {
				if(scoreline[0] != '>')
					getline(scoreFile, scoreline, '\n');
				if (scoreline == seq.seqID) {
					getline(scoreFile, scoreline, '\n');
					while (scoreline[0] != '>' && !scoreFile.eof()) {
						scoreline.erase(remove(scoreline.begin(), scoreline.end(), '\n'), scoreline.end());
						seq.scoreString += scoreline;
						getline(scoreFile, scoreline, '\n');
					}
					break;
				}
			}
			if (seq.seqID.size() > 3) {
				seqences.push_back(seq);
			}
		}
	}
	seqFile.close();
	scoreFile.close();
}
void print(vector<vector<int>> input)
{
	cout << endl<< "potencialne kliki:" << endl;
	for (auto i : input) {
		for (auto j : i) {
			for (auto nuc : vertexes[j]->seq) 
				cout << nuc; 
			cout << "->"<<j<<" ";
		}
		cout << "\n";
	}
}
void printMotive() {

}
void makeVertexes() {
	for (auto seq : seqences) {
		int size = seq.seq.size() / window;
		seq.makeScoreTable();
		for (int i = 0; i <= seq.seq.size() - window; i++) {
			vector<char> seqSubString;
			vector<int> scoreSubString;
			for (int j = i; j < i + window; j++) {
				seqSubString.push_back(seq.seq[j]);
				scoreSubString.push_back(seq.score[j]);
			}
			vertexes.push_back(new SequenceFragment(seq.seqID, seqSubString, scoreSubString, i));
			seqSubString.clear();
			scoreSubString.clear();
		}
	}
	for (int i = 0; i < vertexes.size(); i++) {
		vector<int> temp;
		connections.push_back(temp);
	}
}
void dotPlot() {
	int deletionLevel = ((window / 2) - 1);
	int same = 0, deletion = 0;;
	for (int i = 0; i < vertexes.size() - 1; i++) {
		for (int j = i + 1; j < vertexes.size(); j++) {
			for (int k = 0; k < window; k++) {
				if (vertexes[i]->seq[k] == vertexes[j]->seq[k]) {
					same++;
				}
				else if (vertexes[j]->score[k] <= minScore) {
					deletion++;
				}
			}
			if (same+deletion == window && deletion <= deletionLevel) {
				connections[i].push_back(j);
				connections[j].push_back(i);
			}
			same = 0;
			deletion = 0;
		}
	}
}
vector<int> findPivots(vector<int> nodes, vector<int> skipedNodes) {
	vector<int> pivots(nodes);
	if(!skipedNodes.empty())
		for (auto i : skipedNodes)
			pivots.push_back(i);
	return pivots;
}
vector<int> findCommonNodes(vector<int> nodes, int vertex) {
	vector<int> common;
	for (int n : nodes) {
		for (int i = 0; i < connections[vertex].size(); i++) {
			if (connections[vertex][i] == n) {
				common.push_back(n);
			}
		}
	}
	return common;
}
vector<int> pivotsNeibers(vector<int> nodes, int pivot) {
	vector<int> toReturn;
	int checked = 0;
	for (auto node : nodes) {
		for (int i = 0; i < connections[pivot].size(); i++) {
			if (connections[pivot][i] != node) {
				checked++;
			}
		}
		if (checked == connections[pivot].size())
			toReturn.push_back(node);
		checked = 0;
	}
	return toReturn;
}
int randomizePivot(int size) {	
	uniform_int_distribution<int> dist(0, size-1);
	int test = dist(mt);
	return test;
}
void findCliques(vector<int> potentialClique, vector<int> nodes, vector<int> skipedNodes) {
	if (nodes.size() == 0 && skipedNodes.size() == 0) {
		cliques.push_back(potentialClique);
		cout << "znaleziono klike, szukam dalej..." << endl;
		return;
	}

	vector<int> pivots(findPivots(nodes, skipedNodes));
	vector<int> nodesWithoutPivotsNeibers(pivotsNeibers(nodes, pivots[randomizePivot(pivots.size())]));

	for (int i=0; i < nodesWithoutPivotsNeibers.size();i++) {
		vector<int> newPotentialClique = potentialClique;
		newPotentialClique.push_back(nodesWithoutPivotsNeibers[i]);
		vector<int> newNodes = findCommonNodes(nodesWithoutPivotsNeibers, nodesWithoutPivotsNeibers[i]) ;
		vector<int> newSkippedNodes = findCommonNodes(skipedNodes, nodesWithoutPivotsNeibers[i]);

		findCliques(newPotentialClique, newNodes, newSkippedNodes);

		newPotentialClique.pop_back();
		skipedNodes.push_back(nodesWithoutPivotsNeibers[i]);
		nodesWithoutPivotsNeibers.erase(remove(nodesWithoutPivotsNeibers.begin(), nodesWithoutPivotsNeibers.end(), nodesWithoutPivotsNeibers[i]), nodesWithoutPivotsNeibers.end());
		
		
	}
}