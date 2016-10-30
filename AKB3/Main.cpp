#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "Sequence.cpp"
using namespace std;

void loadFile(string);
void printSeq();
void printGraph();
void makeVertexes();
void dotPlot();

Sequnece seq;
vector<SequenceFragment*> vertexes;
int window, minScore, **graph;

int main() {
	cout << "podaj nazwe pliku do otwarcia bez rozszerzenia" << endl;
	string fileName;
	cin >> fileName;
	cout << endl;

	loadFile(fileName);

	cout << "podaj rozmiar okna(4-7)" << endl;
	cin >> window;
	cout << endl << "podaj minimalna tolerowana ocene nukleotydu" << endl;
	cin >> minScore;

	makeVertexes();
	dotPlot();
	printGraph();
	cout << "zakonczono prace programu, wcisnij dowolny klawisz by wyjsc...." << endl;
	getchar();
	getchar();

	return 0;
}

void printSeq() {
	cout << seq.seqID<<endl << seq.seq << endl << seq.scoreString << endl;
}
void printGraph() {
	for (int i = 0; i < vertexes.size(); i++) {
		for (int j = 0; j < vertexes.size(); j++) {
			cout << graph[i][j] << " ";
		}
		cout << endl;
	}
}
void loadFile(string fileName) {
	fstream seqFile, scoreFile;
	vector<string> seqIds;
	string line;

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
	cout << "prosze podac nr sekwencji" << endl;
	int seqNum;
	cin >> seqNum;
	cout << endl;

	seqFile.clear();
	seqFile.seekg(0);
	while (!seqFile.eof()) {
		getline(seqFile, line, '\n');
		if (line == seqIds[seqNum]) {
				seq.seqID = line;
				getline(seqFile, line, '\n');
				while (line[0] != '>' && !seqFile.eof()) {
					line.erase(remove(line.begin(), line.end(), '\n'), line.end());
					seq.seq += line;
					getline(seqFile, line, '\n');
				}
				break;
			}
		}
	
	while (!scoreFile.eof()){
		getline(scoreFile, line, '\n');
		if (line == seq.seqID) {
				getline(scoreFile, line, '\n');
				while (line[0] != '>' && !scoreFile.eof()) {
					line.erase(remove(line.begin(), line.end(), '\n'), line.end());
					seq.scoreString += line;
					getline(scoreFile, line, '\n');
				}
				break;
			}
		}
	printSeq();
}

void makeVertexes() {
	int size = seq.seq.size() / window;
	seq.makeScoreTable();
	for (int i = 0; i <= seq.seq.size() - window; i++) {
		vector<char> seqSubString;
		vector<int> scoreSubString;
		for (int j = i; j < i+window; j++) {
			seqSubString.push_back(seq.seq[j]);
			scoreSubString.push_back(seq.score[j]);
		}
		vertexes.push_back(new SequenceFragment(seq.seqID, seqSubString, scoreSubString, i));
		seqSubString.clear();
		scoreSubString.clear();
	}
	graph = new int*[vertexes.size()];
	for (int i = 0; i < vertexes.size(); i++) {
		graph[i] = new int[vertexes.size()];
		for (int j = 0; j < vertexes.size(); j++) {
			graph[i][j] = 0;
		}
	}
}

void dotPlot() {
	int same = 0;
	for (int i = 0; i < vertexes.size() - 1; i++) {
		for (int j = i + 1; j < vertexes.size(); j++) {
			for (int k = 0; k < window; k++) {
				if (vertexes[i]->seq[k] == vertexes[j]->seq[k]) {
					same++;
				}
			}
			if (same == window) {
				graph[i][i + 1] = 1;
				graph[i + 1][i] = 1;
			}
			same = 0;
		}
	}
}