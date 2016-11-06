#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <windows.h> 
#include "Sequence.cpp"
using namespace std;

void loadFile(string);
void printPotentialCliques(vector<vector<int>>);
void printSeq();
void printGraph();
void printMotive();
void printSequenceWihtMotiv();
void makeCleanGraph();
vector<int> findPivots(vector<int>, vector<int>);
vector<int> findCommonNodes(vector<int>, int);
void findCliques(vector<int>, vector<int>, vector<int>);
void makeConnections();

void findCliques(vector<int> potentialClique, vector<int> nodes, vector<int> skiped_nodes);

vector<Sequnece> seqences;
vector<SequenceFragment*> vertexes;
vector<vector<int>> connections,cliques;
vector<int> maxClinque;
int window, minScore,motiveCounter=1;
HANDLE  hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

int main() {
	
	cout << "podaj nazwe pliku do otwarcia bez rozszerzenia" << endl;
	string fileName;
	cin >> fileName;
	cout << endl;

	loadFile(fileName);
	//for debuging
	//printSeq();
	cout << "wczytano pliki" << endl;
	cout << "podaj rozmiar okna(4-7)" << endl;
	cin >> window;
	cout << endl << "podaj minimalna tolerowana ocene nukleotydu" << endl;
	cin >> minScore;

	makeCleanGraph();
	cout << "utworzylo wierzcho³chiki, zaczynam tworzyc polaczenia" << endl;
	makeConnections();
	cout << "polaczenia utworzne, zaczynam szukac motywow" << endl;
	//for debuging
	//printGraph();
	
	vector<int> emptyPotentialClique, emptySkipedNodes, vertexList;
	for (int i = 0; i < connections.size(); i++)
		vertexList.push_back(i);
	findCliques(emptyPotentialClique, vertexList, emptySkipedNodes);
	cout << "znaleziono ~" << cliques.size() << " motywów, szukam najdluzszego" << endl;
	//for debuging
	//print(cliques);

	for (auto qc : cliques)
		if (maxClinque.size() < qc.size())
			maxClinque = qc;
	printMotive();
	printSequenceWihtMotiv();
	cout << "zakonczono prace programu, wcisnij dowolny klawisz by wyjsc...." << endl;
	getchar();
	getchar();

	return 0;
}
void loadFile(string fileName) {
	fstream seqFile, scoreFile;
	vector<string> seqIds;
	string line, scoreline;

	seqFile.open(fileName + ".fasta", ios::in);
	scoreFile.open(fileName + ".qual", ios::in);

	while (!seqFile.eof()) {
		getline(seqFile, line, '\n');
		if (line[0] == '>') {
			seqIds.push_back(line);
		}
	}
	cout << "lista identyfikatorow sekwencji w pliku: " << endl;
	for (int i = 0; i < seqIds.size(); i++) {
		cout << i << ". " << seqIds[i] << endl;
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
				if (scoreline[0] != '>')
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
//--------PRINTS
//for debuging
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
void printPotentialCliques(vector<vector<int>> input)
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
//works
void printMotive() {
	for (auto i : maxClinque) {
		for (auto ch : vertexes[i]->seq) {
			cout << ch;
		}
		cout << " ";
	}
	cout << endl;
}
void printSequenceWihtMotiv()
{
	vector<int> startes;
	SequenceFragment *temp;
	SetConsoleTextAttribute(hConsole, 15);

	for (auto seq : seqences) {
		cout << seq.seqID << endl;
		for (auto id: maxClinque) {
			if (vertexes[id]->seqId == seq.seqID) {
				startes.push_back(vertexes[id]->begin);
			}
		}
		for (int i = 0; i < seq.seq.size(); i++) {
			for (auto nn : startes) {
				if (i == nn) {
					SetConsoleTextAttribute(hConsole, 2);
					for (int j = 0; j < window; j++, i++) {
						cout << seq.seq[i];
					}
					SetConsoleTextAttribute(hConsole, 15);
				}
			}
			cout << seq.seq[i];
		}
		cout << endl<< endl;
		startes.clear();
	}
}
//--------TWORZENIE GRAFU
void makeCleanGraph() {
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
void makeConnections() {
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
//--------SZUKANIE MOTYWU
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
void findCliques(vector<int> potentialClique, vector<int> nodes, vector<int> skipedNodes) {
	if (nodes.size() == 0 && skipedNodes.size() == 0) {
		cliques.push_back(potentialClique);
		//cout << "znalezione klike nr " << motiveCounter++ <<endl;
		return;
	}
	
	for (int i=0; i < nodes.size();i++) {
		vector<int> newPotentialClique = potentialClique;
		newPotentialClique.push_back(nodes[i]);
		vector<int> newNodes = findCommonNodes(nodes, nodes[i]) ;
		vector<int> newSkippedNodes = findCommonNodes(skipedNodes, nodes[i]);

		findCliques(newPotentialClique, newNodes, newSkippedNodes);

		newPotentialClique.pop_back();
		skipedNodes.push_back(nodes[i]);
		nodes.erase(remove(nodes.begin(), nodes.end(), nodes[i]), nodes.end());
		
		
	}
}