#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <windows.h> 
#include <map>
#include "Structs.cpp"
using namespace std;


void loadFile(string);

void buildGraph();
void makeConnections();

vector<int> findCommonNodes(vector<int>, int);
void bronKerbosh(vector<int>, vector<int>, vector<int>);
Clinque findClique();

Motive findMotive(Clinque clq);
void maximalizeMotive(Motive motive);

map<string,Sequnece> seqences;
vector<SequenceFragment*> vertexes;
vector<vector<int>> connections;
vector<Clinque> cliques;
const vector<char> Nucleotides = { 'A','C', 'T', 'G','_' };//nukleotydy i _ jako delecja
Motive motive;
int window, minScore,motiveCounter=1, deletionLevel;
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

	buildGraph();
	cout << "utworzylo wierzcho³chiki, zaczynam tworzyc polaczenia" << endl;
	makeConnections();
	cout << "polaczenia utworzne, zaczynam szukac motywow" << endl;
		
	findMotive(findClique());
		
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
				seqences[seq.seqID]=seq;
			}
		}
	}
	seqFile.close();
	scoreFile.close();
}

//--------TWORZENIE GRAFU
void buildGraph() {
	for (auto seq : seqences) {
		int size = seq.second.seq.size() / window;
		seq.second.makeScoreTable();
		for (int i = 0; i <= seq.second.seq.size() - window; i++) {
			vector<char> seqSubString;
			vector<int> scoreSubString;
			for (int j = i; j < i + window; j++) {
				seqSubString.push_back(seq.second.seq[j]);
				scoreSubString.push_back(seq.second.score[j]);
			}
			vertexes.push_back(new SequenceFragment(seq.second.seqID, seqSubString, scoreSubString, i));
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
	deletionLevel = ((window / 2) - 1);
	int same = 0, deletion = 1;;

	for (int i = 0; i < vertexes.size() - 1; i++) {
		for (int j = i+1 ; j < vertexes.size(); j++) {
			for (int k = 0; k < window; k++) {
				if (vertexes[i]->seq[k] == vertexes[j]->seq[k]) {
					same++;
				}
				//A  C  T     A  T	i - sprawdzane okno
				//30 20 20    20 21
				//A  C  T  G  A  T	j - okno sprawdzajace
				//30 20 20 19 20 21
				else if (vertexes[j]->score[k] <= minScore) {
					deletion++;
					int shift = k+1;
					while (deletion <= deletionLevel && shift < window) {
						if (vertexes[i]->seq[shift] == vertexes[j]->seq[k]) {
							same++;
							shift++;
							k++;
						}
						else {
							deletion++;
							shift++;
						}
					}
				}
				//A  C  T  G  A  T	i
				//30 20 20 19 20 21
				//A  C  T     A  T	j
				//30 20 20    20 21
				else if (vertexes[i]->score[k] <= minScore) {
					deletion++;
					int shift = k + 1;
					while (deletion <= deletionLevel && shift < window) {
						if (vertexes[i]->seq[k] == vertexes[j]->seq[shift]) {
							same++;
							shift++;
							k++;
						}
						else {
							deletion++;
							shift++;
						}
					}
				}

			}
			if (same + deletion == window && deletion <= deletionLevel) {
				connections[i].push_back(j);
				connections[j].push_back(i);
			}
			same = 0;
			deletion = 0;
		}
	}
}
//--------SZUKANIE KLIK
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
void bronKerbosh(vector<int> potentialClique, vector<int> nodes, vector<int> skipedNodes) {
	if (nodes.size() == 0 && skipedNodes.size() == 0) {
		cliques.push_back(Clinque(potentialClique));
		//cout << "znalezione klike nr " << motiveCounter++ <<endl;
		return;
	}

	for (int i = 0; i < nodes.size(); i++) {
		vector<int> newPotentialClique = potentialClique;
		newPotentialClique.push_back(nodes[i]);
		vector<int> newNodes = findCommonNodes(nodes, nodes[i]);
		vector<int> newSkippedNodes = findCommonNodes(skipedNodes, nodes[i]);

		bronKerbosh(newPotentialClique, newNodes, newSkippedNodes);

		newPotentialClique.pop_back();
		skipedNodes.push_back(nodes[i]);
		nodes.erase(remove(nodes.begin(), nodes.end(), nodes[i]), nodes.end());


	}
}
Clinque findClique() {
	int maxDeg = 0, maxDegID;
	for (int i = 0; i < connections.size();i++) {
		if (connections[i].size() > maxDeg) {
			maxDeg = connections[i].size();
			maxDegID = i;
		}
	}
	vector<int> potentialClique, skipedNodes;
	bronKerbosh(potentialClique, connections[maxDegID], skipedNodes);
	int maxCliqSize = 0, maxCliqID;
	for (int i = 0; i < cliques.size();i++) {
		if (cliques[i].list.size() > maxCliqSize) {
			maxCliqSize = cliques[i].list.size();
			maxCliqID = i;
		}
	}
	return cliques[maxCliqID];
}
//--------SZUKANIE MOTYWU
Motive findMotive(Clinque clq) {
	Motive motive;

	for (int i = 0; i < window; i++) {
		map<char, int> consensusColumn;
		consensusColumn['A'] = 0;
		consensusColumn['C'] = 0;
		consensusColumn['T'] = 0;
		consensusColumn['G'] = 0;
		consensusColumn['_'] = 0;
		for (int id : clq.list) {
			consensusColumn[vertexes[id]->seq[i]]++;
		}
		bool nucWasAdded = false;
		for (auto nuc : Nucleotides) {
			if (consensusColumn[nuc] >= window / 2) {
				motive.consensus.push_back(nuc);
				nucWasAdded = true;
			}
		}
		if (!nucWasAdded)
			motive.consensus.push_back('_');
	}
	return motive;
}

void maximalizeMotive(Motive motive){
	int misMatchNum = 0;
	while(misMatchNum < deletionLevel){
		Clinque clqLeft, clqRight;
		for(int id : motive.fragments.begin){//lewa strona
			int ext = id - window < 0 ? -1 : id - window;
			if (ext > 0) {
				clqLeft.list.push_back(ext);
			}
		}
		for(auto frag : motive.fragments){//prawa strona
			int ext = frag.end + window <= seqences[frag.seqId].seq.size() ? -1 : frag.end + window;
			if (ext > 0) {
				clqRight.list.push_back(ext);
			}
		}


	}
}
