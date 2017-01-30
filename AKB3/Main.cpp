#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <windows.h>
#include <map>
#include "Structs.cpp"
using namespace std;

//--------METODY
void loadFile(string);

void removeLowScoredNuc();
void buildGraph();
void makeConnections();

vector<Connection> findCommonNodes(vector<Connection>, int);
void bronKerbosh(vector<Connection>, vector<Connection>, vector<Connection>);
vector<Connection> findClique();

Motive findMotive(vector<Connection>);
Motive maximalizeMotive(Motive);

void print(SequenceFragment);

//--------ZMIENNE
map<string, Sequnece> sekwencje;//id sekwencji, sekwencje, punktacja
vector<SequenceVertex*> krawedzie;//sekwencje podzielone na okna
vector<vector<Connection>> polaczenia;//lista nastepnikow z ocena dopasowania
vector<vector<Connection>> kliki;
const vector<char> Nukleotydy = {'A','C', 'T', 'G','_'};//nukleotydy i _ jako delecja
Motive motyw;
int windowSize, minScore, motiveCounter = 1, deletionLevel;

HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

int main()
{
	cout << "podaj nazwe pliku do otwarcia bez rozszerzenia" << endl;
	string fileName;
	cin >> fileName;
	cout << endl;

	loadFile(fileName);
	cout << "wczytano pliki" << endl;
	cout << "podaj rozmiar okna(4-7)" << endl;
	cin >> windowSize;
	deletionLevel = windowSize * 0.4f;
	cout << endl << "podaj minimalna tolerowana ocene nukleotydu" << endl;
	cin >> minScore;

	removeLowScoredNuc();
	buildGraph();
	cout << "utworzylo wierzcho³chiki, zaczynam tworzyc polaczenia" << endl;
	makeConnections();
	cout << "polaczenia utworzne, zaczynam szukac motywu" << endl;
		
	motyw = findMotive(findClique());
	motyw = maximalizeMotive(motyw);
	for (auto seq : motyw.fragments)
		print(seq);
	cout << "zakonczono prace programu, wcisnij dowolny klawisz by wyjsc...." << endl;
	getchar();
	getchar();

	return 0;
}

void loadFile(string fileName)
{
	fstream seqFile, scoreFile;
	vector<string> seqIds;
	string line, scoreline;

	seqFile.open(fileName + ".fasta", ios::in);
	scoreFile.open(fileName + ".qual", ios::in);

	while (!seqFile.eof())
	{
		getline(seqFile, line, '\n');
		if (line[0] == '>')
		{
			seqIds.push_back(line);
		}
	}
	cout << "lista identyfikatorow sekwencji w pliku: " << endl;
	for (int i = 0; i < seqIds.size(); i++)
	{
		cout << i << ". " << seqIds[i] << endl;
	}

	//rewinding file to begining
	seqFile.clear();
	seqFile.seekg(0);

	getline(seqFile, line, '\n');
	while (!seqFile.eof())
	{
		if (line[0] == '>')
		{
			Sequnece seq;
			seq.seqID = line;
			getline(seqFile, line, '\n');
			while (line[0] != '>')
			{
				line.erase(remove(line.begin(), line.end(), '\n'), line.end());
				seq.seq += line;
				if (seqFile.eof()) break;
				getline(seqFile, line, '\n');
				
			}
			while (!scoreFile.eof())
			{
				if (scoreline[0] != '>')
					getline(scoreFile, scoreline, '\n');
				if (scoreline == seq.seqID)
				{
					getline(scoreFile, scoreline, '\n');
					while (scoreline[0] != '>' )
					{
						scoreline.erase(remove(scoreline.begin(), scoreline.end(), '\n'), scoreline.end());
						seq.scoreString += scoreline;
						if (scoreFile.eof()) break;
						getline(scoreFile, scoreline, '\n');

					}
					break;
				}
			}
			if (seq.seqID.size() > 3)
			{
				sekwencje[seq.seqID] = seq;
			}
		}
	}
	seqFile.close();
	scoreFile.close();
}

//--------TWORZENIE GRAFU
void removeLowScoredNuc()
{
	for(auto sekwencja : sekwencje)
	{
		sekwencja.second.makeScoreTable();
		for(int i=0;i<sekwencja.second.score.size()-1;i++)
		{
			if(sekwencja.second.score[i]<minScore)
			{
				sekwencja.second.score.erase(remove(sekwencja.second.score.begin(), sekwencja.second.score.end(), i), sekwencja.second.score.end());
				sekwencja.second.seq.erase(remove(sekwencja.second.seq.begin(), sekwencja.second.seq.end(),i), sekwencja.second.seq.end());
			}
		}
		sekwencje[sekwencja.second.seqID] = sekwencja.second;
	}
//		for(auto sekwencja : sekwencje)
//	{
//		sekwencja.second.makeScoreTable();
//		for(int i=0;i<sekwencja.second.score.size()-1;i++)
//		{
//			if(sekwencja.second.score[i]<minScore)
//			{
//				sekwencja.second.score.erase(sekwencja.second.score.begin()+(i-1));
//				sekwencja.second.seq.erase(sekwencja.second.seq.begin() + (i - 1));
//			}
//		}
//		sekwencje[sekwencja.second.seqID] = sekwencja.second;
//	}
}
void buildGraph()
{
	for (auto seq : sekwencje)
	{
		int size = seq.second.seq.size() - windowSize;
		for (int i = 0; i <= size; i++)
		{
			string seqSubString;
			
			for (int j = i; j < i + windowSize; j++)
			{
				seqSubString.push_back(seq.second.seq[j]);
				
			}
			krawedzie.push_back(new SequenceVertex(seq.second.seqID, seqSubString, i));
			seqSubString.clear();
			
		}
	}
	for (int i = 0; i < krawedzie.size(); i++)
	{
		vector<Connection> temp;
		polaczenia.push_back(temp);
	}
}
void makeConnections()
{
	int same = 0, deletion = 0;;
	
	for (int i = 0; i < krawedzie.size() - 1; i++)
	{
		for (int j = i+1; j < krawedzie.size(); j++)
		{
			if(krawedzie[i]->seqId == krawedzie[j]->seqId)
				continue;
			for (int k = 0; k < windowSize; k++)
			{
				//z pominieciem 
				//A  C  T  C  A  T
				//|	 |	|	  |	 |
				//A  C  T  G  A  T	
				if (krawedzie[i]->seq[k] == krawedzie[j]->seq[k])
				{
					same++;
				}
			}

			if (same >= windowSize - deletionLevel )
			{
				Connection connectI(i,same);
				Connection connectJ(j,same);
				polaczenia[i].push_back(connectJ);
				polaczenia[j].push_back(connectI);
			}
//			else
//			{
//				//przesuniecie w lewo
//				//A  C  T  A  T	 G
//				//|	 |	|	 \	\	 
//				//A  C  T  G  A  T
//				deletion = 0;
//				same = 0;
//				for (int k = 0; k < windowSize; k++) {
//					if (krawedzie[i]->seq[k] == krawedzie[j]->seq[k]) {
//						same++;
//					}
//					else {
//						deletion++;
//						int shift = k + 1;
//						while (deletion <= deletionLevel && shift < windowSize) {
//							if (krawedzie[i]->seq[shift] == krawedzie[j]->seq[k]) {
//								same++;
//								shift++;
//								k++;
//							}
//							else {
//								deletion++;
//								shift++;
//							}
//						}
//					}
//				}
//				if (same + deletion == windowSize && deletion <= deletionLevel)
//				{
//					Connection connectI(i, same);
//					Connection connectJ(j, same);
//					polaczenia[i].push_back(connectJ);
//					polaczenia[j].push_back(connectI);
//				}
//				else {
//					//przesuniecie w prawo
//					//A  C  G  T  G	 A
//					//|	 |	  /  /  /	 
//					//A  C  T  G  A  T
//					deletion = 0;
//					same = 0;
//					for (int k = 0; k < windowSize; k++) {
//						if (krawedzie[i]->seq[k] == krawedzie[j]->seq[k]) {
//							same++;
//						}
//						else {
//							deletion++;
//							int shift = k + 1;
//							while (deletion <= deletionLevel && shift < windowSize) {
//								if (krawedzie[j]->seq[shift] == krawedzie[i]->seq[k]) {
//									same++;
//									shift++;
//									k++;
//								}
//								else {
//									deletion++;
//									shift++;
//								}
//							}
//						}
//					}
//					if (same + deletion == windowSize && deletion <= deletionLevel)
//					{
//						Connection connectI(i, same);
//						Connection connectJ(j, same);
//						polaczenia[i].push_back(connectJ);
//						polaczenia[j].push_back(connectI);
//					}
//				}
//			}
			same = 0;
			deletion = 0;
		}
	}
}

//--------SZUKANIE KLIK
vector<Connection> findCommonNodes(vector<Connection> nodes, int vertex)
{
	vector<Connection> common;
	for (Connection n : nodes)
	{
		for (int i = 0; i < polaczenia[vertex].size(); i++)
		{
			if (polaczenia[vertex][i].vertexId == n.vertexId)
			{
				common.push_back(n);
			}
		}
	}
	return common;
}
void bronKerbosh(vector<Connection> potentialClique, vector<Connection> nodes, vector<Connection> skipedNodes)
{
	if (nodes.size() == 0 && skipedNodes.size() == 0)
	{
		kliki.push_back(potentialClique);
		//cout << "znalezione klike nr " << motiveCounter++ <<endl;
		return;
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		vector<Connection> newPotentialClique = potentialClique;
		newPotentialClique.push_back(nodes[i]);
		vector<Connection> newNodes = findCommonNodes(nodes, nodes[i].vertexId);
		vector<Connection> newSkippedNodes = findCommonNodes(skipedNodes, nodes[i].vertexId);

		bronKerbosh(newPotentialClique, newNodes, newSkippedNodes);

		newPotentialClique.pop_back();
		skipedNodes.push_back(nodes[i]);
		nodes.erase(remove(nodes.begin(), nodes.end(), nodes[i]), nodes.end());
	}
}
vector<Connection> findBestClique()
{
	int maxCliqSize = -1, maxCliqID = -1, score = 0;
	vector<int> wynikKliki;
	for(auto clq : kliki)
	{
		int sum = 1;
		for(auto x : clq)
		{
			sum += x.score;
		}
		wynikKliki.push_back(sum);
	}
	for (int i = 0; i < kliki.size(); i++)
	{
		
			int newScore = 0;
			for(int k=0;k<  kliki[i].size()-1;k++)
			{
				if (kliki[i][k + 1].vertexId - kliki[i][k].vertexId > windowSize &&
					krawedzie[kliki[i][k].vertexId]->seqId != krawedzie[kliki[i][k + 1].vertexId]->seqId)
					newScore++;
				else
					newScore--;
			}
			int klikScore = maxCliqID != -1 ? wynikKliki[maxCliqID] : 0;
			if (score < newScore && klikScore < wynikKliki[i]) {
				score = newScore;
				maxCliqSize = kliki[i].size();
				maxCliqID = i;
			}
		
	}
	if(maxCliqID == -1)
	{
		cout << "nie znaleziono nic sensownego" << endl << "zakonczenie pracy programu...." << endl;
		getchar();
		getchar();
		exit(0);
	}
	return kliki[maxCliqID];
}
vector<Connection> findClique()
{
	int maxDeg = 0, maxDegID;
	for (int i = 0; i < polaczenia.size(); i++)
	{
		if (polaczenia[i].size() > maxDeg)
		{
			maxDeg = polaczenia[i].size();
			maxDegID = i;
		}
	}
	vector<Connection> potentialClique, skipedNodes, wszystkiePolaczenia = polaczenia[maxDegID];
	wszystkiePolaczenia.push_back(Connection(maxDegID, windowSize));
	
	bronKerbosh(potentialClique, wszystkiePolaczenia, skipedNodes);
	//bronKerbosh(potentialClique, polaczenia[maxDegID], skipedNodes);
	
	return findBestClique();
}

//--------SZUKANIE MOTYWU
Motive findMotive(vector<Connection> clq)
{
	Motive motive;
	int idGenerator = 0;
	for (auto frag : clq)
	{
		motive.fragments.push_back(SequenceFragment(krawedzie[frag.vertexId]->seqId, krawedzie[frag.vertexId]->begin,
			krawedzie[frag.vertexId]->begin + windowSize, idGenerator++, krawedzie[frag.vertexId]->seq));
	}

	for (int i = 0; i < windowSize; i++)
	{
		map<char, int> consensusColumn;
		consensusColumn['A'] = 0;
		consensusColumn['C'] = 0;
		consensusColumn['T'] = 0;
		consensusColumn['G'] = 0;
		consensusColumn['_'] = 0;
		for (auto id : clq)
		{
			consensusColumn[krawedzie[id.vertexId]->seq[i]]++;
		}

		bool nucWasAdded = false;
		char nukleotyd='_';
		int iloscNukl = 0;
		for (auto nuc : Nukleotydy)
		{
			if (consensusColumn[nuc] >= clq.size() / 2 && consensusColumn[nuc] > iloscNukl)
			{
				iloscNukl = consensusColumn[nuc];
				nukleotyd = nuc;
			}
		}
		motive.consensus.push_back(nukleotyd);
	}
	cout << "znaleziono motyw: " << motive.consensus << endl;
	return motive;
}//do przerobeinia
Motive extendMotive(vector<SequenceFragment> fragments,bool directionLeft)
{
	int extendedWindow = 3;
	Motive motive;
	vector<string> clq;
	for (auto frag : fragments)
	{
		int ext=0;
		if (directionLeft)
		{
			ext = frag.begining - extendedWindow < 0 ? -1 : frag.begining - extendedWindow;
			if(ext>=0)
			{
				SequenceFragment extendedFragment(frag.seqId, ext, frag.begining-1, frag.id,"");
				for(int i= ext; i < frag.begining; i++)
				{
					extendedFragment.seq += sekwencje[frag.seqId].seq[i];
				}
				motive.fragments.push_back(extendedFragment);
			}
		}
		else {
			ext = frag.end + extendedWindow >= sekwencje[frag.seqId].seq.size() ? -1 : frag.end + extendedWindow;
			if (ext>=0)
			{
				SequenceFragment extendedFragment(frag.seqId, frag.end + 1, ext, frag.id,"");
				for (int i = frag.end + 1; i <= ext; i++)
				{
					extendedFragment.seq += sekwencje[frag.seqId].seq[i];
				}
				motive.fragments.push_back(extendedFragment);
			}
		}
				
	}
	if (motive.fragments.size() >= (fragments.size() / 2))
	{
		for (int i = 0; i < extendedWindow; i++)
		{
			map<char, int> consensusColumn;

			consensusColumn['A'] = 0;
			consensusColumn['C'] = 0;
			consensusColumn['T'] = 0;
			consensusColumn['G'] = 0;
			consensusColumn['_'] = 0;

			for (auto frag : motive.fragments)
			{
				consensusColumn[frag.seq[i]]++;
			}

			char nukleotyd = '_';
			int iloscNukl = fragments.size() *0.4f;

			for (auto nuc : Nukleotydy)
			{
				if ( consensusColumn[nuc] > iloscNukl)
				{
					iloscNukl = consensusColumn[nuc];
					nukleotyd = nuc;
				}
			}

			motive.consensus.push_back(nukleotyd);
			
		}
		return motive;
	}
	
	return Motive();	
}
Motive maximalizeMotive(Motive motive)
{
	int score = 0;
	cout << "poszerzanie motywu w lewo..." << endl;
	while (true)//lewo
	{
		score = 0;
		Motive extendedMotive = extendMotive(motive.fragments, true);
		if (extendedMotive.fragments.size() != 0) {
			for (int i = 0; i < 3; i++)
			{
				if (extendedMotive.consensus[i] != '_')score++;
			}
			if (score >= 2)
			{
				motive.consensus = extendedMotive.consensus + motive.consensus;
				for (int i=0;i<extendedMotive.fragments.size();i++)
				{
					for (int o = 0; o < motive.fragments.size();o++) {
						if (motive.fragments[o].id == extendedMotive.fragments[i].id) {
							motive.fragments[o].begining = extendedMotive.fragments[i].begining;
							motive.fragments[o].seq = extendedMotive.fragments[i].seq + motive.fragments[o].seq;
							
						}
					}
				}
			}
			else break;
		}
		else break;
	}
	score = 0;
	cout << "poszerzanie motywu w prawo...."<< endl;
	while (true)//prawo
	{
		score = 0;
		Motive extendedMotive = extendMotive(motive.fragments, false);
		if (extendedMotive.fragments.size() != 0) {
			for (int i = 0; i < 3; i++)
			{
				if (extendedMotive.consensus[i] != '_')score++;
			}
			if (score >= 2)
			{
				motive.consensus = motive.consensus + extendedMotive.consensus;
				for (int i = 0; i<extendedMotive.fragments.size(); i++)
				{
					for (int o = 0; o < motive.fragments.size(); o++) {
						if (motive.fragments[o].id == extendedMotive.fragments[i].id) {
							motive.fragments[o].end = extendedMotive.fragments[i].end;
							motive.fragments[o].seq += extendedMotive.fragments[i].seq;
							
						}
					}
				}

			}
			else break;
		}
		else break;
	}
	cout << motive.consensus<<"\t<-nowy motyw " << endl;
	SetConsoleTextAttribute(hConsole, 15);
	for (auto frag : motive.fragments) {
		for (int i = 0; i < frag.seq.size(); i++) {
			if(frag.seq[i] != motive.consensus[i] && motive.consensus[i] != '_')
				SetConsoleTextAttribute(hConsole, 4);
			else 
				if(motive.consensus[i] == '_')
					SetConsoleTextAttribute(hConsole, 11);
				else
					SetConsoleTextAttribute(hConsole, 15);
			cout << frag.seq[i];
		}
		cout << endl;
	}
	return motive;
}

//--------WYPISANIE MOTYWU
void print(SequenceFragment seq)
{
	SetConsoleTextAttribute(hConsole, 15);
	cout << endl<<seq.seqId << endl;
	for(int i=0;i<sekwencje[seq.seqId].seq.size();i++)
	{			
		if(seq.begining <= i && seq.end >= i)
		{
			SetConsoleTextAttribute(hConsole, 2);

			for(int j=0;j<seq.seq.size();j++)
			{
//				if(motyw.consensus[j] == '_' || motyw.consensus[j] != seq.seq[j]) SetConsoleTextAttribute(hConsole, 4);
//				else SetConsoleTextAttribute(hConsole, 2);
				cout << seq.seq[j];
				i++;
			}
		}

		SetConsoleTextAttribute(hConsole, 15);
		cout << sekwencje[seq.seqId].seq[i];
	}
	cout << endl;
}
