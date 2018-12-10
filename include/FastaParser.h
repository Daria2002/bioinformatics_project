#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <string>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
using namespace std;

class FastaParser
{

public:
	FastaParser(string input_file, int K);
	~FastaParser();

	unordered_set<string> parseKmers();
	vector<string> parseKmersToVector();
	vector<string> hittingSetSparsification();
	vector<string> edgeKmers();
	vector<string> bestIndexSparsification();
	
	void printKmers();

private:
	fstream fastaStream;
	int K;
	unordered_set<string> kmerSet;
	vector<string> kmerVector;
};

#endif
