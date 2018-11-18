#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <string>
#include <fstream>
#include <unordered_set>
using namespace std;

class FastaParser
{

public:
	FastaParser(string input_file, int K);
	~FastaParser();

	unordered_set<string> parseKmers();
	void printKmers();

private:
	fstream fastaStream;
	int K;
	unordered_set<string> kmerSet;

};

#endif
