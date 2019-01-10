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

	unordered_set<string> ParseKmers();
	vector<string> ParseKmersToVector();
	vector<string> HittingSetSparsification();
	vector<string> EdgeKmers();
	vector<string> BestIndexSparsification();
	
	void PrintKmers();

private:
	fstream fasta_stream;
	int K;
	unordered_set<string> kmer_set;
	vector<string> kmer_vector;
};

#endif
