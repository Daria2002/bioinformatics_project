#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

#include <string>
#include <fstream>
#include <unordered_set>

class FastaParser
{

public:
	FastaParser(std::string input_file, int K);
	~FastaParser();

	void parseKmers();
	void printKmers();

private:
	std::fstream fastaStream;
	int K;
	std::unordered_set<char*> kmerSet;

};

#endif
