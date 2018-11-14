#include "FastaParser.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

FastaParser::FastaParser(std::string input_file, int K)
	: fastaStream(input_file),
	  K(K)
{
	cout << "FastaParser() - Initialized stream: " << input_file << endl;
}

FastaParser::~FastaParser()
{
}

void FastaParser::parseKmers()
{
	cout << "FastaParser::Parse() - starting to parse file..." << endl;

	string current_line;
	int counter = 0;
	while (getline(fastaStream, current_line))
	{
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>')
		{
			cout << "Skipping header: " << current_line << endl;
			continue;
		}

		// ... if not header start to parse line
		if (counter % 500 == 0)
			cout << "Parsing line " << counter << endl;

		// TODO: consider characters on the end of each line
		const char* c_line = current_line.c_str();
		for (int i = 0; i < current_line.size(); i += K)
		{
			char* newKmer = (char*)malloc(K * sizeof(char));
			strncpy(newKmer, c_line+i, K);
			//cout << newKmer << endl;

			// Try to insert it, if not successful, free the memory
			if (!kmerSet.insert(newKmer).second)
				free(newKmer);
		}

		counter++;
	}

	cout << "FastaParser::parse() - Finished parsing!" << endl;
}

void FastaParser::printKmers()
{
	cout << K << "-mers are :" << endl;
	int i = 0, kmer_size = kmerSet.size();
	for(char* kmer : kmerSet)
	{
		i++;
		cout << i << "/" << kmer_size << ":\t" << kmer << endl;
	}
}
