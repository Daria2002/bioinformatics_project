#include "FastaParser.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

FastaParser::FastaParser(std::string input_file, int K)
	: fastaStream(input_file),
	  K(K)
{
	cout << "FastaParser() - Initialized stream: " << input_file << endl;
}

FastaParser::~FastaParser()
{
}

std::unordered_set<string> FastaParser::parseKmers()
{
	cout << "FastaParser::Parse() - starting to parse file..." << endl;
	char* help = (char*)calloc((K-1), sizeof(char));
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
		if (counter % 3500 == 0)
			cout << "Parsing line... " << counter << endl;

		const char* c_line = current_line.c_str();
		
		if(counter>0)
		{
			for (int j = 1; j<K; j+=1)
			{
				char* newKmer = (char*)calloc(K, sizeof(char));
				strncpy(newKmer, help+j-1, K-j);
				strncat(newKmer, c_line, j);
				string* newKmerString = new string(newKmer);
				if (!kmerSet.insert(*newKmerString).second)
				{
					free(newKmer);
					free(newKmerString);
				}
			}
		}
		

		for (int i = 0; i < current_line.size()-K+1; i += 1)
		{ 
			char* newKmer2 = (char*)calloc(K, sizeof(char));
			strncpy(newKmer2, c_line+i, K);
			//cout << newKmer2 << endl;
			string* newKmerString2 = new string(newKmer2);
			// Try to insert it, if not successful, free the memory
			//cout<<"kmer:"<<*newKmerString2<<endl;
			if (!kmerSet.insert(*newKmerString2).second){
				free(newKmer2);
				free(newKmerString2);
			}
		}

		strncpy(help, c_line+current_line.size()-K+1, K-1);
		counter++;
	}
	cout << "FastaParser::parse() - Finished parsing!" << endl;
	return kmerSet;
}

void FastaParser::printKmers()
{
	cout << K << "-mers are :" << endl;
	int i = 0, kmer_size = kmerSet.size();
	for(string kmer : kmerSet)
	{
		i++;
		cout << i << "/" << kmer_size << ":\t" << kmer << endl;
	}
}
