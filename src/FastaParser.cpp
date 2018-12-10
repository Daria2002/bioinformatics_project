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
		
		for (int i = 0; i < current_line.size()-K+1; i += 1)
		{ 
			char* newKmer2 = (char*)calloc(K, sizeof(char));
			strncpy(newKmer2, c_line+i, K);
			string* newKmerString2 = new string(newKmer2);
			// Try to insert it, if not successful, free the memory
			if (!kmerSet.insert(*newKmerString2).second){
				free(newKmer2);
				free(newKmerString2);
			}
		}

		counter++;
	}
	cout << "FastaParser::parse() - Finished parsing!" << endl;
	return kmerSet;
}

vector<string> FastaParser::parseKmersToVector()
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
		if (counter % 3500 == 0)
			cout << "Parsing line... " << counter << endl;

		const char* c_line = current_line.c_str();
		
		for (int i = 0; i < current_line.size()-K+1; i += 1)
		{ 
			char* newKmer = (char*)calloc(K, sizeof(char));
			strncpy(newKmer, c_line+i, K);
			string* newKmerString = new string(newKmer);
			// Try to insert it, if not successful, free the memory
			kmerVector.push_back(*newKmerString);
		}

		counter++;
	}
	cout << "FastaParser::parse() - Finished parsing!" << endl;
	return kmerVector;
}

std::vector<string> FastaParser::edgeKmers()
{
	vector<string> edgeNeighbours;
	string current_line;
	int counter = 0;
	while (getline(fastaStream, current_line))
	{
		unordered_map<std::string, int> kmerMap;
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>')
			continue;

		const char* c_line = current_line.c_str();

		//Counting how many times each kmer accurate in sequence
		for (int i = 0; i < current_line.size()-K+1; i++)
		{ 
			char* kmer = (char*)calloc(K, sizeof(char));
			strncpy(kmer, c_line+i, K);
			string* kmerString = new string(kmer);
			if(i==current_line.size()-K || i==0)
				edgeNeighbours.push_back(*kmerString);
		}
		
		counter++;
	}

	return edgeNeighbours;
}	

std::vector<string> FastaParser::bestIndexSparsification()
{
	string current_line;
	int counter = 0;
	int s = 1;
	std::vector<string> kmersToSave;
	while(getline(fastaStream, current_line))
	{
		if (current_line.at(0) == '>')
			continue;
		// save all kmers from first line
		if (counter==0)
		{
			const char* c_line = current_line.c_str();
			for (int i = 0; i <= current_line.size()-K; i=i+1+s)
			{ 
				char* newKmer = (char*)calloc(K, sizeof(char));
				strncpy(newKmer, c_line+i, K);

				string* newKmerString = new string(newKmer);

				kmersToSave.push_back(*newKmerString);
			}
		}
		//save kmers from best Kr
		else if(counter>0)
		{	
			unordered_map<int, int> bestIndexMap;
			const char* c_line = current_line.c_str();
			for(int startIndex = 0; startIndex<=current_line.size()-K;	startIndex=startIndex+1)
			{
				std::vector<string> kmersAtBestIndex;
				for (int i = startIndex; i <= current_line.size()-K; i=i+1+s)
				{ 
					char* newKmer = (char*)calloc(K, sizeof(char));
					strncpy(newKmer, c_line+i, K);

					string* newKmerString = new string(newKmer);

					kmersAtBestIndex.push_back(*newKmerString);
				}

				for(string kmer : kmersAtBestIndex)
				{
					for(string savedKmers : kmersToSave)
					{
						if(kmer==savedKmers)
						{
							bestIndexMap[startIndex]++;
							break;
						}
					}
				}
			}

			int maxFreq = 0;
			int bestIndex = 0;
			for (auto index = bestIndexMap.begin(); index != bestIndexMap.end(); index++)
			{
				if(index->second>maxFreq)
				{
					maxFreq = index->second;
					bestIndex = index->first;
				}
			}

			std::vector<string> kmersAtBestIndex;
			for(int i=bestIndex; i<=current_line.size()-K; i=i+1)
			{
				char* newKmer = (char*)calloc(K, sizeof(char));
				strncpy(newKmer, c_line+i, K);

				string* newKmerString = new string(newKmer);

				kmersAtBestIndex.push_back(*newKmerString);
			}
		}

		counter++;
	}
}

vector<string> FastaParser::hittingSetSparsification()
{
	vector<string> kmersToSave;
	string current_line;
	int counter = 0;
	while (getline(fastaStream, current_line))
	{
		unordered_map<std::string, int> kmerMap;
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

		//Counting how many times each kmer accurate in sequence
		for (int i = 0; i < current_line.size()-K+1; i++)
		{ 
			char* kmer = (char*)calloc(K, sizeof(char));
			strncpy(kmer, c_line+i, K);
			string* kmerString = new string(kmer);
			kmerMap[*kmerString]++;
		}
		bool skip = false;
		vector<string> kmersToDelete;
		for ( auto kmer = kmerMap.begin(); kmer != kmerMap.end(); kmer++)
		{
			if(kmer->second>1)
				kmersToDelete.push_back(kmer->first);
		}

		for (int i = 0; i < current_line.size()-K+1; i++)
		{ 
			skip = false;
			char* newKmer = (char*)calloc(K, sizeof(char));
			strncpy(newKmer, c_line+i, K);

			string* newKmerString = new string(newKmer);

			for(string kmerToDelete : kmersToDelete)
			{
				if(*newKmerString == kmerToDelete)
				{
					skip = true;
					break;
				}
			}
			if(skip)
				continue;
			kmersToSave.push_back(*newKmerString);
		}
		counter++;
	}

	return kmersToSave;
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
