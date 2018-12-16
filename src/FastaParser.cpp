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
		if (counter == 0)
			cout << "Parsing fasta file"<< endl;

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
		if (counter == 0)
			cout << "Parsing fasta file" << endl;

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
	}
	return edgeNeighbours;
}	

vector<string> FastaParser::bestIndexSparsification()
{
	string current_line;
	int counter = 0;
	int s = 1;

	vector<string> kmersToSave;
	bool alreadyInSet;
	while(getline(fastaStream, current_line))
	{
		if (current_line.at(0) == '>')
			continue;
		cout<<"doslo je do linije u fasta file koja je br:"<<counter<<endl;
		// save all kmers from first line
		if (counter==0)
		{
			const char* c_line = current_line.c_str();
			for (int i = 0; i <= current_line.size()-K; i=i+1+s)
			{ 
				alreadyInSet = false;
				char* newKmer = (char*)calloc(K, sizeof(char));
				strncpy(newKmer, c_line+i, K);
				string* newKmerString = new string(newKmer);
				//check if kmers are already saved
				for(string kmer:kmersToSave)
				{
					if(kmer==*newKmerString)
					{
						alreadyInSet = true;
						break;
					}
				}
				if(alreadyInSet==true)
					continue;
				kmersToSave.push_back(*newKmerString);
			}
		}
		//save kmers from best Kr
		//Kr is kmer set started at index r
		else if(counter>0)
		{	
			unordered_map<int, int> bestIndexMap;
			const char* c_line = current_line.c_str();
			int maxFreq = 0;
			int bestIndex = 0;

			for(int startIndex = 0; startIndex<=K+1;	startIndex=startIndex+1)
			{
				vector<string> kmersAtBestIndex;
				for (int i = startIndex; i <= current_line.size()-K; i=i+1+s)
				{ 
					char* newKmer = (char*)calloc(K, sizeof(char));
					strncpy(newKmer, c_line+i, K);
					string* newKmerString = new string(newKmer);
					kmersAtBestIndex.push_back(*newKmerString);
					for(string savedKmers : kmersToSave)
					{
						if(*newKmerString==savedKmers)
						{
							bestIndexMap[startIndex]++;
							break;
						}
					}
				}
				
				if(bestIndexMap[startIndex]>=maxFreq)
				{
					maxFreq = bestIndexMap[startIndex];
					bestIndex = startIndex;
				}
			}

			for(int i=bestIndex; i<=current_line.size()-K; i=i+s+1)
			{
				bool alreadyInSet = false;
				char* newKmer = (char*)calloc(K, sizeof(char));
				strncpy(newKmer, c_line+i, K);
				string* newKmerString = new string(newKmer);
				for(auto kmer : kmersToSave)
				{
					if(kmer == *newKmerString)
					{
						alreadyInSet = true;
						break;
					}
				}
				if(!alreadyInSet)
					kmersToSave.push_back(*newKmerString);
			}
		}
		counter++;
	}
	return kmersToSave;
}

//returns neighbours for every kmer in sequence
//key in map is kmer, value in map are neighbours for kmer
unordered_map<string, int> makeMapWithKmersAndNeighbours(string sequence, int K)
{
	unordered_map<string, int> neighbours;
	string leftNeighbour;
	string kmerHelp;
	const char* c_line = sequence.c_str();
	//first kmer doesn't have left neighbour
	char* kmer = (char*)calloc(K, sizeof(char));
	strncpy(kmer, c_line, K);
	string* kmerString = new string(kmer);
	leftNeighbour = *kmerString;

	//loop that adds left neighbours, starting from second kmer because first kmer doesn't have left neighbour
	for(int i=1; i<sequence.size()-K; i++)
	{
		char* kmer = (char*)calloc(K, sizeof(char));
		strncpy(kmer, c_line+i, K);
		string* kmerString = new string(kmer);
		//cout<<"leftNeighbour:"<<leftNeighbour<<" and current kmer: "<<*kmerString<<endl;
		neighbours[*kmerString]++;
	}

	//loop that adds right neighbours, starting from last kmers 
	for(int j=sequence.size()-K; j>0; j--)
	{
		char* rightNeighbour = (char*)calloc(K, sizeof(char));
		//last kmer
		strncpy(rightNeighbour, c_line+j, K);
		string* kmerString = new string(rightNeighbour);
		neighbours[*kmerString]++;
	}
	return neighbours;
}

vector<string> FastaParser::hittingSetSparsification()
{
	vector<string> kmersToSave;
	string current_line;
	while (getline(fastaStream, current_line))
	{
		unordered_map<string, int> kmersAndNeighbours;
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>')
			continue;

		kmersAndNeighbours = makeMapWithKmersAndNeighbours(current_line, K);
		const char* c_line = current_line.c_str();
		bool skip = false;
		vector<string> kmersToDelete;

		//kmers that appeared more than once are saved in kmersToDelete
		for ( auto kmer = kmersAndNeighbours.begin(); kmer != kmersAndNeighbours.end(); kmer++)
		{
			if(kmer->second>1)
				kmersToDelete.push_back(kmer->first);
		}

		//adding kmers that appeared only once in kmersToSave vector
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
