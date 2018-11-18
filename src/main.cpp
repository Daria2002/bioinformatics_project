#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <bf/all.hpp>
#include "FastaParser.h"
#include <stdlib.h>

using namespace std;

char randomKmerChar ()
{
	int randomACGT = rand()%4;
	switch (randomACGT)
	{
		case 0:
			return 'A';
	    case 1:
	        return 'C';
	    case 2:
	        return 'G';
	    case 3:
	    	return 'T';
	}
}

vector<string> generateTestKmerSet(unordered_set<string> kmerSet, int testSize, int K)
{
	vector<string> kmerSetVecor;
	int index = 0;
	cout << "Adding kmers to vector" << endl;
	for (string kmer : kmerSet)
	{
		kmerSetVecor.push_back(kmer);
		index++;
	}
	cout << "Finished adding kmers to vector" << endl;

	int randomKmerSetIndex, randomKmerIndex;
	std::vector<string> testKmers;

	cout << "kmerSetSize: " << kmerSet.size() << endl;
	cout << "K: " << K << endl;

	for (int i=0; i<testSize; i++)
	{
		randomKmerSetIndex = rand()%(kmerSet.size());
		randomKmerIndex = rand()%(K);

		cout<<"randomKmerSetIndex:"<<randomKmerSetIndex<<endl;
		cout<<"randomKmerIndex:"<<randomKmerIndex<<endl;

		string randomSelectedKmer = kmerSetVecor[randomKmerSetIndex];

		cout<<"randomSelectedKmer:"<<randomSelectedKmer<<endl;

		char randomChar;
		while(1)
		{
			randomChar = randomKmerChar();
			if(randomChar!=randomSelectedKmer[randomKmerIndex])
				break;
		}

		randomSelectedKmer[randomKmerIndex] = randomChar;
		cout<<"promjenjeni randomSelectedKmer:"<<randomSelectedKmer<<" char:"<<randomChar<<endl;
		testKmers.push_back(randomSelectedKmer);
		cout << randomSelectedKmer << endl;
	}
	cout << "Finished generating test kmers" << endl;
	return testKmers;
}

int main (int argc, char *argv[]) 
{
	srand(time(NULL));
	size_t numOfCells = 1024*1024*20;
	int numOfHashes = 2;
	int testSetSize = 20;
	//bool bloomFilterResult = (bool)malloc(testSetSize*sizeof(bool));
	std::unordered_set<string> kmerSet;
	vector<string> kmerSetTest;
	if(argc < 3) {
		cerr << "Please write two arguments: path to fasta file an number of k-mers"
			<< endl;
		return -1;
	}

	int K = atoi(argv[2]);
	string fastaFile = argv[1];

	FastaParser fp(fastaFile, K);
	kmerSet = fp.parseKmers();
	

	cout<<"*******Classic Bloom Filter*******"<<endl;

	bf::basic_bloom_filter bloomFilter(
		bf::make_hasher(numOfHashes), numOfCells);
	for (auto kmer : kmerSet)
		bloomFilter.add(kmer);
	cout << "Starting generating test kmers." << endl;
	kmerSetTest = generateTestKmerSet(kmerSet, testSetSize, K);
	size_t result;
	
	std::vector<bool> bloomFilterResult;

	for (auto kmer : kmerSetTest){
		result = bloomFilter.lookup(kmer);
		bloomFilterResult.push_back((bool)result);
	}
	
	int br = 0;
	for(bool res : bloomFilterResult)
	{
		br++;
		cout << br<<": " << res << endl;
	}


	return 0;
}