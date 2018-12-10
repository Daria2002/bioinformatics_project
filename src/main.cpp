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

// Single sequence sparsification
vector<string> singleSequenceSparsification(vector<string> kmerSet, int s)
{
	vector<string> kmerSetSparse;
	int i=0;
	for(auto kmer : kmerSet)
	{
		if(i%(s+1)==0)
		{
			kmerSetSparse.push_back(kmer);
		}
		i++;
	}

	for (auto kmer : kmerSet)
		cout<<"kmer:"<<kmer<<endl;
	for(auto kmer : kmerSetSparse)
		cout<<"Sparse:"<<kmer<<endl;

	return kmerSetSparse; 
}

vector<string> makeRightNeighbours(string kmer)
{
	string firstCharCombinations = "ACTG";
	vector<string> rightNeighbours;
	string makeNeighbour;
	for(char possibleChar : firstCharCombinations)
	{	
		makeNeighbour = kmer[1];
		for (int i=2; i<kmer.size(); i++)
		{
			makeNeighbour+=kmer[i];
		}

		makeNeighbour += possibleChar;
		rightNeighbours.push_back(makeNeighbour);
	}
	return rightNeighbours;
}

vector<string> makeLeftNeighbours(string kmer)
{
	string firstCharCombinations = "ACTG";
	vector<string> leftNeighbours;
	string makeNeighbour;
	for(char possibleChar : firstCharCombinations)
	{	
		makeNeighbour = possibleChar;
		for (int i=0; i<kmer.size()-1; i++)
		{
			makeNeighbour+=kmer[i];
		}
		leftNeighbours.push_back(makeNeighbour);
	}/*
	cout<<"kmer"<<kmer<<endl;
	for (string el : leftNeighbours)
	{
		cout<<"element u lijevim susjedima:"<<el<<endl;
	}*/
	return leftNeighbours;
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
/*
		cout<<"randomKmerSetIndex:"<<randomKmerSetIndex<<endl;
		cout<<"randomKmerIndex:"<<randomKmerIndex<<endl;
*/
		string randomSelectedKmer = kmerSetVecor[randomKmerSetIndex];
		//cout<<"randomSelectedKmer:"<<randomSelectedKmer<<endl;

		char randomChar;
		while(1)
		{
			randomChar = randomKmerChar();
			if(randomChar!=randomSelectedKmer[randomKmerIndex])
				break;
		}

		randomSelectedKmer[randomKmerIndex] = randomChar;
		//cout<<"promjenjeni randomSelectedKmer:"<<randomSelectedKmer<<" char:"<<randomChar<<endl;
		testKmers.push_back(randomSelectedKmer);
		//cout << randomSelectedKmer << endl;
	}
	cout << "Finished generating test kmers" << endl;

	return testKmers;
}


float falsePositiveRate(vector<bool> bloomFilterResult, vector<bool> bloomFilterResultReal)
{
	float FPrate = 0;

	for (int i = 0; i < bloomFilterResult.size(); i++)
	{
		if(!(bloomFilterResult[i] == bloomFilterResultReal[i]))
			FPrate++;
	}

	FPrate = (FPrate / bloomFilterResult.size())*100 ;
	return FPrate;
}

vector<bool> compareTestKmerWithSavedKmers(unordered_set<string> kmerSet, vector<string> kmerSetTest)
{
	vector<bool> bloomFilterResultReal;
	bool setted = false;
	for(auto kmer : kmerSetTest)
	{
		for(auto kmerSaved : kmerSet)
		{
			if(kmer == kmerSaved)
			{
				bloomFilterResultReal.push_back((bool)true);
				setted = true;
				break;
			}
		}
		if(!setted)
		{
			bloomFilterResultReal.push_back((bool)false);
		}
		else
			setted = false;
	}
	return bloomFilterResultReal;
}

int main (int argc, char *argv[]) 
{
	srand(time(NULL));
	size_t numOfCells = 500;
	//size_t numOfCells = 100*10000;
	int numOfHashes = 2;
	int testSetSize = 20;

	//bool bloomFilterResult = (bool)malloc(testSetSize*sizeof(bool));  
	unordered_set<string> kmerSet;
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

	bf::basic_bloom_filter bloomFilter(
		bf::make_hasher(numOfHashes), numOfCells);

	for (auto kmer : kmerSet)
		bloomFilter.add(kmer);


	cout << "Starting generating test kmers." << endl;
	kmerSetTest = generateTestKmerSet(kmerSet, testSetSize, K);

	cout<<"*******Classic Bloom filter*******"<<endl;
	size_t result;
	vector<bool> bloomFilterResult;
	vector<bool> bloomFilterResultReal;

	int start_s=clock();

	for (auto kmer : kmerSetTest)
	{
		result = bloomFilter.lookup(kmer);
		bloomFilterResult.push_back((bool)result);
	}

	int stop_s=clock();

	float timeClassicBloomFilter = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;
	cout << "Classic Bloom Filter-time: " << timeClassicBloomFilter << " s" << endl;

	bloomFilterResultReal = compareTestKmerWithSavedKmers(kmerSet, kmerSetTest);

	float FPrateClassicBloomFilter;
	FPrateClassicBloomFilter = falsePositiveRate(bloomFilterResult, bloomFilterResultReal);
	cout << "Classic Bloom filter-fp rate:" << FPrateClassicBloomFilter << "%" << endl;

	cout<<"*******One-sided k-mer Bloom filter*******"<<endl;
	vector<bool> oneSidedResult;
	string *kmerToChange;
	vector<string> leftNeighbours;
	vector<string> rightNeighbours;
	bool neighbourInSet;

	start_s=clock();
	
	for(int i=0; i<bloomFilterResult.size(); i++)
	{
		if(bloomFilterResult[i])
		{
			kmerToChange = &kmerSetTest[i];
			//cout<<"kmerToChange:"<<*kmerToChange<<endl;
			
			leftNeighbours = makeLeftNeighbours(*kmerToChange);
			rightNeighbours = makeRightNeighbours(*kmerToChange);
			
			for(int j=0; j<4; j++)
			{
				if(bloomFilter.lookup(leftNeighbours[j]) || bloomFilter.lookup(rightNeighbours[j]))
				{
					oneSidedResult.push_back(true);
					neighbourInSet = true;
					break;
				}
			}
			if(!neighbourInSet)
			{
				oneSidedResult.push_back(false);
			}	
			else
				neighbourInSet = false;
		}
		else
			oneSidedResult.push_back(false);	
	}
	stop_s=clock();

	float timeOneSided = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 + timeClassicBloomFilter;
	cout << "One-sided Bloom Filter-time: " << timeOneSided << " s" << endl;

	float FPrateOneSided;
	FPrateOneSided = falsePositiveRate(oneSidedResult, bloomFilterResultReal);
	cout << "One-sided Bloom filter-fp rate:" << FPrateOneSided << "%" << endl;

	// Edge kmers
	std::vector<string> edgeKmersSet;
	FastaParser fastaParserEdgeKmers(fastaFile, K);
	edgeKmersSet = fastaParserEdgeKmers.edgeKmers();

	cout<<"*******Two-sided k-mer Bloom filter*******"<<endl;
	vector<bool> twoSidedResult;

	neighbourInSet = false;
	start_s=clock();
	for(int i=0; i<bloomFilterResult.size(); i++)
	{
		if(bloomFilterResult[i])
		{
			kmerToChange = &kmerSetTest[i];
			leftNeighbours = makeLeftNeighbours(*kmerToChange);
			rightNeighbours = makeRightNeighbours(*kmerToChange);
			
			for(auto leftNeighbour : leftNeighbours)
			{
				for(auto rightNeighbour : rightNeighbours)
				{
					if(bloomFilter.lookup(rightNeighbour) && bloomFilter.lookup(leftNeighbour))
					{
						twoSidedResult.push_back(true);
						neighbourInSet = true;
						break;	
					}

					else if(bloomFilter.lookup(rightNeighbour) || bloomFilter.lookup(leftNeighbour))
					{
						for(string edgeKmer : edgeKmersSet)
						{
							if(*kmerToChange==edgeKmer)
							{
								twoSidedResult.push_back(true);
								neighbourInSet = true;
								break;
							}
						}
					}
				}
				if(neighbourInSet)
					break;
			}
			if(!neighbourInSet)
				twoSidedResult.push_back(false);
			else
				neighbourInSet = false;
		}
		else
			twoSidedResult.push_back(false);	
	}
	stop_s=clock();

	float timeTwoSided = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 + timeClassicBloomFilter;
	cout << "Two-sided Bloom Filter-time: " << timeTwoSided << " s" << endl;

	float FPrateTwoSided;
	FPrateTwoSided = falsePositiveRate(twoSidedResult, bloomFilterResultReal);
	cout << "Two-sided Bloom filter-fp rate:" << FPrateTwoSided << "%" << endl;
	
	cout<<"*******Sparse: Best index match per reat sparsification*******"<<endl;
	vector<bool> bestIndexSetResult;
	//test set doesn't change
	std::vector<string> indexSet;
	FastaParser fpIndexSet(fastaFile, K);
	indexSet = fpIndexSet.bestIndexSparsification();
	
	return 0;
}