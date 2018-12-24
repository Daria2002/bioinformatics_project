#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <bf/all.hpp>
#include "FastaParser.h"
#include <stdlib.h>
#include <math.h>

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
// Select every sth k-mer
vector<string> singleSequenceSparsification(vector<string> kmerSet, int s)
{
	vector<string> kmerSetSparse;
	int i=0;
	for(auto kmer : kmerSet)
	{
		if(i%(s+1)==0)
			kmerSetSparse.push_back(kmer);
		i++;
	}
	return kmerSetSparse; 
}

bool checkEdgeKmer (vector<string> edgeKmersSet, string kmer)
{
	for(auto edgeKmer:edgeKmersSet)
	{
		if(edgeKmer == kmer)
			return true;
	}
	return false;
}

bool decidePresent(string kmer, bool contains_left, bool contains_right, vector<string> edgeKmersSet)
{
	if(contains_left || contains_right)
	{
		if(contains_right && contains_left)
			return true;
		else
			return checkEdgeKmer(edgeKmersSet, kmer);
	}
	return false;
}

bool containsSet(std::vector<string> kmerSet, const bf::basic_bloom_filter &bloomFilter)
{
	for(string kmer:kmerSet)
	{
		if(bloomFilter.lookup(kmer))
			return true; 
	}
	return false;
}

string convertTo4Base(int numberInDecimal, int length)
{
	string numberIn4BaseInverse;
	string numberIn4Base;
	int remainder;
	for(int i = 0; i<length; i++)
	{
		while(numberInDecimal>0)
		{
			remainder = numberInDecimal%4;
			numberInDecimal = numberInDecimal/4;
			numberIn4BaseInverse = numberIn4BaseInverse + std::to_string(remainder);
		}
		for(int j=numberIn4BaseInverse.size(); j<length; j++)
			numberIn4BaseInverse = numberIn4BaseInverse + '0';
	}
	
	for(int i=numberIn4BaseInverse.size()-1; i>=0; i--)
		numberIn4Base += numberIn4BaseInverse.at(i);
	return numberIn4Base;
}

std::vector<string> prefixOrSuffixCombinations(int length)
{
	int numberOfCombinations = pow(4, length);
	std::vector<string> combinations;
	string base4Number;
	for(int i = 0; i<numberOfCombinations; i++)
	{
		base4Number = convertTo4Base(i, length);
		string combination;
		for(auto el:base4Number)
		{
			if(el=='0')
				combination += 'A';
			else if(el == '1')
				combination += 'C';
			else if(el == '2')
				combination += 'T';
			else if(el == '3')
				combination += 'G';
		}
		combinations.push_back(combination);
	}
	return combinations;
}

std::vector<string> sDistantLeftNeighbourSet(string kmer, int leftDist)
{
	std::vector<string> neighbours;

	for(int i = 1; i<=leftDist; i++)
	{
		// make prefixSet and base
		// prefix set are all combinations of prefixes 
		// base consists of first k-i element in kmer
		std::vector<string> prefixSet = prefixOrSuffixCombinations(i);
		string base = kmer.substr(0, kmer.size()-i);
		for(auto prefix:prefixSet)
			neighbours.push_back(prefix + base);
	}
	return neighbours;
}

std::vector<string> sDistantRightNeighbourSet(string kmer, int rightDist)
{
	std::vector<string> neighbours;
	for(int i = 1; i<=rightDist; i++) 
	{
		std::vector<string> suffixSet = prefixOrSuffixCombinations(i);
		string base = kmer.substr(i, kmer.size()-i);
		for(auto suffix:suffixSet)
			neighbours.push_back(base+suffix);
	}
	return neighbours;
}

bool relaxedContainsNeighbours(string kmer, int leftDist, int rightDist, const bf::basic_bloom_filter &bloomFilter, std::vector<string> edgeKmersSet)
{
	bool containsLeft=false;
	bool containsRight=false;
	for (int i = 0; i <=leftDist; i++)
	{
		// check if bf contains at least one left neighbour
		if(containsSet(sDistantLeftNeighbourSet(kmer, i), bloomFilter))
		{
			containsLeft = true;
			break;
		}
	}
	for(int j=0; j<=rightDist; j++)
	{
		// check if bf contains at least on right neighbour
		if(containsSet(sDistantRightNeighbourSet(kmer, j), bloomFilter))
		{
			containsRight = true;
			break;
		}
	}
	return decidePresent(kmer, containsLeft, containsRight, edgeKmersSet);
}

vector<bool> relaxedContains(vector<string> kmerTestSet, vector<string> edgeKmersSet, const bf::basic_bloom_filter &bloomFilter, int s)
{
	vector<bool> relaxedResults;
	bool kmerSaved=false;

	for(auto kmer : kmerTestSet)
	{
		kmerSaved = false;
		// if kmer is saved, search if kmer is in edgeKmersSet and check if neighbour/s are saved 
		if(bloomFilter.lookup(kmer))
			relaxedResults.push_back(relaxedContainsNeighbours(kmer, s+1, s+1, bloomFilter, edgeKmersSet));
		else
		{
			for(int i=0; i <= s; i++)
			{
				if(relaxedContainsNeighbours(kmer, i, s-i+1, bloomFilter,edgeKmersSet))
				{
					relaxedResults.push_back(true);
					kmerSaved=true;
					break;
				}
			}
			if(!kmerSaved)
				relaxedResults.push_back(false);
		}
	}
	return relaxedResults;
}

// Returns all combinations for right neighbours
// e.g. if kmer is ACTG, right neighbours can be:
// CTGX, where x={A,G,C,T}  
vector<string> makeRightNeighbours(string kmer)
{
	string firstCharCombinations = "ACTG";
	vector<string> rightNeighbours;
	string makeNeighbour;
	for(char possibleChar : firstCharCombinations)
	{	
		makeNeighbour = kmer[1];
		for (int i=2; i<kmer.size(); i++)
			makeNeighbour+=kmer[i];

		makeNeighbour += possibleChar;
		rightNeighbours.push_back(makeNeighbour);
	}
	return rightNeighbours;
}

// Returns all combinations for left neighbour
vector<string> makeLeftNeighbours(string kmer)
{
	string firstCharCombinations = "ACTG";
	vector<string> leftNeighbours;
	string makeNeighbour;
	for(char possibleChar : firstCharCombinations)
	{	
		makeNeighbour = possibleChar;
		for (int i=0; i<kmer.size()-1; i++)
			makeNeighbour+=kmer[i];
		
		leftNeighbours.push_back(makeNeighbour);
	}
	return leftNeighbours;
}

// returns test set
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
		string randomSelectedKmer = kmerSetVecor[randomKmerSetIndex];

		char randomChar;
		while(1)
		{
			randomChar = randomKmerChar();
			if(randomChar!=randomSelectedKmer[randomKmerIndex])
				break;
		}

		randomSelectedKmer[randomKmerIndex] = randomChar;
		testKmers.push_back(randomSelectedKmer);
	}
	cout << "Finished generating test kmers" << endl;

	return testKmers;
}

// calculate false positive rate
float falsePositiveRate(vector<bool> bloomFilterResult, vector<bool> bloomFilterResultReal)
{
	float FPrate = 0;

	for (int i = 0; i < bloomFilterResult.size(); i++)
	{
		if(!(bloomFilterResult[i] == bloomFilterResultReal[i]))
		{
			if(bloomFilterResult[i]==0 and bloomFilterResultReal[i]==1)
				cout<<"real=1, bf=0"<<endl;
			FPrate++;
		}
	}

	FPrate = (FPrate / bloomFilterResult.size())*100 ;
	return FPrate;
}

//check if kmers in kmerSetTest are in kmerSet
//kmerSet is set of all kmers that really exists in fasta file
//kmerSet is seto of test kmers
//function returns bool vector, where false means that kmer doesn't exist in fasta file and true that it exists
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
			bloomFilterResultReal.push_back((bool)false);
		else
			setted = false;
	}
	return bloomFilterResultReal;
}

int main (int argc, char *argv[]) 
{
	srand(time(NULL));
	size_t numOfCells = 90000;
	//size_t numOfCells = 100*10000;
	int numOfHashes = 2;
	int testSetSize = 100;

	//bool bloomFilterResult = (bool)malloc(testSetSize*sizeof(bool));  
	unordered_set<string> kmerSet;
	vector<string> kmerSetTest;
	if(argc < 4) {
		cerr << "Please write three arguments: path to fasta file, number of k-mers and output file"
			<< endl;
		return -1;
	}

	int K = atoi(argv[2]);
	string fastaFile = argv[1];
	string outputFile = argv[3];

	ofstream outputFileStream;
	outputFileStream.open(outputFile);

	FastaParser fp(fastaFile, K);
	kmerSet = fp.parseKmers();

	bf::basic_bloom_filter *bloomFilter;
	bloomFilter = new bf::basic_bloom_filter(bf::make_hasher(numOfHashes), numOfCells);

	for (auto kmer : kmerSet)
		bloomFilter->add(kmer);

	cout << "Starting generating test kmers." << endl;
	kmerSetTest = generateTestKmerSet(kmerSet, testSetSize, K);
	cout<<"**********************************"<<endl;
	cout<<"*******Classic Bloom filter*******"<<endl;
	cout<<"**********************************"<<endl;
	size_t result;
	vector<bool> bloomFilterResult;
	vector<bool> bloomFilterResultReal;

	int start_s=clock();

	for (auto kmer : kmerSetTest)
		bloomFilterResult.push_back((bool)bloomFilter->lookup(kmer));

	int stop_s=clock();

	float timeClassicBloomFilter = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;

	bloomFilterResultReal = compareTestKmerWithSavedKmers(kmerSet, kmerSetTest);

	float FPrateClassicBloomFilter;
	FPrateClassicBloomFilter = falsePositiveRate(bloomFilterResult, bloomFilterResultReal);

	cout << "Classic Bloom Filter-time: " << timeClassicBloomFilter << " s" << endl;
	cout << "Classic Bloom filter-fp rate:" << FPrateClassicBloomFilter << "%" << endl;

	outputFileStream << "Classic Bloom Filter-time: " << timeClassicBloomFilter << " s" << endl;
	outputFileStream << "Classic Bloom filter-fp rate:" << FPrateClassicBloomFilter << "%" << endl;

	cout<<"******************************************"<<endl;
	cout<<"*******One-sided k-mer Bloom filter*******"<<endl;
	cout<<"******************************************"<<endl;
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
				if(bloomFilter->lookup(leftNeighbours[j]) || bloomFilter->lookup(rightNeighbours[j]))
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

	outputFileStream << "One-sided Bloom Filter-time: " << timeOneSided << " s" << endl;
	outputFileStream << "One-sided Bloom filter-fp rate:" << FPrateOneSided << "%" << endl;

	cout<<"******************************************"<<endl;
	cout<<"*******Two-sided k-mer Bloom filter*******"<<endl;
	cout<<"******************************************"<<endl;
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
					if(bloomFilter->lookup(rightNeighbour) && bloomFilter->lookup(leftNeighbour))
					{
						twoSidedResult.push_back(true);
						neighbourInSet = true;
						break;	
					}
					// if only one neighbour is in set, check if kmer is in edgeKmersSet
					// kmers in edgeKmersSet should have only one neighbour saved
					else if(bloomFilter->lookup(rightNeighbour) || bloomFilter->lookup(leftNeighbour))
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
						if(neighbourInSet)
							break;
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
	int s = 1;

	outputFileStream << "Two-sided Bloom Filter-time: " << timeTwoSided << " s" << endl;
	outputFileStream << "Two-sided Bloom filter-fp rate:" << FPrateTwoSided << "%" << endl;

	cout<<"**************************************************************"<<endl;
	cout<<"*******Sparse: Best index match per read sparsification*******"<<endl;
	cout<<"**************************************************************"<<endl;

	vector<bool> bestIndexSetResult;
	//test set doesn't change
	std::vector<string> indexSet;
	FastaParser fpIndexSet(fastaFile, K);
	indexSet = fpIndexSet.bestIndexSparsification();
	bf::basic_bloom_filter *bloomFilterSparse;
	bloomFilterSparse = new bf::basic_bloom_filter(bf::make_hasher(numOfHashes), numOfCells);
	
	for (auto kmer : indexSet)
		bloomFilterSparse->add(kmer);

	std::vector<bool> sparseResults;
	start_s=clock();

	sparseResults = relaxedContains(kmerSetTest, edgeKmersSet, *bloomFilterSparse, s);

	stop_s=clock();

	float timeSparse = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;
	cout << "Sparse Bloom Filter-time: " << timeSparse << " s" << endl;
	// kmerSet contains kmers parsed from fasta file
	//checking which kmers from kmerSetTest are in kmerSet
	bloomFilterResultReal = compareTestKmerWithSavedKmers(kmerSet, kmerSetTest);

	float FPrateSparse;
	FPrateSparse = falsePositiveRate(sparseResults, bloomFilterResultReal);
	cout << "Sparse Bloom filter-fp rate:" << FPrateSparse << "%" << endl;

	outputFileStream << "Sparse Bloom Filter-time: " << timeSparse << " s" << endl;
	outputFileStream << "Sparse Bloom filter-fp rate:" << FPrateSparse << "%" << endl;

	cout<<"***************************************************************"<<endl;
	cout<<"***Sparse:Spasification via approximate hitting set, Relaxed***"<<endl;
	cout<<"***************************************************************"<<endl;

	//test set doesn't change
	std::vector<string> hittingSet;
	FastaParser fpHittingSet(fastaFile, K);
	hittingSet = fpHittingSet.hittingSetSparsification();

	bf::basic_bloom_filter *bloomFilterHittingSet;
	bloomFilterHittingSet = new bf::basic_bloom_filter(bf::make_hasher(numOfHashes), numOfCells);

	for (auto kmer : hittingSet)
		bloomFilterHittingSet->add(kmer);

	vector<bool> hittingSetRelaxedResults;
	vector<bool> hittingSetStrictResults;
	start_s=clock();

	hittingSetRelaxedResults = relaxedContains(kmerSetTest, edgeKmersSet, *bloomFilterHittingSet, s);
	//hittingSetStrictResults = strictContains(kmerSetTest, edgeKmersSet, bloomFilterHittingSet);

	stop_s=clock();

	float timeHittingSet = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;
	cout << "Hitting set Bloom Filter-time: " << timeHittingSet << " s" << endl;

	bloomFilterResultReal = compareTestKmerWithSavedKmers(kmerSet, kmerSetTest);

	float FPrateHittingSet;

	FPrateHittingSet = falsePositiveRate(hittingSetRelaxedResults, bloomFilterResultReal);
	cout << "Hitting set Bloom filter-fp rate:" << FPrateHittingSet << "%" << endl;
	
	outputFileStream << "Hitting set Bloom Filter-time: " << timeHittingSet << " s" << endl;
	outputFileStream << "Hitting set Bloom filter-fp rate:" << FPrateHittingSet << "%" << endl;

	cout<<"***************************************************************"<<endl;
	cout<<"********Sparse: Single Sequence Sparsification, Relaxed********"<<endl;
	cout<<"***************************************************************"<<endl;
	vector<bool> singleSequenceSparsificationSetResult;
	//test set doesn't change
	vector<string> singleSequenceSparsificationSet;
	FastaParser fpSequenceSparsificationSet(fastaFile, K);
	vector<string> kmerSetVecor;
	kmerSetVecor = fpSequenceSparsificationSet.parseKmersToVector();
	singleSequenceSparsificationSet = singleSequenceSparsification(kmerSetVecor, s);
	bf::basic_bloom_filter *bloomFilterSequenceSparsificationSet;
	bloomFilterSequenceSparsificationSet = new bf::basic_bloom_filter(bf::make_hasher(numOfHashes), numOfCells);

	for (auto kmer : singleSequenceSparsificationSet)
		bloomFilterSequenceSparsificationSet->add(kmer);

	start_s=clock();

	singleSequenceSparsificationSetResult = relaxedContains(kmerSetTest, edgeKmersSet, *bloomFilterSequenceSparsificationSet, s);
	
	stop_s=clock();

	float timeSequenceSparsification = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;

	float FPrateSequenceSparsification;
	FPrateSequenceSparsification = falsePositiveRate(singleSequenceSparsificationSetResult, bloomFilterResultReal);
	
	cout << "Sequence Sparsification Bloom Filter-time: " << timeSequenceSparsification << " s" << endl;
	cout << "Sequence sparsification Bloom filter-fp rate:" << FPrateSequenceSparsification << "%" << endl;

	outputFileStream << "Sequence Sparsification Bloom Filter-time: " << timeSequenceSparsification << " s" << endl;
	outputFileStream << "Sequence sparsification Bloom filter-fp rate:" << FPrateSequenceSparsification << "%" << endl;

	outputFileStream.close();
	return 0;
}