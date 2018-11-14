#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <bf/all.hpp>
#include "FastaParser.h"

using namespace std;

int main (int argc, char *argv[]) 
{
	if(argc < 3) {
		cerr << "Please write two arguments: path to fasta file an number of k-mers"
			<< endl;
		return -1;
	}
	cerr << "Classic Bloom Filter" << endl;

	int K = atoi(argv[2]);
	string fastaFile = argv[1];

	FastaParser fp(fastaFile, K);
	fp.parseKmers();
	fp.printKmers();
}