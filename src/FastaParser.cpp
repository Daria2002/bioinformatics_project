#include "FastaParser.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

FastaParser::FastaParser(std::string input_file, int K)
	: fasta_stream(input_file), K(K) {
	cout << "FastaParser() - Initialized stream: " << input_file << endl;
}

FastaParser::~FastaParser()
{
}

// parse kmers from fasta file
unordered_set<string> FastaParser::ParseKmers() {
	cout << "FastaParser::Parse() - starting to parse file..." << endl;
	string current_line;
	int counter = 0;
	while (getline(fasta_stream, current_line)) {
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>') {
			cout << "Skipping header: " << current_line << endl;
			continue;
		}

		// ... if not header start to parse line
		if (counter == 0)
			cout << "Parsing fasta file"<< endl;

		const char* kCLine = current_line.c_str();
		
		for (int i = 0; i < current_line.size() - K + 1; i += 1) { 
			char* new_kmer = (char*)calloc(K, sizeof(char));
			strncpy(new_kmer, kCLine + i, K);
			string* new_kmer_string = new string(new_kmer);
			// Try to insert it, if not successful, free the memory
			if (!kmer_set.insert(*new_kmer_string).second) {
				free(new_kmer);
				free(new_kmer_string);
			}
		}

		counter++;
	}
	cout << "FastaParser::parse() - Finished parsing!" << endl;
	return kmer_set;
}

// parse kmers from fasta file and returns kmer set in vector
vector<string> FastaParser::ParseKmersToVector() {
	cout << "FastaParser::Parse() - starting to parse file..." << endl;
	string current_line;
	int counter = 0;
	while (getline(fasta_stream, current_line)) {
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>') {
			cout << "Skipping header: " << current_line << endl;
			continue;
		}

		// ... if not header start to parse line
		if (counter == 0)
			cout << "Parsing fasta file" << endl;

		const char* kCLine = current_line.c_str();
		
		for (int i = 0; i < current_line.size() - K + 1; i += 1) { 
			char* new_kmer = (char*)calloc(K, sizeof(char));
			strncpy(new_kmer, kCLine + i, K);
			string* new_kmer_string = new string(new_kmer);
			// Try to insert it, if not successful, free the memory
			kmer_vector.push_back(*new_kmer_string);
		}

		counter++;
	}
	cout << "FastaParser::parse() - Finished parsing!" << endl;
	return kmer_vector;
}

// save edge kmers in vector
// if kmer doesn't have right and left neighbour saved in bf, 
// check if kmer is in egde kmers
std::vector<string> FastaParser::EdgeKmers()
{
	vector<string> edge_neighbours;
	string current_line;
	while (getline(fasta_stream, current_line)) {
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>')
			continue;

		const char* kCLine = current_line.c_str();

		//Counting how many times each kmer accurate in sequence
		for (int i = 0; i < current_line.size() - K + 1; i++) { 
			char* kmer = (char*)calloc(K, sizeof(char));
			strncpy(kmer, kCLine + i, K);
			string* kmer_string = new string(kmer);
			if (i == current_line.size() - K || i == 0)
				edge_neighbours.push_back(*kmer_string);
		}	
	}
	return edge_neighbours;
}	

// best index sparsification takes kmers from best index
// best index has largest number of kmers as already saved in kmer set
vector<string> FastaParser::BestIndexSparsification()
{
	string current_line;
	int counter = 0;
	int s = 1;

	vector<string> kmers_to_save;
	bool already_in_set;
	while(getline(fasta_stream, current_line)) {
		if (current_line.at(0) == '>')
			continue;
		// save all kmers from first line
		if (counter == 0) {
			const char* kCLine = current_line.c_str();
			for (int i = 0; i <= current_line.size() - K; i = i + 1 + s) { 
				already_in_set = false;
				char* new_kmer = (char*)calloc(K, sizeof(char));
				strncpy(new_kmer, kCLine + i, K);
				string* new_kmer_string = new string(new_kmer);
				//check if kmers are already saved
				for (string kmer : kmers_to_save) {
					if (kmer == *new_kmer_string) {
						already_in_set = true;
						break;
					}
				}
				if (already_in_set==false)
					kmers_to_save.push_back(*new_kmer_string);
			}
		}
		//save kmers from best Kr
		//Kr is kmer set started at index r
		else if (counter > 0) {	
			unordered_map<int, int> best_index_map;
			const char* kCLine = current_line.c_str();
			int max_freq = 0;
			int best_index = 0;
			//take different indexes
			
			for (int start_index = 0; start_index <= K + 1; start_index = start_index + 1) {
				vector<string> kmers_at_best_index;
				for (int i = start_index; i <= current_line.size() - K; i = i + 1 + s) { 
					char* new_kmer = (char*)calloc(K, sizeof(char));
					strncpy(new_kmer, kCLine+i, K);
					string* new_kmer_string = new string(new_kmer);
					kmers_at_best_index.push_back(*new_kmer_string);
					// compare kmer with already saved kmers
					// but only with 50 kmers because it's faster

					int count = 0;
					for (string saved_kmers : kmers_to_save) {
						if (count > 50)
							break;
						if (*new_kmer_string == saved_kmers) {
							best_index_map[start_index]++;
							break;
						}
						count++;
					}
				}
				// compare number of already saved kmers with last maximum of same kmers
				if (best_index_map[start_index] >= max_freq) {
					max_freq = best_index_map[start_index];
					best_index = start_index;
				}
			}
			
			//add kmers from best index in kmers set
			for (int i = best_index; i <= current_line.size() - K; i = i + s + 1) {
				bool already_in_set = false;
				char* new_kmer = (char*)calloc(K, sizeof(char));
				strncpy(new_kmer, kCLine + i, K);
				string* new_kmer_string = new string(new_kmer);
				// it's much faster if duplicates are not checked
				/*
				for (auto kmer : kmers_to_save)
				{
					if (kmer == *new_kmer_string)
					{
						already_in_set = true;
						break;
					}
				}
				if (!already_in_set)*/
					kmers_to_save.push_back(*new_kmer_string);
			}
		}
		counter++;
	}
	return kmers_to_save;
}

//returns neighbours for every kmer in sequence
//key in map is kmer, value in map are neighbours for kmer
unordered_map<string, int> MakeMapWithKmersAndNeighbours(string sequence, int K) {
	unordered_map<string, int> neighbours;
	string left_neighbour;
	string kmer_help;
	const char* kCLine = sequence.c_str();
	//first kmer doesn't have left neighbour
	char* kmer = (char*)calloc(K, sizeof(char));
	strncpy(kmer, kCLine, K);
	string* kmer_string = new string(kmer);
	left_neighbour = *kmer_string;
	//loop that adds left neighbours, starting from second kmer because first kmer doesn't have left neighbour
	for (int i = 1; i < sequence.size() - K; i = i + 1) {
		char* kmer = (char*)calloc(K, sizeof(char));
		strncpy(kmer, kCLine + i, K);
		string* kmer_string = new string(kmer);
		neighbours[*kmer_string]++;
	}
	//provjeri j>0
	//loop that adds right neighbours, starting from last kmers 
	for (int j = sequence.size() - K; j >= 0; j = j - 1) {
		char* rightNeighbour = (char*)calloc(K, sizeof(char));
		//last kmer
		strncpy(rightNeighbour, kCLine + j, K);
		string* kmer_string = new string(rightNeighbour);
		neighbours[*kmer_string]++;
	}
	return neighbours;
}

// delete kmers that appears more then once in neighbour sets of other kmers
vector<string> FastaParser::HittingSetSparsification()
{
	vector<string> kmers_to_save;
	string current_line;
	int count = 0;
	while (getline(fasta_stream, current_line)) {
		count=count+1;
		unordered_map<string, int> kmers_and_neighbours;
		// If first character is '>' then omit the line
		if (current_line.at(0) == '>')
			continue;
		// returns map with kmers and neighbours for each kmers
		// kmer appeared as neighbours are map key and
		// map value is how much each kmer appeared as neighbour
		kmers_and_neighbours = MakeMapWithKmersAndNeighbours(current_line, K);
		const char* kCLine = current_line.c_str();
		bool skip = false;
		vector<string> kmers_to_delete;
		int max_freq = 0;
		//kmers that appeared more than once are saved in kmers_to_delete
		for (auto kmer = kmers_and_neighbours.begin(); kmer != kmers_and_neighbours.end(); kmer++) {
			if (kmer -> second > max_freq)
				max_freq = kmer -> second;
		}

		for (auto kmer = kmers_and_neighbours.begin(); kmer != kmers_and_neighbours.end(); kmer++) {
			if (kmer -> second < max_freq)
				kmers_to_delete.push_back(kmer -> first);
		}
			
		//adding kmers that appeared only once in kmers_to_save vector
		for (int i = 0; i < current_line.size() - K + 1; i++) {
			skip = false;
			char* new_kmer = (char*)calloc(K, sizeof(char));
			strncpy(new_kmer, kCLine + i, K);

			string* new_kmer_string = new string(new_kmer);

			for (string kmer_to_Delete : kmers_to_delete) {
				if (*new_kmer_string == kmer_to_Delete) {
					skip = true;
					break;
				}
			}
			if (skip)
				continue;
			kmers_to_save.push_back(*new_kmer_string);
		}
	}
	return kmers_to_save;
}	

void FastaParser::PrintKmers()
{
	cout << K << "-mers are :" << endl;
	int i = 0, kmer_size = kmer_set.size();
	for (string kmer : kmer_set) {
		i++;
		cout << i << "/" << kmer_size << ":\t" << kmer << endl;
	}
}
