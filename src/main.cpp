#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <ctime>
#include <ctime>
#include <bf/all.hpp>
#include <thread>
#include "FastaParser.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <stdio.h>  /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/ 
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#include <iostream>
#include <libgen.h>
#endif
#include<iostream>
using namespace std;

char RandomKmerChar() {
	int random_acgt = rand() % 4;
	switch (random_acgt) {
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
vector<string> SingleSequenceSparsification(vector<string> kmer_set, int s) {
	vector<string> kmer_set_sparse;
	int i = 0;
	for (auto kmer : kmer_set) {
		if (i % (s + 1) == 0)
			kmer_set_sparse.push_back(kmer);
		i++;
	}
	return kmer_set_sparse; 
}

bool CheckEdgeKmer(vector<string> edge_kmers_set, string kmer) {
	if (std::find(edge_kmers_set.begin(), edge_kmers_set.end(), kmer) != edge_kmers_set.end())
		return true;
	return false;
}

bool DecidePresent(
	string kmer, 
	bool contains_left, 
	bool contains_right, 
	vector<string> edge_kmers_set) {
	if (contains_left || contains_right) {
		if (contains_right && contains_left)
			return true;
		else
			return CheckEdgeKmer(edge_kmers_set, kmer);
	}
	return false;
}

bool ContainsSet(vector<string> kmer_set, const bf::basic_bloom_filter &kBloomFilter) {
	for (string kmer : kmer_set)
		if (kBloomFilter.lookup(kmer))
			return true;
	return false;
}

string ConvertToBase4(int decimal_number, int length) {
	string number_in_base4_inverse;
	string number_in_base4;
	int remainder;
	for (int i = 0; i < length; i++) {
		while(decimal_number > 0) {
			remainder = decimal_number % 4;
			decimal_number = decimal_number / 4;
			number_in_base4_inverse = number_in_base4_inverse + to_string(remainder);
		}
		for (int j = number_in_base4_inverse.size(); j < length; j++)
			number_in_base4_inverse = number_in_base4_inverse + '0';
	}
	for (int i = number_in_base4_inverse.size() - 1; i >= 0; i--)
		number_in_base4 += number_in_base4_inverse.at(i);
	return number_in_base4;
}

vector<string> PrefixOrSuffixCombinations(int length) {
	int number_of_combinations = pow(4, length);
	vector<string> combinations;
	string number_base4;
	for (int i = 0; i < number_of_combinations; i++) {
		number_base4 = ConvertToBase4(i, length);
		string combination;
		for (auto el : number_base4) {
			if (el == '0')
				combination += 'A';
			else if (el == '1')
				combination += 'C';
			else if (el == '2')
				combination += 'T';
			else if (el == '3')
				combination += 'G';
		}
		combinations.push_back(combination);
	}
	return combinations;
}

// this function is used only for relaxed method
vector<string> SDistantLeftNeighbourSet(string kmer, int left_dist) {
	vector<string> neighbours;
	for (int i = 1; i <= left_dist; i++) {
		// make prefix_set and base
		// prefix set are all combinations of prefixes 
		// base consists of first k-i element in kmer
		vector<string> prefix_set = PrefixOrSuffixCombinations(i);
		string base = kmer.substr(0, kmer.size() - i);
		for (auto prefix : prefix_set)
			neighbours.push_back(prefix + base);
	}
	return neighbours;
}

// this function is used only for relaxed mothod
vector<string> SDistantRightNeighbourSet(string kmer, int right_dist) {
	vector<string> neighbours;
	for (int i = 1; i <= right_dist; i++) {
		vector<string> suffix_set = PrefixOrSuffixCombinations(i);
		string base = kmer.substr(i, kmer.size() - i);
		for (auto suffix:suffix_set)
			neighbours.push_back(base + suffix);
	}
	return neighbours;
}

bool StrictContainsNeighbours(
	string kmer,
	int left_dist, 
	int right_dist, 
	const bf::basic_bloom_filter &kBloomFilter, 
	vector<string> edge_kmers_set) {
	vector<string> left_neighbours;
	vector<string> right_neighbours;
	vector<string> prefix_set = PrefixOrSuffixCombinations(left_dist);
	vector<string> suffix_set = PrefixOrSuffixCombinations(right_dist);
	string base_left = kmer.substr(0, kmer.size() - left_dist);
	string base_right = kmer.substr(right_dist, kmer.size() - right_dist);
	for (auto prefix : prefix_set)
		left_neighbours.push_back(prefix + base_left);
	for (auto suffix : suffix_set)
		right_neighbours.push_back(base_right + suffix);
	return DecidePresent(kmer,
		ContainsSet(left_neighbours, kBloomFilter),
		ContainsSet(right_neighbours, kBloomFilter),
		edge_kmers_set);
}

bool RelaxedContainsNeighbours(
	string kmer,
	int left_dist, 
	int right_dist, 
	const bf::basic_bloom_filter &kBloomFilter, 
	vector<string> edge_kmers_set) {
	return DecidePresent(kmer,
		ContainsSet(SDistantLeftNeighbourSet(kmer, left_dist), kBloomFilter),
		ContainsSet(SDistantRightNeighbourSet(kmer, right_dist), kBloomFilter),
		edge_kmers_set);
}

vector<bool> StrictContains(vector<string> kmer_test_set,
 vector<string> edge_kmers_set, const bf::basic_bloom_filter &kBloomFilter, int s) {
	vector<bool> strict_results;
	bool kmer_saved = false;
	for (auto kmer : kmer_test_set) {
		kmer_saved = false;
		// if kmer is saved, search if kmer is in edge_kmers_set and check if neighbour/s are saved 
		if (kBloomFilter.lookup(kmer)) {
			strict_results.push_back(StrictContainsNeighbours(kmer, s + 1, s + 1, kBloomFilter, edge_kmers_set));
		}
		else {
			for (int i = 0; i <= s; i++) {
				if (StrictContainsNeighbours(kmer, i, s - i + 1, kBloomFilter,edge_kmers_set)) {
					strict_results.push_back(true);
					kmer_saved = true;
					break;
				}
			}
			if (!kmer_saved)
				strict_results.push_back(false);
		}
	}
	return strict_results;
}

vector<bool> RelaxedContains(
	vector<string> kmer_test_set, 
	vector<string> edge_kmers_set, 
	const bf::basic_bloom_filter &kBloomFilter, 
	int s) {
	vector<bool> relaxed_results;
	bool kmer_saved = false;
	for (auto kmer : kmer_test_set) {
		kmer_saved = false;
		// if kmer is saved, search if kmer is in edge_kmers_set and check if neighbour/s are saved 
		if (kBloomFilter.lookup(kmer)) {
			relaxed_results.push_back(RelaxedContainsNeighbours(kmer, s + 1, s + 1, kBloomFilter, edge_kmers_set));
		}
		else {
			for (int i = 0; i <= s; i++) {
				if (RelaxedContainsNeighbours(kmer, i, s - i + 1, kBloomFilter, edge_kmers_set)) {
					relaxed_results.push_back(true);
					kmer_saved = true;
					break;
				}
			}
			if (!kmer_saved)
				relaxed_results.push_back(false);
		}
	}
	return relaxed_results;
}

// Returns all combinations for right neighbours
// e.g. if kmer is ACTG, right neighbours can be:
// CTGX, where x={A,G,C,T}  
vector<string> MakeRightNeighbours(string kmer) {
	string char_combinations = "ACTG";
	vector<string> right_neighbours;
	string make_neighbour;
	for (char possibleChar : char_combinations) {	
		make_neighbour = kmer[1];
		for (int i = 2; i < kmer.size(); i++)
			make_neighbour += kmer[i];
		make_neighbour += possibleChar;
		right_neighbours.push_back(make_neighbour);
	}
	return right_neighbours;
}

// Returns all combinations for left neighbour
vector<string> MakeLeftNeighbours(string kmer) {
	string char_combinations = "ACTG";
	vector<string> left_neighbours;
	string make_neighbour;
	for (char possibleChar : char_combinations) {	
		make_neighbour = possibleChar;
		for (int i = 0; i < kmer.size() - 1; i++)
			make_neighbour += kmer[i];
		left_neighbours.push_back(make_neighbour);
	}
	return left_neighbours;
}

// returns test set
vector<string> GenerateTestKmerSet(unordered_set<string> kmer_set, int testSize, int K) {
	vector<string> kmer_set_vecor;
	int index = 0;
	for (string kmer : kmer_set) {
		kmer_set_vecor.push_back(kmer);
		index++;
	}
	int random_kmer_set_index, random_kmer_index;
	vector<string> test_kmers;
	cout << "kmerSetSize: " << kmer_set.size() << endl;
	cout << "K: " << K << endl;
	for (int i = 0; i < testSize; i++) {
		random_kmer_set_index = rand() % (kmer_set.size());
		random_kmer_index = rand() % (K);
		string random_selected_kmer = kmer_set_vecor[random_kmer_set_index];
		char random_char;
		while(1) {
			random_char = RandomKmerChar();
			if (random_char != random_selected_kmer[random_kmer_index])
				break;
		}
		random_selected_kmer[random_kmer_index] = random_char;
		test_kmers.push_back(random_selected_kmer);
	}
	cout << "Finished generating test kmers" << endl;
	return test_kmers;
}

// calculate false positive rate
float FalsePositiveRate(vector<bool> bloom_filter_result, vector<bool> bloom_filter_result_real) {
	float fp_rate = 0;
	for (int i = 0; i < bloom_filter_result.size(); i++) {
		if (!(bloom_filter_result[i] == bloom_filter_result_real[i])) {
			fp_rate++;
		}
	}
	fp_rate = (fp_rate / bloom_filter_result.size())*100 ;
	return fp_rate;
}

static std::string timePointAsString(const std::chrono::system_clock::time_point& tp) {
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::string ts = std::ctime(&t);
    ts.resize(ts.size()-1);
    return ts;
}

//check if kmers in kmer_set_test are in kmer_set
//kmer_set is set of all kmers that really exists in fasta file
//kmer_set is seto of test kmers
//function returns bool vector, where false means that kmer doesn't exist in fasta file and true that it exists
vector<bool> CompareTestKmerWithSavedKmers(unordered_set<string> kmer_set, vector<string> kmer_set_test) {
	vector<bool> bloom_filter_result_real;
	bool setted = false;
	for (auto kmer : kmer_set_test) {
		if (kmer_set.find(kmer) != kmer_set.end())
			bloom_filter_result_real.push_back((bool)true);
		else
			bloom_filter_result_real.push_back((bool)false);
	}
	return bloom_filter_result_real;
}

int main (int argc, char *argv[]) {
	srand(time(NULL));
	size_t num_of_cells;
	int num_of_hashes = 2;
	int test_set_size = 100000;
	unordered_set<string> kmer_set;
	vector<string> kmer_set_test;
	if (argc < 3 || argc > 4) {
		cerr << "Please write three arguments: path to fasta file, k-mer length and optionally output file"
			<< endl;
		return -1;
	}
	int K = atoi(argv[2]);
	string fasta_file = argv[1];
	string output_file;
	if (argc == 4) {
		output_file = argv[3];
	}
	else if (argc == 3) {
		char buff[FILENAME_MAX];
	  	GetCurrentDir( buff, FILENAME_MAX );
	  	string current_working_dir(buff);
	  	string results_dir = current_working_dir + "/../results";
	  	string output_file_name = "results.txt";
	  	output_file = results_dir + "/" + output_file_name;
	}

	ofstream output_file_stream;
	output_file_stream.open(output_file);

	cout << "**********************************" << endl;
	cout << "*******Classic Bloom filter*******" << endl;
	cout << "**********************************" << endl;
	std::chrono::duration<double> duration;
	auto memory = 0;
	FastaParser fp(fasta_file, K);
	kmer_set = fp.ParseKmers();
	int number_of_kmers = kmer_set.size();
	num_of_cells = 1024 * 1024 * 32;
	bf::basic_bloom_filter *kBloomFilter;
	kBloomFilter = new bf::basic_bloom_filter(bf::make_hasher(num_of_hashes), num_of_cells);
	for (auto kmer : kmer_set) {
		memory += kmer.size() * sizeof(char);
		kBloomFilter -> add(kmer);
	}
	cout << "Kmers added in bloom filter" << endl;
	size_t result;
	vector<bool> bloom_filter_result;
	vector<bool> bloom_filter_result_real;
	cout << "Starting generating test kmers." << endl;
	kmer_set_test = GenerateTestKmerSet(kmer_set, test_set_size, K);
	auto start_s = chrono::system_clock::now();
	for (auto kmer : kmer_set_test)
		bloom_filter_result.push_back((bool)kBloomFilter -> lookup(kmer));
	cout << "Test kmers checked" << endl;
	auto stop_s = chrono::system_clock::now();
	duration = stop_s - start_s;
	auto time_classic_bloom_filter = duration.count();
	bloom_filter_result_real = CompareTestKmerWithSavedKmers(kmer_set, kmer_set_test);
	cout << "Test results compared with real" << endl;
	float fp_rate_classic_bloom_filter;
	fp_rate_classic_bloom_filter = FalsePositiveRate(bloom_filter_result, bloom_filter_result_real);
	cout << "Size of classic Bloom filter: " << memory << " Bytes" << endl;
	cout << "Classic Bloom Filter-time: " << time_classic_bloom_filter << " s" << endl;
	cout << "Classic Bloom filter-fp rate: " << fp_rate_classic_bloom_filter << " %" << endl;
	output_file_stream << "Size of classic Bloom filter: " << memory << " Bytes" << endl;
	output_file_stream << "Classic Bloom Filter-time: " << time_classic_bloom_filter << " s" << endl;
	output_file_stream << "Classic Bloom filter-fp rate: " << fp_rate_classic_bloom_filter << " %" << endl;

	cout << "******************************************" << endl;
	cout << "*******One-sided k-mer Bloom filter*******" << endl;
	cout << "******************************************" << endl;
	vector<bool> one_sided_result;
	string *kmer_to_change;
	vector<string> left_neighbours;
	vector<string> right_neighbours;
	bool neighbour_in_set;
	start_s = chrono::system_clock::now();
	for (int i = 0; i < bloom_filter_result.size(); i++) {
		if (bloom_filter_result[i]) {
			kmer_to_change = &kmer_set_test[i];
			left_neighbours = MakeLeftNeighbours(*kmer_to_change);
			right_neighbours = MakeRightNeighbours(*kmer_to_change);
			
			for (int j = 0; j < 4; j++) {
				if (kBloomFilter -> lookup(left_neighbours[j]) || kBloomFilter -> lookup(right_neighbours[j])) {
					one_sided_result.push_back(true);
					neighbour_in_set = true;
					break;
				}
			}
			if (!neighbour_in_set)
				one_sided_result.push_back(false);
			else
				neighbour_in_set = false;
		}
		else
			one_sided_result.push_back(false);	
	}
	stop_s = chrono::system_clock::now();
	float fp_rate_one_sided;
	fp_rate_one_sided = FalsePositiveRate(one_sided_result, bloom_filter_result_real);
	duration = stop_s - start_s;
	auto time_one_sided = duration.count() + time_classic_bloom_filter;
	cout << "Size of one-sided Bloom filter: " << memory << " Bytes" << endl;
	cout << "One-sided Bloom Filter-time: " << time_one_sided << " s" << endl;
	cout << "One-sided Bloom filter-fp rate: " << fp_rate_one_sided << " %" << endl;
	output_file_stream << "Size of one-sided Bloom filter: " << memory << " Bytes" << endl;
	output_file_stream << "One-sided Bloom Filter-time: " << time_one_sided << " s" << endl;
	output_file_stream << "One-sided Bloom filter-fp rate: " << fp_rate_one_sided << " %" << endl;

	cout << "******************************************" << endl;
	cout << "*******Two-sided k-mer Bloom filter*******" << endl;
	cout << "******************************************" << endl;
	// Edge kmers
	vector<string> edge_kmers_set;
	FastaParser fasta_parser_edge_kmers(fasta_file, K);
	edge_kmers_set = fasta_parser_edge_kmers.EdgeKmers();
	vector<bool> two_sided_result;
	neighbour_in_set = false;
	start_s = chrono::system_clock::now();
	for (int i = 0; i < bloom_filter_result.size(); i++) {
		if (bloom_filter_result[i]) {
			kmer_to_change = &kmer_set_test[i];
			left_neighbours = MakeLeftNeighbours(*kmer_to_change);
			right_neighbours = MakeRightNeighbours(*kmer_to_change);
			for (auto left_neighbour : left_neighbours) {
				for (auto right_neighbour : right_neighbours) {
					if (kBloomFilter -> lookup(right_neighbour) && kBloomFilter -> lookup(left_neighbour)) {
						two_sided_result.push_back(true);
						neighbour_in_set = true;
						break;	
					}
					// if only one neighbour is in set, check if kmer is in edge_kmers_set
					// kmers in edge_kmers_set should have only one neighbour saved
					else if (kBloomFilter -> lookup(right_neighbour) || kBloomFilter -> lookup(left_neighbour)) {
						for (string edge_kmer : edge_kmers_set) {
							if (*kmer_to_change == edge_kmer) {
								two_sided_result.push_back(true);
								neighbour_in_set = true;
								break;
							}
						}
						if (neighbour_in_set)
							break;
					}
				}
				if (neighbour_in_set)
					break;
			}
			if (!neighbour_in_set)
				two_sided_result.push_back(false);
			else
				neighbour_in_set = false;
		}
		else
			two_sided_result.push_back(false);	
	}
	stop_s = chrono::system_clock::now();
	duration = stop_s - start_s;
	auto time_two_sided = duration.count();
	float fp_rate_two_sided;
	fp_rate_two_sided = FalsePositiveRate(two_sided_result, bloom_filter_result_real);
	int s = 1;
	cout << "Size of two-sided Bloom filter: " << memory << " Bytes" << endl;
	cout << "Two-sided Bloom Filter-time: " << time_two_sided << " s" << endl;
	cout << "Two-sided Bloom filter-fp rate: " << fp_rate_two_sided << " %" << endl;
	output_file_stream << "Size of two-sided Bloom filter: " << memory << " Bytes" << endl;
	output_file_stream << "Two-sided Bloom Filter-time: " << time_two_sided << " s" << endl;
	output_file_stream << "Two-sided Bloom filter-fp rate: " << fp_rate_two_sided << " %" << endl;

	cout << "**************************************************************" << endl;
	cout << "***Sparse: Best index match per read sparsification, strict***" << endl;
	cout << "**************************************************************" << endl;
	vector<bool> best_index_set_result;
	//test set doesn't change
	vector<string> best_index_set;
	FastaParser fp_index_set(fasta_file, K);
	memory = 0;
	best_index_set = fp_index_set.BestIndexSparsification();
	bf::basic_bloom_filter *bloom_filter_best_index;
	bloom_filter_best_index = new bf::basic_bloom_filter(bf::make_hasher(num_of_hashes), num_of_cells);
	for (auto kmer : best_index_set) {
		memory += kmer.size() * sizeof(char);
		bloom_filter_best_index -> add(kmer);
	}
	vector<bool> best_index_relaxed_results;
	vector<bool> best_index_strict_results;
	start_s = chrono::system_clock::now();
	best_index_strict_results = StrictContains(kmer_set_test, edge_kmers_set, *bloom_filter_best_index, s);
	stop_s = chrono::system_clock::now();
	duration = stop_s - start_s;
	auto time_sparse = duration.count();
	// kmer_set contains kmers parsed from fasta file
	//checking which kmers from kmer_set_test are in kmer_set
	bloom_filter_result_real = CompareTestKmerWithSavedKmers(kmer_set, kmer_set_test);
	float fp_rate_best_index_relaxed;
	float fp_rate_best_index_strict;
	fp_rate_best_index_relaxed = FalsePositiveRate(best_index_relaxed_results, bloom_filter_result_real);
	fp_rate_best_index_strict = FalsePositiveRate(best_index_strict_results, bloom_filter_result_real);
	cout << "Size of Bloom filter using best index match: " << memory << " Bytes" << endl;
	cout << "Best index Bloom Filter-time: " << time_sparse << " s" << endl;
	cout << "Best index strict Bloom filter-fp rate: " << fp_rate_best_index_strict << " %" << endl;
	output_file_stream << "Size of Bloom filter using best index match: " << memory << " Bytes" << endl;
	output_file_stream << "Best index Bloom Filter-time: " << time_sparse << " s" << endl;
	output_file_stream << "Best index strict Bloom filter-fp rate: " << fp_rate_best_index_strict << " %" << endl;

	cout << "***************************************************************" << endl;
	cout << "***Sparse:Spasification via approximate hitting set, Relaxed***" << endl;
	cout << "***************************************************************" << endl;
	//test set doesn't change
	vector<string> hitting_set;
	memory = 0;
	FastaParser fp_hitting_set(fasta_file, K);
	hitting_set = fp_hitting_set.HittingSetSparsification();
	bf::basic_bloom_filter *bloom_filter_hitting_set;
	bloom_filter_hitting_set = new bf::basic_bloom_filter(bf::make_hasher(num_of_hashes), num_of_cells);
	for (auto kmer : hitting_set) {
		memory += kmer.size() * sizeof(char);
		bloom_filter_hitting_set -> add(kmer);
	}
	vector<bool> hitting_set_selaxed_results;
	vector<bool> hitting_set_strict_results;
	start_s = chrono::system_clock::now();
	hitting_set_selaxed_results = RelaxedContains(kmer_set_test, edge_kmers_set, *bloom_filter_hitting_set, s);
	stop_s = chrono::system_clock::now();
	duration = stop_s - start_s;
	float timeHittingSet = duration.count();
	cout << "Hitting set Bloom Filter-time: " << timeHittingSet << " s" << endl;
	bloom_filter_result_real = CompareTestKmerWithSavedKmers(kmer_set, kmer_set_test);
	float fp_rate_hitting_set_relaxed;
	float fp_rate_hitting_set_strict;
	fp_rate_hitting_set_relaxed = FalsePositiveRate(hitting_set_selaxed_results, bloom_filter_result_real);
	cout << "Size of Bloom filter using hitting set: " << memory << " Bytes" << endl;
	cout << "Hitting set relaxed Bloom filter-fp rate: " << fp_rate_hitting_set_relaxed << " %" << endl;
	output_file_stream << "Size of Bloom filter using hitting set: " << memory << " Bytes" << endl;
	output_file_stream << "Hitting set Bloom Filter-time: " << timeHittingSet << " s" << endl;
	output_file_stream << "Hitting set relaxed Bloom filter-fp rate: " << fp_rate_hitting_set_relaxed << " %" << endl;

	cout << "***************************************************************" << endl;
	cout << "********Sparse: Single Sequence Sparsification, Relaxed********" << endl;
	cout << "***************************************************************" << endl;
	vector<bool> single_sequence_sparsification_set_relaxed_result;
	vector<bool> single_sequence_sparsification_set_strict_result;
	//test set doesn't change
	vector<string> single_sequence_sparsification_set;
	FastaParser fp_sequence_sparsification_set(fasta_file, K);
	vector<string> kmer_set_vecor;
	kmer_set_vecor = fp_sequence_sparsification_set.ParseKmersToVector();
	single_sequence_sparsification_set = SingleSequenceSparsification(kmer_set_vecor, s);
	bf::basic_bloom_filter *bloom_filter_sequence_sparsification_set;
	bloom_filter_sequence_sparsification_set = new bf::basic_bloom_filter(bf::make_hasher(num_of_hashes), num_of_cells);
	memory = 0;
	for (auto kmer : single_sequence_sparsification_set) {
		memory += kmer.size() * sizeof(char);
		bloom_filter_sequence_sparsification_set -> add(kmer);
	}
	start_s = chrono::system_clock::now();
	single_sequence_sparsification_set_relaxed_result = RelaxedContains(
		kmer_set_test, edge_kmers_set, *bloom_filter_sequence_sparsification_set, s);
	stop_s = chrono::system_clock::now();
	single_sequence_sparsification_set_strict_result = StrictContains(
		kmer_set_test, edge_kmers_set, *bloom_filter_sequence_sparsification_set, s);
	duration = stop_s - start_s;
	auto time_sequence_sparsification = duration.count();
	float fp_rate_sequence_sparsification_relaxed;
	fp_rate_sequence_sparsification_relaxed = FalsePositiveRate(
		single_sequence_sparsification_set_relaxed_result, bloom_filter_result_real);
	float fp_rate_sequence_sparsification_strict;
	fp_rate_sequence_sparsification_strict = FalsePositiveRate(
		single_sequence_sparsification_set_strict_result, bloom_filter_result_real);
	cout << "Size of Bloom filter using single sequence sparsification: " << memory << " Bytes" << endl;
	cout << "Sequence Sparsification Bloom Filter-time: " << time_sequence_sparsification << " s" << endl;
	cout << "Sequence sparsification relaxed Bloom filter-fp rate: " << fp_rate_sequence_sparsification_relaxed << " %" << endl;
	cout << "Sequence sparsification strict Bloom filter-fp rate: " << fp_rate_sequence_sparsification_strict << " %" << endl;
	output_file_stream << "Size of Bloom filter using single sequence sparsification: " << memory << " Bytes" << endl;
	output_file_stream << "Sequence Sparsification Bloom Filter-time: " << time_sequence_sparsification << " s" << endl;
	output_file_stream << "Sequence sparsification relaxed Bloom filter-fp rate: " << fp_rate_sequence_sparsification_relaxed << " %" << endl;
	output_file_stream << "Sequence sparsification strict Bloom filter-fp rate: " << fp_rate_sequence_sparsification_strict << " %" << endl;
	output_file_stream.close();
	return 0;
}