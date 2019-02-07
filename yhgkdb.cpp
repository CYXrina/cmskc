#include <zlib.h>
#include <stdio.h>
#include <cctype>
#include "ketopt.h"
#include "kseq.h"
#include "nthash.hpp"
#include "fmm_sketch.h"
extern "C" {
#include "count_min_sketch.h"
}

KSEQ_INIT(gzFile, gzread)

using namespace std;

/*
char safe_toupper(char ch)
{
	return static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
}
*/

uint8_t get_length_of_leading_zeros(int32_t max_freq)
{
	const uint32_t mask = 0x80000000;
	uint8_t n = 0;
	while(!(max_freq & mask) && n < 33)
	{
		++n;
		max_freq <<= 1;
	}
	return n;
}

uint64_t generate_mask(size_t tail_len)
{
	return (~0) >> (sizeof(uint64_t) - tail_len);
}

void print_help()
{
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "p\tnumber of bits of the frequencies to be used for bucketing\n");
	fprintf(stderr, "l\tnumber of heading bits of each hash for bucketing inside the HyperMinHash\n");
	fprintf(stderr, "b\tlength of the b-bit signature to keep for the hyperMinHash in the interval [0, 10]\n");
	fprintf(stderr, "s\tinitial width step. The first counter interval goes from 0 to s excluded\n");
	fprintf(stderr, "q\tratio between the width of counter interval [i+1] and the width of counter interval [i]\n");
	fprintf(stderr, "c\tpath to the pre-builded Count-Min-Sketch\n");
	fprintf(stderr, "o\toutput file for the sketch. If not specified stdout is used.\n");
	fprintf(stderr, "i\tinput (gzipped) file. stdin as default if not specified\n");
	//fprintf(stderr, "\nExample usage: \n");
	//fprintf(stderr, "\n");
	//fprintf(stderr, "\n");

}

int main(int argc, char* argv[])
{
	size_t k, t, s, seq_len, sk_len, *pars;
	int c, p, l, b;
	double q;
	gzFile fp;
	kseq_t *seq;
	ketopt_t opt = KETOPT_INIT;
	FILE *instream, *outstream, *cmsstream;
	CountMinSketch cms;
	hllmh *sk_vec;
	int32_t freq_max, freq_min;

	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};
	/*
	if(argc < 9) {
		fprintf(stderr, "Not enough arguments\n");
		print_help();
		if(argc == 1) return 0;
		return -1;
	}
	*/
	//fprintf(stderr, "Initialization and Checking done\n");

	//reading command-line options
	instream = nullptr;
	outstream = nullptr;
	p = 10;
	l = 10;
	b = 5;
	while((c = ketopt(&opt, argc, argv, 1, "p:l:b:s:q:c:o:i:", longopts)) >= 0) {
		//fprintf(stderr, "opt = %c | arg = %s\n", c, opt.arg);
		if(c == 'p') {
			p = atoi(opt.arg);
			if(p < 0) {
				fpritnf(stderr, "p must be positive\n");
				return -3;
			}
		}
		else if(c == 'l') {
			l = atoi(opt.arg);
			if(l < 0 || l > 64) {
				fprintf(stderr, "l cannot be negative or greater than the number of bits of the hashes\n");
			}
		}
		else if(c == 'b') {
			b = atoi(opt.arg);
			if(b > 10 || b < 0) {
				fprintf(stderr, "b must be in the interval [0, 10]\n");
				return -4;
			}
		}
		else if(c == 's') {
			s = strtoull(opt.arg, 10, nullptr);
		}
		else if(c == 'q') {
			q = atof(opt.arg);
			if(q < 1.0) {
				fprintf(stderr, "q must be strictly greater than 1\n");
				return -5;
			}	
		}
		else if(c == 'c') {
			cmsstream = fopen(opt.arg, "r+b");
			if(cmsstream == nullptr) {
				fprintf(stderr, "Unable to open count-min sketch file\n");
				return -2
			}
			if(cms_read_from_file(&cms, opt.arg, pars) != 2) {
				fprintf(stderr, "This count-min sketch has something else other than k and t as additional information\n");
				free(pars);
				return -2;
			}
			fclose(cmsstream);
			k = pars[0];
			t = pars[1];
			free(pars);
		}
		else if(c == 'o') {
			if(!(instream = fopen(opt.arg, "r"))) {
				fprintf(stderr, "Unable to open input file\n");
				return -2;
			}
		}
		else if(c == 'i') {
			if(!(outstream = fopen(opt.arg, "w"))) {
				fprintf(stderr, "Unable to open output file\n");
				return -2;
			}
		}
		else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_help();
			return -1;
		}
	}
	if(instream == nullptr) instream = stdin;
	if(outstream == nullptr) outstream = stdout;
	if(k == 0) {
		fprintf(stderr, "Use -k to specify the k-mer length (or if you insterted 0 it is not a proper setting\n");
	}
	fp = gzdopen(fileno(instream), "r");
	seq = kseq_init(fp);

	
	freq_max = 0;
	freq_min = numeric_limits<int32_t>::max();
	for(size_t i = 0; i < cms.depth * cms.width; ++i)
	{
		if(cms.bins[i] > freq_max) freq_max = cms.bins[i];
		if(cms.bins[i] < freq_min) freq_min = cms.bins[i];
	}

	/*
	size_t zeros_header_len = get_length_of_leading_zeros(freq_max);

	if(32 - zeros_header_len < p) {
		fprintf(stderr, "Warning! p value greater than maximum suffix frequency length\n");
		fprintf(stderr, "Consequence: The sketch degenerates to a simple minHash\n");
	}
	*/
	
	
	//Create the HyperMinHash vector
	sk_len = static_cast<size_t>(ceil(log(static_cast<double>(freq_max)/static_cast<double>(s))/log(q)));
	sk_vec = malloc(sizeof(hllmh), sk_len);
	for(size_t i = 0; i < sk_len; ++i) 
		hllmh_init(&sk_vec[i], static_cast<uint8_t>(l), static_cast<uint8_t>(b), 1); //for now the number of minimums is 1

	//The actual algorithm loop
	uint64_t hVec[cms.depth];
	bool first = true;
	while(kseq_read(seq) >= 0)
	{
		seq_len = seq->seq.l;
		bool add_to_sketch = false;
		if(seq_len >= k) 
		{
			for(size_t i = 0; i < seq_len - k; ++i)
			{
				if(first) //find the first good kmer
				{
					first = false;
					for(;seedTab[seq->seq.s[i]] == 0 && i < seq_len; ++i); //lock i into good position
					size_t j = i;
					size_t len = 0;
					while(len < k && j < seq_len) //check the next k bases if they are good
					{
						if(seedTab[seq->seq.s[j]] != 0)
						{
							++j;
							++len;
						}
						else {
							i = ++j;
							len = 0;
						}
					}
					if(i < seq_len - k)
					{
						NTM64(&seq->seq.s[i], k, h, hVec);
						if((hVec[0] & mask) == 0) add_to_sketch = true;
						else add_to_sketch = false;
					} else {
						add_to_sketch = false;
					}
				} else { //roll
					if(seedTab[seq->seq.s[i+k-1]] != 0)
					{
						NTM64(seq->seq.s[i-1], seq->seq.s[i+k-1], k, h, hVec);
						if((hVec[0] & mask) == 0) add_to_sketch = true;
						else add_to_sketch = false;
					} else {
						i = i+k-1;
						first = true;
						add_to_sketch = false;
					}
				}
				if(add_to_sketch) //ok, process the frequency and the k-mer
				{
					size_t get_idx = [] (size_t s, double q, int32_t f) {
						return static_cast<size_t>(ceil( log(static_cast<double>(f) / static_cast(s)) / log(q) )) - 1;
					};

					int32_t freq = cms_check_alt(&cms, hVec, cms.depth);
					hllmh_add(&sk_vec[get_idx(freq)], hVec[j]);
				}
			}
		} else {
			fprintf(stderr, "SEQUENCE TOO SHORT\n");
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	fclose(instream);
	
	//For the moment the function saves everything, it might be useful for non homogeneous minHashes
	dump_hllmh_vector(sk_vec, sk_len, s, q, outstream);
	fclose(outstream);
}
