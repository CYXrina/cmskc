#include <zlib.h>
#include <stdio.h>
#include <limits>
#include <unordered_map>
#include "ketopt.h"
#include "kseq.h"
#include "nthash.hpp"
extern "C" {
#include "count_min_sketch.h"
}

KSEQ_INIT(gzFile, gzread)

using namespace std;

uint64_t generate_mask(size_t head_len)
{
	if(head_len == 0) return 0;
	return static_cast<uint64_t>(~0) << (sizeof(uint64_t) * 8 - head_len);
}

void print_help()
{
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "h\thash row index\n");
	fprintf(stderr, "d\thash column index\n");
	fprintf(stderr, "c\tcms sketch\n");
	fprintf(stderr, "i\tinput (gzipped) file. stdin as default if not specified\n");
	fprintf(stderr, "Output to stdout\n");
}

int main(int argc, char** argv)
{
	uint64_t k, h, d, t, seq_len, mask, *buf;
	int c, streaming;
	char *input_file;
	FILE *instream;
	bool first;
	gzFile fp;
	kseq_t *seq;
	ketopt_t opt = KETOPT_INIT;
	CountMinSketch cms;

	char *debug_kmer;
	
	static ko_longopt_t longopts[] = {
		{NULL, 0, 0}
	};

	if(argc < 9) {
		fprintf(stderr, "Not enough arguments\n");
		print_help();
		if(argc == 1) return 0;
		return -1;
	}

	//fprintf(stderr, "Initialization and Checking done\n");

	//reading command-line options
	input_file = nullptr;
	while((c = ketopt(&opt, argc, argv, 1, "h:d:c:i:", longopts)) >= 0) {
		//fprintf(stderr, "opt = %c | arg = %s\n", c, opt.arg);
		if(c == 'h') h = strtoull(opt.arg, nullptr, 10);
		else if(c == 'd') d = strtoull(opt.arg, nullptr, 10);
		else if(c == 'c') {
			instream = fopen(opt.arg, "r");
			if(cms_read_from_file(instream, &cms, &buf) != 2) {
				fprintf(stderr, "Unable to read k and t\n");
				free(buf);
				fclose(instream);
				return -2;
			}
			fclose(instream);
			k = buf[0];
			t = buf[1];
			free(buf);
		}
		else if(c == 'i') {
			input_file = new char[strlen(opt.arg)];
			strcpy(input_file, opt.arg);
		}
		else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_help();
			return -1;
		}
	}
	
	if(k > numeric_limits<int>::max()) {
		fprintf(stderr, "Error, I cannot print k-mers longer than (2^31 - 1)\n");
		return -1;
	}

	instream = nullptr;
	if(input_file) {
		instream = fopen(input_file, "r");
	}
	else instream = stdin;
	if(!instream) {
		fprintf(stderr, "Unable to open the input file\n");
		return 1;
	}
	fp = gzdopen(fileno(instream), "r");
	seq = kseq_init(fp);
	
	unordered_map<string, size_t> kmer_bag;
	uint64_t hVec[h];
	debug_kmer = new char[k];
	if(t > sizeof(uint64_t) * 8) t = sizeof(uint64_t) * 8;
	mask = generate_mask(t);
	//fprintf(stderr, "mask = %lu\n", mask);
	while(kseq_read(seq) >= 0)
	{
		seq_len = seq->seq.l;
		bool add_to_sketch = false;
		if(seq_len >= k) {
			first = true;
			for(size_t i = 0; i < seq_len - k; ++i)
			{
				if(first)
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
				} else {
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
				if(add_to_sketch) //ok, add the sketches to the vector
				{
					if(cms_check_cell_binning_alt(&cms, hVec, h, d)) 
					{
						memcpy(debug_kmer, &seq->seq.s[i], k);
						string buff(debug_kmer);
						if(kmer_bag.count(buff) == 0) kmer_bag[buff] = 1;
						else ++kmer_bag[buff];
						//printf("%.*s", static_cast<int>(k), seq->seq.s + i);
					}	
				}
			}
		} else {
			fprintf(stderr, "SEQUENCE TOO SHORT\n");
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	fclose(instream);

	for(const auto& p : kmer_bag)
	{
		printf("%s %lu\n", p.first.c_str(), p.second);
	}
	return 0;
}
