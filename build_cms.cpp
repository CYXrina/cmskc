#include <zlib.h>
#include <stdio.h>
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
	fprintf(stderr, "k\tk-mer length\n");
	fprintf(stderr, "h\tnumber of hashes\n");
	fprintf(stderr, "d\tnumber of cells for each hash\n");
	fprintf(stderr, "t\tlength of heading zeros for a hash to be considered (ntCard sampling)\n");
	fprintf(stderr, "o\toutput file for the saved count-min sketch\n");
	fprintf(stderr, "i\tinput (gzipped) file. stdin as default if not specified\n");
	fprintf(stderr, "\nExample usage: \n");
	fprintf(stderr, "cat file.gz | cmskc -k 15 -h 7 -d 1000 -o sketch.cms\n");
	fprintf(stderr, "cmskc -i file.gz -k 15 -h 7 -d 1000 -o sketch.cms\n");
}

int main(int argc, char** argv)
{
	uint64_t k, h, d, t, seq_len, mask;
	int c, streaming;
	FILE *instream, *outstream;
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
	t = 0;
	instream = nullptr;
	outstream = nullptr;
	while((c = ketopt(&opt, argc, argv, 1, "k:h:d:t:o:i:", longopts)) >= 0) {
		//fprintf(stderr, "opt = %c | arg = %s\n", c, opt.arg);
		if (c == 'k') k = strtoull(opt.arg, nullptr, 10);
		else if (c == 'h') h = strtoull(opt.arg, nullptr, 10);
		else if (c == 'd') d = strtoull(opt.arg, nullptr, 10);
		else if (c == 't') t = strtoull(opt.arg, nullptr, 10);
		else if (c == 'o') {
			if(!(outstream = fopen(opt.arg, "w+b"))) {
				fprintf(stderr, "Unable to open output file\n");
				return -2;
			//output_file = static_cast<char*>(malloc(strlen(opt.arg) + 1));
			//strcpy(output_file, opt.arg);
			}
		} else if (c == 'i') {
			if(!(instream = fopen(opt.arg, "r"))) {
				fprintf(stderr, "Unable to open input file\n");
				return -2;
			}
		} else {
			fprintf(stderr, "Option (%c) not available\n", c);
			print_help();
			return -1;
		}
	}

	//count-min set-up
	cms_init(&cms, d, h);
	if(cms.width != d || cms.depth != h) {
		fprintf(stderr, "Error initializing the sketch\n");
		return -2;
	}

	//Check if stdin is used
	if(instream == nullptr) {
		instream = stdin;
	}
	fp = gzdopen(fileno(instream), "r");
	seq = kseq_init(fp);
	
	uint64_t *hVec = static_cast<uint64_t*>(malloc(sizeof(uint64_t) * h));
	if(t > sizeof(uint64_t) * 8) t = sizeof(uint64_t) * 8;
	mask = generate_mask(t);
	while(kseq_read(seq) >= 0)
	{
		seq_len = seq->seq.l;
		//fprintf(stderr, "Sequence length = %lu\n", seq->seq.l);
		//printf("%s\n%s\n", seq->seq.name, seq->seq.s);
		bool add_to_sketch = false;
		if(seq_len >= k) {
			first = true;
			for(size_t i = 0; i < seq_len - k; ++i)
			{
				if(first) //find the first good kmer
				{
					//fprintf(stderr, "Finding good k-mer\n");
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
					//fprintf(stderr, "Good k-mer found\n");
					if(i < seq_len - k)
					{
						//fprintf(stderr, "Ok, add it to the sketch\n");
						NTM64(&seq->seq.s[i], k, h, hVec);
						if((hVec[0] & mask) == 0) add_to_sketch = true;
						else add_to_sketch = false;
					
					} else {
						//fprintf(stderr, "The length of the k-mer is not enough\n");
						add_to_sketch = false;
					}
				} else {
					//fprintf(stderr, "Rolling %lu: \n", i);
					if(seedTab[seq->seq.s[i+k-1]] != 0)
					{
						//fprintf(stderr, "ok, the char is good\n");
						NTM64(seq->seq.s[i-1], seq->seq.s[i+k-1], k, h, hVec);
						if((hVec[0] & mask) == 0) add_to_sketch = true;
						else add_to_sketch = false;
					
					} else {
						//fprintf(stderr, "INTERRUPT\n");
						i = i+k-1;
						first = true;
						add_to_sketch = false;
					}
				}
				if(add_to_sketch) //ok, add the sketches to the vector
				{
					fprintf(stderr, "Adding hash to sketch\n");
					//fprintf(stdout, "[");
					//for(int i = 0; i < h; ++i) fprintf(stdout, "%lu, ", hVec[i]);
					//fprintf(stdout, "]\n");
					int32_t freq = cms_add_alt(&cms, hVec, h);
					//fprintf(stderr, "c[i+k] = %c, c[i] = %c at i = %lu\n", seq->seq.s[i+k], seq->seq.s[i], i);
					//fprintf(stderr, "Element added, current estimated frequency = %i\n", freq);
				}
			}
			} else {
			fprintf(stderr, "SEQUENCE TOO SHORT\n");
		}
	}
	free(hVec);
	kseq_destroy(seq);
	gzclose(fp);
	fclose(instream);

	//cms_print(&cms);

	cms_write_to_file(&cms, outstream, 2, k, t);
	fclose(outstream);
}
