#include <zlib.h>
#include <stdio.h>
#include <cctype>
#include "ketopt.h"
#include "kseq.h"
#include "nthash.hpp"
extern "C" {
#include "count_min_sketch.h"
}

KSEQ_INIT(gzFile, gzread)

using namespace std;

unsigned char seq_table[256] = {
	  0,   1,   2,   3,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255,   0, 255,   1,  255, 255, 255,   2,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,    3,   3, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255,   0, 255,   1,  255, 255, 255,   2,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,    3,   3, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,
	255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255
};

char safe_toupper(char ch)
{
	return static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
}

void print_help()
{
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "k\tk-mer length\n");
	fprintf(stderr, "h\tnumber of hashes\n");
	fprintf(stderr, "d\tnumber of cells for each hash\n");
	fprintf(stderr, "o\toutput file for the saved count-min sketch\n");
	fprintf(stderr, "i\tinput (gzipped) file. stdin as default if not specified\n");
	fprintf(stderr, "\nExample usage: \n");
	fprintf(stderr, "cat file.gz | cmskc -k 15 -h 7 -d 1000 -o sketch.cms\n");
	fprintf(stderr, "cmskc -i file.gz -k 15 -h 7 -d 1000 -o sketch.cms\n");
}

int main(int argc, char** argv)
{
	size_t k, h, d, seq_len;
	int c, streaming;
	char output_file[100];
	char *input_file;
	gzFile fp;
	kseq_t *seq;
	ketopt_t opt = KETOPT_INIT;
	CountMinSketch cms;
	FILE *instream;
	
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
	while((c = ketopt(&opt, argc, argv, 1, "k:h:d:o:i:", longopts)) >= 0) {
		//fprintf(stderr, "opt = %c | arg = %s\n", c, opt.arg);
		if(c == 'k') k = strtoull(opt.arg, nullptr, 10);
		else if(c == 'h') h = strtoull(opt.arg, nullptr, 10);
		else if(c == 'd') d = strtoull(opt.arg, nullptr, 10);
		else if(c == 'o') strcpy(output_file, opt.arg);
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
	
	//fprintf(stderr, "Input parameters read\n");

	//count-min set-up
	cms_init(&cms, d, h);
	if(cms.width != d || cms.depth != h) {
		fprintf(stderr, "Error initializing the sketch\n");
		return -2;
	}
	
	//fprintf(stderr, "Count min setted up\n");

	//Opening the file
	instream = nullptr;
	if(input_file) {
		fprintf(stderr, "%s\n", input_file);
		instream = fopen(input_file, "r");
		fprintf(stderr, "check 0\n");
		delete [] input_file;
	}
	else instream = stdin;
	fprintf(stderr, "check 1\n");
	fp = gzdopen(fileno(instream), "r");
	fprintf(stderr, "check 2\n");
	seq = kseq_init(fp);
	fprintf(stderr, "check 3\n");

	//fprintf(stderr, "Input file opened\n");

	//Input handling
	uint64_t hVec[h];
	bool first = true;
	while(kseq_read(seq) >= 0)
	{
		seq_len = seq->seq.l;
		//fprintf(stderr, "Sequence length = %lu\n", seq->seq.l);
		bool add_to_sketch = false;
		//sketch the sequences
		if(seq_len >= k) {
		for(size_t i = 0; i < seq_len - k; ++i)
		{
			if(first) //find the first good kmer
			{
				//fprintf(stderr, "Initializing first hash\n");
				first = false;
				for(;seq_table[seq->seq.s[i]] == 255 && i < seq_len; ++i); //lock i into good position
				seq->seq.s[i] = safe_toupper(seq->seq.s[i]);
				size_t j = i;
				size_t len = 0;
				while(len < k && j < seq_len) //check the next k bases if they are good
				{
					seq->seq.s[j] = safe_toupper(seq->seq.s[j]);
					if(seq_table[seq->seq.s[j]] != 255)
					{
						++j;
						++len;
					}
					else {
						i = ++j;
						len = 0;
					}
				}
				//fprintf(stderr, "First hash initialized: ");
				if(i < seq_len - k)
				{
					//fprintf(stderr, "Ok, add it to the sketch\n");
					add_to_sketch = true;
					NTM64(&seq->seq.s[i], k, h, hVec);
				} else {
					//fprintf(stderr, "The length of the k-mer is not enough\n");
					add_to_sketch = false;
				}
			} else {
				//fprintf(stderr, "Rolling %lu: \n", i);
				if(seq_table[seq->seq.s[i+k]] != 255)
				{
					//fprintf(stderr, "ok, the char is good\n");
					NTM64(seq->seq.s[i], safe_toupper(seq->seq.s[i+k]), k, h, hVec);
					add_to_sketch = true;
				} else {
					//fprintf(stderr, "INTERRUPT\n");
					first = true;
					add_to_sketch = false;
				}
			}
			if(add_to_sketch) //ok, add the sketches to the vector
			{
				//fprintf(stderr, "Adding hash to sketch\n");
				int32_t freq = cms_add_alt(&cms, hVec, h);
				//fprintf(stderr, "Element added, current estimated frequency = %i\n", freq);
			}
		}
		} else {
			fprintf(stderr, "SEQUENCE TOO SHORT\n");
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

	//fprintf(stderr, "Exporting count-min\n");
	cms_export(&cms, output_file);
	//fprintf(stderr, "Count-Min exported\n");
}
