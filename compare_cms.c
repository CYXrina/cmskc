#include <stdio.h>
#include <limits.h>
#include "count_min_sketch.h"

void denoise(CountMinSketch *cms)
{
    int32_t min = INT_MAX;
    for(long i = 0; i < cms->depth * cms->width; ++i) if(min > cms->bins[i]) min = cms->bins[i];
	for(long i = 0; i < cms->depth * cms->width; ++i) cms->bins[i] -= min;
}

int main(int argc, char** argv)
{
	double num, den;
	CountMinSketch cms1, cms2;
	cms_import(&cms1, argv[1]);
	cms_import(&cms2, argv[2]);
	if(cms1.depth != cms2.depth || cms1.width != cms2.width) {
		fprintf(stderr, "The two sketches are not compatible!");
		return 1;
	}
	num = den = 0;
	for(long i = 0; i < cms1.depth * cms1.width; ++i)
	{
		if(cms1.bins[i] < cms2.bins[i]){
			num += cms1.bins[i];
			den += cms2.bins[i];
		} else if(cms1.bins[i] > cms2.bins[i]) {
			num += cms2.bins[i];
			den += cms1.bins[i];
		}
	}
	if(num == den) printf("1.0");
	else printf("%f", num/den);
	return 0;	
}
