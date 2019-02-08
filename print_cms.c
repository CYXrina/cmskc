#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "count_min_sketch.h"

int main(int argc, char** argv)
{
    if(argc == 2) 
    {
	    CountMinSketch cms;
	    cms_import(&cms, argv[1]);
	    cms_print(&cms);
	} 
	else 
	{
	    if(!strcmp(argv[1], "-m"))
	    {
	        CountMinSketch cms;
	        cms_import(&cms, argv[2]);
	        int32_t freq_max = 0;
	        int32_t freq_min = INT_MAX;
	        size_t idx_max, idx_min;
	        for(size_t i = 0; i < cms.depth * cms.width; ++i)
	        {
		        if(cms.bins[i] > freq_max) 
		        {
		            freq_max = cms.bins[i];
		            idx_max = i;
		        }
		        if(cms.bins[i] < freq_min) 
		        {
		            freq_min = cms.bins[i];
		            idx_min = i;
		        }
	        }
	        printf("min = %i at [%lu, %lu]\nMAX = %i at [%lu, %lu]\n", 
	            freq_min, idx_min / cms.width, idx_min % cms.width, freq_max, idx_max / cms.width, idx_max % cms.width);
	    }
	}
	return 0;
}
