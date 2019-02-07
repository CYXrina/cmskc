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
	    printf("argv[1]=%s\n", argv[1]);
	    if(!strcmp(argv[1], "-m"))
	    {
	        CountMinSketch cms;
	        cms_import(&cms, argv[2]);
	        int32_t freq_max = 0;
	        int32_t freq_min = INT_MAX;
	        for(size_t i = 0; i < cms.depth * cms.width; ++i)
	        {
		        if(cms.bins[i] > freq_max) freq_max = cms.bins[i];
		        if(cms.bins[i] < freq_min) freq_min = cms.bins[i];
	        }
	        printf("min = %i, MAX = %i\n", freq_min, freq_max);
	    }
	}
	return 0;
}
