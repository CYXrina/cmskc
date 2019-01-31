#include <stdio.h>
#include "count_min_sketch.h"

int main(int argc, char** argv)
{
	CountMinSketch cms;
	cms_import(&cms, argv[1]);
	for(int i = 0; i < cms.width; ++i)
	{
	    for(int j = 0; j < cms.depth; ++j) printf("%i ", cms.bins[i * cms.depth + j]);
	    printf("\n");
	}
	printf("width = %i | depth = %i\n", cms.width, cms.depth);
	return 0;	
}
