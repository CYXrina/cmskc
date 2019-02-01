#include <stdio.h>
#include "count_min_sketch.h"

int main(int argc, char** argv)
{
	CountMinSketch cms;
	cms_import(&cms, argv[1]);
	cms_print(&cms);
	return 0;	
}
