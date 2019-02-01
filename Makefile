all:builder compare print

builder:build_cms count_min_sketch
	g++ -o build_cms build_cms.o count_min_sketch.o -lz -lm

build_cms:
	g++ -c build_cms.cpp -o build_cms.o

compare:count_min_sketch
	gcc compare_cms.c count_min_sketch.o -lm -o compare_cms

print:count_min_sketch
	gcc print_cms.c count_min_sketch.o -lm -o print_cms

count_min_sketch:
	gcc -c count_min_sketch.c -o count_min_sketch.o
