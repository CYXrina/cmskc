all:builder sketcher compare print

builder:build_cms count_min_sketch
	g++ -o build_cms build_cms.o count_min_sketch.o -lz -lm

sketcher:yhgkdb fmm_sketch hllmh count_min_sketch
	g++ -o sketcher yhgkdb.o count_min_sketch.o fmm_sketch.o hllmh.o -lz -lm

build_cms:
	g++ -c build_cms.cpp -o build_cms.o

compare:count_min_sketch
	gcc compare_cms.c count_min_sketch.o -lm -o compare_cms

print:count_min_sketch
	gcc print_cms.c count_min_sketch.o -lm -o print_cms

count_min_sketch:
	gcc -c count_min_sketch.c -o count_min_sketch.o

yhgkdb:
	g++ -c yhgkdb.cpp -o yhgkdb.o

fmm_sketch:
	gcc -c fmm_sketch.c -o fmm_sketch.o

hllmh:
	gcc -c hllmh.c -o hllmh.o
