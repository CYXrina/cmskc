all:main compare print

main:cmskc count_min_sketch
	g++ -o cmskc cmskc.o count_min_sketch.o -lz -lm
cmskc:
	g++ -c cmskc.cpp -o cmskc.o

compare:count_min_sketch
	gcc compare_cms.c count_min_sketch.o -lm -o compare_cms

print:count_min_sketch
	gcc print_cms.c count_min_sketch.o -lm -o print_cms

count_min_sketch:
	gcc -c count_min_sketch.c -o count_min_sketch.o
