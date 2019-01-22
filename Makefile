all:sketcher

sketcher:main count_min_sketch
	g++ -o cmskc main.o count_min_sketch.o -lz

main:
	g++ -c main.cpp -o main.o

count_min_sketch:
	gcc -c count_min_sketch.c -o count_min_sketch.o
