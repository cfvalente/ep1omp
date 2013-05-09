prog: main.cpp
	g++ -Wall -ansi -O3 main.cpp -march=native -lm -fopenmp -std=c++0x -o prog