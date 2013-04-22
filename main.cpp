#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <list>

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

struct particle
{
	double radius;
	double pos_x;
	double pos_y;
	double charge;
};

struct particle *parray;
list <double> **particle_position;
double iteration_time, charge_const,epsilon,size_x,size_y,polynomial[2][6];
int partition_x,partition_y,particles, iterations;


void calculate_partitions()
{
	partition_x = size_x/epsilon;
	partition_y = size_y/epsilon;
	printf("%d\n",partition_x);
}



void init()
{
	calculate_partitions();
	parray = new particle[particles];
	particle_position = new list<double>*[partition_x];
	for(volatile int i=0;i<partition_x;i++)
	{
		particle_position[i] = new list<double>[partition_y];
	}
}

void read_particles(char *fname)
{
	volatile int id = 0;
	string line;
	ifstream file (fname);
	getline(file,line);
	sscanf_s(line.c_str(),"%d %lf %lf",&particles,&size_x,&size_y);
	init();
	while(file.good())
	{
		getline(file,line);
		sscanf(line.c_str(),"%lf %lf %lf %lf",&parray[id].pos_x,&parray[id].pos_y,&parray[id].charge,&parray[id].radius);
		id++;
	}
	return;
}

void read_polynomial(char *fname)
{
	volatile int id = 0;
	string line;
	ifstream file (fname);
	getline(file,line);
	sscanf_s(line.c_str(),"%lf %lf %lf %lf %lf %lf",&polynomial[0][0],&polynomial[0][1],&polynomial[0][2],&polynomial[0][3],&polynomial[0][4],&polynomial[0][5]);
	getline(file,line);
	sscanf_s(line.c_str(),"%lf %lf %lf %lf %lf %lf",&polynomial[1][0],&polynomial[1][1],&polynomial[1][2],&polynomial[1][3],&polynomial[1][4],&polynomial[1][5]);
	return;
}



double fx(double px, double py)
{
	return 0;
}


int main(int argc, char *argv[])
{
	epsilon = atof(argv[4]);
	charge_const = atof(argv[3]);
	iterations = atoi(argv[6]);
	iteration_time = atof(argv[5]);
	read_particles(argv[1]);
	read_polynomial(argv[2]);
	particle_position[0][0].push_front(42);
	printf("%d -- Tamanhp lista[0][0]\n", particle_position[0][0].size());
	for(int i=0;i<particles;i++)
	{
		parray[i].charge = 10+i;
	}
	printf("%f %d\n", parray[2].charge,iterations);
	return 0;
}