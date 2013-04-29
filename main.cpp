#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <list>
#include <string>
#include <fstream>
#include <array>

#define DEBUG 1
#ifdef DEBUG
# define print(x) printf x
#else
# define print(x) do {} while (0)
#endif

using namespace std;

struct particle
{
	double radius;
	double pos_x;
	double pos_y;
	double vx;
	double vy;
	double charge;
};

struct point
{
	double px;
	double py;
};


struct particle *parray[2];
list <int> **particle_position;
list <array<int,2>> used;
double iteration_time, charge_const,epsilon,size_x,size_y,polynomial[2][6], ***force_vector,half_square_iterations_time;
int partition_x,partition_y,particles, iterations, thread_num;
bool turn;


void calculate_partitions()
{
	partition_x = (int)ceil(size_x/epsilon)+2;
	partition_y = (int)ceil(size_y/epsilon)+2;
	print(("%d %d\n",partition_x,partition_y));

}



void init()
{
	half_square_iterations_time = (iteration_time*iteration_time)/2.0;
	turn = 0;
	calculate_partitions();
	parray[0] = new particle[particles];
	parray[1] = new particle[particles];
	particle_position = new list<int>*[partition_x];
	for(volatile int i=0;i<partition_x;i++)
	{
		particle_position[i] = new list<int>[partition_y];
	}
	force_vector = new double**[thread_num];
	for(int i=0;i<thread_num;i++)
	{
		force_vector[i] = new double*[particles+1];
		for(int j=0;j<particles+1;j++)
		{
			force_vector[i][j] = new double[2];
		}
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
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",&parray[0][id].pos_x,&parray[0][id].pos_y,&parray[0][id].vx, &parray[0][id].vy, &parray[0][id].charge,&parray[0][id].radius);
		id++;
	}
	for(int i=0;i<particles;i++)
	{
		parray[1][i].pos_x = parray[0][i].pos_x;
		parray[1][i].pos_y = parray[0][i].pos_y;
		parray[1][i].vx = parray[0][i].vx;
		parray[1][i].vy = parray[0][i].vy;
		parray[1][i].charge = parray[0][i].charge;
		parray[1][i].radius = parray[0][i].radius;
	}
	return;
}


array<int,2> find_particle_position_in_matrix(int id)
{
	array<int,2> var;
	var[0] = (int)floor(parray[turn][id].pos_x/epsilon)+1;
	var[1] = (int)floor(parray[turn][id].pos_y/epsilon)+1;
	return var;
}


void calculate_particle_position(int id)
{
	array<int,2> var = find_particle_position_in_matrix(id);
#pragma omp critical
	{
	used.push_front(var);
	particle_position[var[0]][var[1]].push_front(id);
	}
}


/* Usando a matrix extendida */
void print_particle_position()
{
	for(int i = 0;i<partition_x;i++)
	{
		for(int j=0;j<partition_y;j++)
		{
			if(!particle_position[i][j].empty())
			{
				for(list<int>::iterator it=particle_position[i][j].begin(); it!=particle_position[i][j].end(); it++)
				{
					printf("M[%d][%d] --> %d\n", i,j,*it);
				}
			}
		}
	}
}

void print_parray()
{
	for(int i=0;i<particles;i++)
	{
		printf("Particula[%d] --> X = %lf      Y = %lf\n",i,parray[turn][i].pos_x,parray[turn][i].pos_y);
	}
}

void print_used()
{
	int i = 0;
	list<array<int,2>> uni;
	used.unique();
	used.sort();
	used.unique();
	print(("Size: %d\n",used.size()));
	for(list<array<int,2>>::iterator it=used.begin(); it!=used.end(); it++)
	{
		printf("Used:  %d  %d\n", ((array<int,2>)*it)[0],((array<int,2>)*it)[1]);
		i++;
	}
	print(("i: %d\n",i));
}




double distance(double p1x,double p1y, double p2x, double p2y)
{
	double px = p1x-p2x;
	double py = p1y-p2y;
	return sqrt(px*px+py*py);
}

/* Note que se as particulas tiverem a mesma carga, entao a força resultante será positiva */
double force(double dist, double q1, double q2)
{
	return (charge_const*q1*q2)/(dist*dist);
}

int collision_detection(int id,int thread_id, int position_x, int position_y, int force_size)
{
	particle p,paux;
	list<int>::iterator it;
	list<int> aux;
	double dist,ang,f;
	aux = particle_position[position_x][position_y];
	p = parray[turn][id];
	for(it = aux.begin(); it != aux.end(); it++)
	{
		paux = parray[turn][*it];
		dist = distance(p.pos_x,p.pos_y,paux.pos_x,paux.pos_y);
		if(dist < epsilon && dist > (p.radius+paux.radius))
		{
			f = force(dist,p.charge,paux.charge);
			ang = atan2(p.pos_y-paux.pos_y,p.pos_x-paux.pos_x);
			force_vector[thread_id][force_size][0] = f*cos(ang);
			force_vector[thread_id][force_size][1] = f*sin(ang);
			force_size++;
			print(("Forca resultante em %d com dist = %lf\nX: %lf Y: %lf\n",id,dist,force_vector[thread_id][force_size-1][0],force_vector[thread_id][force_size-1][1]));
		}
	}
	return force_size;
}

void final_position(int id, int thread_id)
{
	int divx,divy;
	double sx,sy,colx,coly;
	sx = parray[turn][id].pos_x+parray[turn][id].vx*iteration_time+force_vector[thread_id][0][0]*half_square_iterations_time;
	sy = parray[turn][id].pos_y+parray[turn][id].vy*iteration_time+force_vector[thread_id][0][1]*half_square_iterations_time;
	/* Atualiza a força no final */
	divx = 0;
	divy = 1;
	while((sx > size_x || sx < 0) || (sy > size_y || sy < 0))
	{

	}
}

void calculate_forces(int id)
{
	int thread_id,force_size = 0;
	array<int,2> position = find_particle_position_in_matrix(id);
	thread_id = omp_get_thread_num();
	force_vector[thread_id][0][0] = 0;
	force_vector[thread_id][0][1] = 0;
	print(("Thread ID %d\n",thread_id));
	force_size = collision_detection(id,thread_id, position[0],position[1], force_size);
	force_size = collision_detection(id,thread_id, position[0]-1,position[1], force_size);
	force_size = collision_detection(id,thread_id, position[0]-1,position[1]+1, force_size);
	force_size = collision_detection(id,thread_id, position[0],position[1]+1, force_size);
	force_size = collision_detection(id,thread_id, position[0]+1,position[1]+1, force_size);
	force_size = collision_detection(id,thread_id, position[0]+1,position[1], force_size);
	force_size = collision_detection(id,thread_id, position[0]+1,position[1]-1, force_size);
	force_size = collision_detection(id,thread_id, position[0],position[1]-1, force_size);
	force_size = collision_detection(id,thread_id, position[0]-1,position[1]-1, force_size);
	for(int i=1;i<force_size;i++)
	{
		force_vector[thread_id][0][0] += force_vector[thread_id][id][0];
		force_vector[thread_id][0][1] += force_vector[thread_id][id][1];
	}
	final_position(id,thread_id);
}

int main(int argc, char *argv[])
{
	//int teste;
	thread_num = omp_get_num_procs();
	if(argc == 7) thread_num = atoi(argv[6]);
	omp_set_num_threads(thread_num);
	epsilon = atof(argv[3]);
	charge_const = atof(argv[2]);
	iterations = atoi(argv[5]);
	iteration_time = atof(argv[4]);
	read_particles(argv[1]);

	for(int i=0;i<particles;i++) calculate_particle_position(i);
	print_particle_position();
	print_parray();
	print_used();
	for(int i=0;i<iterations;i++)
	{
		#pragma omp parallel for schedule(static) 
		for(int id=0;id<particles;id++)
		{
			calculate_forces(id);
		}
		/* Talvez nao seja necessario coloca single e barrier pois ja esta fora da pare paralela? */
		//#pragma omp single
		//{
			turn = !turn;
		//}
		//#pragma omp barrier
	}
	return 0;
}