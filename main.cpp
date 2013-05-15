// Caio de Freitas Valente 6552442
// Geraldo Castro Zampoli  6552380



#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <list>
#include <string>
#include <fstream>
#include <array>
#include <math.h>


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
list <int> **particle_position[2];
double iteration_time, charge_const,epsilon,size_x,size_y,polynomial[2][6], ***force_vector,half_square_iterations_time;
int partition_x,partition_y,particles, iterations, thread_num;
bool turn;


/* Calcula o tamanho da matriz de colisoes que devera ser alocada */
void calculate_partitions()
{
	partition_x = (int)ceil(size_x/epsilon)+2;
	partition_y = (int)ceil(size_y/epsilon)+2;
	print(("%d %d\n",partition_x,partition_y));

}


/* Inicializacao das estruturas de dados */
void init()
{
	half_square_iterations_time = (iteration_time*iteration_time)/2.0;
	turn = 0;
	calculate_partitions();
	parray[0] = new particle[particles];
	parray[1] = new particle[particles];
	particle_position[0] = new list<int>*[partition_x];
	particle_position[1] = new list<int>*[partition_x];
	for(volatile int i=0;i<partition_x;i++)
	{
		particle_position[0][i] = new list<int>[partition_y];
		particle_position[1][i] = new list<int>[partition_y];
	}
	force_vector = new double**[thread_num];
	for(int i=0;i<thread_num;i++)
	{
		force_vector[i] = new double*[particles];
		for(int j=0;j<particles;j++)
		{
			force_vector[i][j] = new double[2];
		}
	}
}

/* Funcao de entrada */
void read_particles(char *fname)
{
	volatile int id = 0;
	string line;
	ifstream file (fname);
	getline(file,line);
	sscanf(line.c_str(),"%d %lf %lf",&particles,&size_x,&size_y);
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

/* Determina a posicao de uma particula na matriz de colisoes */
array<int,2> find_particle_position_in_matrix(int id,bool turn)
{
	array<int,2> var;
	if(parray[turn][id].pos_x == size_x) var[0] = partition_x-2;
	else var[0] = (int)floor(parray[turn][id].pos_x/epsilon)+1;
	if(parray[turn][id].pos_y == size_y) var[1] = partition_y-2;
	else var[1] = (int)floor(parray[turn][id].pos_y/epsilon)+1;
	return var;
}


/* Preenche a matriz de lista, que é usada para facilitar a deteccao de colisoes */
void calculate_particle_position(int id,bool turn)
{
	array<int,2> var = find_particle_position_in_matrix(id,turn);
	print(("Calculate particle position -- "));
	#pragma omp critical
	{
		particle_position[turn][var[0]][var[1]].push_front(id);
	}
	print(("calculate particle OK!\n"));
}


/* Usando a matriz extendida */
void print_particle_position()
{
	for(int i = 0;i<partition_x;i++)
	{
		for(int j=0;j<partition_y;j++)
		{
			if(!particle_position[turn][i][j].empty())
			{
				for(list<int>::iterator it=particle_position[turn][i][j].begin(); it!=particle_position[turn][i][j].end(); it++)
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

void print_f()
{
	printf("%d %lf %lf\n",particles,size_x,size_y);
	for(int i=0;i<particles;i++)
	{
		printf("%lf %lf %lf %lf %lf %lf\n",parray[turn][i].pos_x,parray[turn][i].pos_y,parray[turn][i].vx,parray[turn][i].vy,parray[turn][i].charge,parray[turn][i].radius);
	}
}



/* Distancia entre duas particulas */
double distance(double p1x,double p1y, double p2x, double p2y)
{
	double px = p1x-p2x;
	double py = p1y-p2y;
	return sqrt(px*px+py*py);
}

/* Calcula uma forca sobre uma particula */
double force(double dist, double q1, double q2)
{
	return (charge_const*q1*q2)/(dist*dist);
}


/* Verifica se ha colisao entre a particula id, e as particulas que estao em position_x, position_y preenche o vetor de forcas 
e retorna o numero de forcas atualizado agindo sobre essa particula */
int collision_detection(int id,int thread_id, int position_x, int position_y, int force_size)
{
	particle p,paux;
	list<int>::iterator it;
	list<int> aux;
	double dist,ang,f;
	print(("Collision detection -- "));
	print(("id %d - tid %d - posx %d - posy %d - fsize %d  ",id,thread_id,position_x,position_y,force_size));
	aux = particle_position[turn][position_x][position_y];
	p = parray[turn][id];
	int b = -5;
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
		b = *it;
	}
	print(("collision OK!\n"));
	return force_size;
}

/* Calcula a posicao e velocidade final de uma particula */
void final_position(int id, int thread_id)
{
	bool divx,divy;
	double sx,sy,vx,vy;
	print(("Final position -- "));
	sx = parray[turn][id].pos_x+parray[turn][id].vx*iteration_time+force_vector[thread_id][0][0]*half_square_iterations_time;
	sy = parray[turn][id].pos_y+parray[turn][id].vy*iteration_time+force_vector[thread_id][0][1]*half_square_iterations_time;
	divx = divy = false;
	while((sx > size_x || sx < 0) || (sy > size_y || sy < 0))
	{
		if(sx>=size_x)
		{
			sx = (2*size_x - sx)/2.0;
			divx=true;
		}
		else if(sx<=0)
		{
			sx = (-sx)/2.0;
			divx=true;
		}
		if(sy>=size_y)
		{
			sy = (2*size_y -sy)/2.0;
			divy=true;
		}
		else if(sy<=0)
		{
			sy = (-sy)/2.0;
			divy=true;
		}
	}
	vx =  (parray[turn][id].vx+force_vector[thread_id][0][0]*iteration_time);
	vy =  (parray[turn][id].vy+force_vector[thread_id][0][1]*iteration_time);
	if(divx)
	{
		vx = -vx/2.0;
	}
	if(divy) 
	{
		vy = -vy/2.0; 
	}
	parray[!turn][id].pos_x = sx;
	parray[!turn][id].pos_y = sy;
	parray[!turn][id].vx = vx;
	parray[!turn][id].vy = vy;
	print(("final Ok!\n"));
	calculate_particle_position(id,!turn);
}


/* Chama funcoes para verificar forcas que agem numa particula e sua posicao final */
void calculate_forces(int id)
{
	int thread_id,force_size = 0;
	array<int,2> position = find_particle_position_in_matrix(id,turn);
	print(("Calculate forces"));
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
	print(("calculate Ok!\n"));
}

int main(int argc, char *argv[])
{
	thread_num = omp_get_num_procs();
	if(argc == 7) thread_num = atoi(argv[6]);
	omp_set_num_threads(thread_num);
	epsilon = atof(argv[3]);
	charge_const = atof(argv[2]);
	iterations = atoi(argv[5]);
	iteration_time = atof(argv[4]);
	read_particles(argv[1]);

	#pragma omp parallel for schedule(static)
	for(int i=0;i<particles;i++) calculate_particle_position(i,turn);
	for(int i=0;i<iterations;i++)
	{
		#pragma omp parallel for schedule(static)
		for(int id=0;id<particles;id++)
		{
			calculate_forces(id);
		}
		#pragma omp parallel for schedule(static)
		for(int m = 0;m<partition_x;m++)
		{
			for(int n=0;n<partition_y;n++)
			{
				while(!particle_position[turn][m][n].empty())
				{
					particle_position[turn][m][n].pop_front();
				}
			}
		}
		turn = !turn;
	}
	print_f();
	return 0;
}