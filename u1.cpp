#include <iostream>
#include <math.h> 
#include <time.h> 
#include <random>
#include <vector> 
#include <fstream>
#include <omp.h>

#include "timer.h"

using namespace std;


int DIRS;
int TDir;
int volume;
int half_volume;
int spatial_volume;
int Grid[4];
int iter = 0;
double accept_ratio = 0.;

 void indexEO(int id, int parity, int x[4]){
	int za = (id / (Grid[0]/2));
	int zb =  (za / Grid[1]);
	x[1] = za - zb * Grid[1];
	x[3] = (zb / Grid[2]);
	x[2] = zb - x[3] * Grid[2];
	int xodd = (x[1] + x[2] + x[3] + parity) & 1;
	x[0] = (2 * id + xodd )  - za * Grid[0];
 }

int indexId(int x[4]){
 return  ((((x[3] * Grid[2] + x[2]) * Grid[1]) + x[1] ) * Grid[0] + x[0]);
}
int indexId(int x[4], int parity, int dir){
	int id = ((((x[3] * Grid[2] + x[2]) * Grid[1]) + x[1] ) * Grid[0] + x[0]) >> 1;
 return id + parity * half_volume + volume * dir;
}
int indexId(int x[4], int dir){
	int id = ((((x[3] * Grid[2] + x[2]) * Grid[1]) + x[1] ) * Grid[0] + x[0]) >> 1;
	int parity = (x[0] + x[1] + x[2] +x[3]) & 1;
 return id + parity * half_volume + volume * dir;
}

int indexEO_neg(const int id, int parity, int mu, int lmu){
	int x[4];
	indexEO(id, parity, x);
	x[mu] = (x[mu]+lmu+Grid[mu]) % Grid[mu];
	
	int pos = ((((x[3] * Grid[2] + x[2]) * Grid[1]) + x[1] ) * Grid[0] + x[0]) >> 1;
	int oddbit = (x[0] + x[1] + x[2] +x[3]) & 1;
	pos += oddbit  * half_volume;
	return pos;
}
int indexEO_neg(const int id, int parity, int mu, int lmu, int nu, int lnu){
	int x[4];
	indexEO(id, parity, x);
	x[mu] = (x[mu]+lmu+Grid[mu]) % Grid[mu];
	x[nu] = (x[nu]+lnu+Grid[nu]) % Grid[nu];

	int pos = ((((x[3] * Grid[2] + x[2]) * Grid[1]) + x[1] ) * Grid[0] + x[0]) >> 1;
	int oddbit = (x[0] + x[1] + x[2] +x[3]) & 1;
	pos += oddbit  * half_volume;
	return pos;
}


int Index_4D_Neig_EO(const int x[], const int dx[], const int X[4]) {
	int y[4];
	for (int i=0; i<4; i++) y[i] = (x[i] + dx[i] + X[i]) % X[i];
	int idx = (((y[3]*X[2] + y[2])*X[1] + y[1])*X[0] + y[0]) >> 1;
	return idx;
}



/*
default_random_engine seed;
uniform_real_distribution<double> rand02(0., 2.);
uniform_real_distribution<double> rand01(0,1);*/

std::random_device device;
//std::mt19937 generator(device());
std::mt19937 *generator;//(device());
uniform_real_distribution<double> rand02(0., 2.);
uniform_real_distribution<double> rand01(0,1);




void staple(const double *lat, const int id, const int parity, const int mu, double &stapleRe, double &stapleIm){
	stapleRe = 0., stapleIm = 0.;	
	int idmu1 = indexEO_neg(id, parity, mu, 1);			
	for(int nu = 0; nu < DIRS; nu++)  if(mu != nu) {
		double upperStaple = lat[idmu1 + volume * nu];
		upperStaple -= lat[indexEO_neg(id, parity, nu, 1) + volume * mu];
		upperStaple -= lat[id + parity * half_volume + nu * volume];
		
		double lowerStaple = -lat[indexEO_neg(id, parity, mu, 1, nu, -1) + volume * nu];	
		lowerStaple -= lat[indexEO_neg(id, parity, nu, -1) + volume * mu];	
		lowerStaple += lat[indexEO_neg(id, parity, nu, -1) + volume * nu];	
		
		stapleRe += cos(upperStaple) + cos(lowerStaple);
		stapleIm += sin(upperStaple) + sin(lowerStaple);
	}
}



void metropolis(double *lat, double beta){
	for(int parity = 0; parity < 2; ++parity)
	for(int mu = 0; mu < DIRS; mu++){
		#pragma omp parallel for
		for(int id = 0; id < volume/2; ++id){
			double phase_old = lat[id + parity * half_volume + mu * volume];
			int idmu1 = indexEO_neg(id, parity, mu, 1);
			double stapleRe = 0., stapleIm = 0.;
			staple(lat, id, parity, mu, stapleRe, stapleIm);			
			double r = std::sqrt( stapleRe*stapleRe + stapleIm*stapleIm );
			double t2 = atan2(stapleIm, stapleRe);

			double new_phase = M_PI * rand02(generator[omp_get_thread_num()]);
			double b = rand01(generator[omp_get_thread_num()]);

			double S1 = cos(phase_old + t2);
			double S2 = cos(new_phase + t2);
			double dS = exp(beta*r*(S2-S1));
			if(dS > b){
				lat[id + parity * half_volume + mu * volume] = new_phase;
				accept_ratio += 1.;
			}
		}
	}
}




void overrelaxation(double *lat, double beta){
	for(int parity = 0; parity < 2; ++parity)
	for(int mu = 0; mu < DIRS; mu++){
		#pragma omp parallel for
		for(int id = 0; id < volume/2; ++id){
			double stapleRe = 0., stapleIm = 0.;
			staple(lat, id, parity, mu, stapleRe, stapleIm);
			int pos = id + parity * half_volume + mu * volume;
			double phase_old = lat[pos];
			double t2 = atan2(stapleIm, stapleRe);
			double new_phase = fmod(6.* M_PI - phase_old - 2. * t2, 2.* M_PI);
			lat[pos] = new_phase;
		}
	}
}








void plaquette(double *lat, double *plaq){
	for(int i = 0; i < 2; ++i) plaq[i] = 0.;
	for(int parity = 0; parity < 2; ++parity){
		#pragma omp parallel for reduction(+:plaq[:2])
		for(int id = 0; id < volume/2; ++id){
			for(int mu = 0; mu < DIRS - 1; mu++){	
				double tmp = lat[id + parity * half_volume + mu * volume];
				int idmu1 = indexEO_neg(id, parity, mu, 1);
				for (int nu = (mu+1); nu < DIRS; nu++){			
					double plaqi = tmp;
					plaqi += lat[idmu1 + volume * nu];
					plaqi -= lat[indexEO_neg(id, parity, nu, 1) + volume * mu];
					plaqi -= lat[id + parity * half_volume + nu * volume];
					
					plaq[0] += cos(plaqi);
					plaq[1] += sin(plaqi);	
				}
			}
		}
	}
	int numplaqs = 6; //DIRS=4 3D+1
	if(DIRS==2) numplaqs = 1.;
	else if(DIRS==3) numplaqs = 3.;
	double norm = 1. / double(volume * numplaqs);
	for(int i = 0; i < 2; ++i) plaq[i] *= norm;
}

void polyakov(double *lat, double *poly){
	for(int i = 0; i < 2; ++i) poly[i] = 0.;
	for(int parity = 0; parity < 2; ++parity){
		#pragma omp parallel for reduction(+:poly[:2])
		for(int id = 0; id < spatial_volume/2; ++id){
			int x[4];
			indexEO(id, parity, x);
			double tmp = 0.;
			for(x[TDir] = 0; x[TDir] < Grid[TDir]; ++x[TDir])
				tmp += lat[ indexId(x, TDir) ];
			poly[0] += cos(tmp);
			poly[1] += sin(tmp);
		}
	}
	double norm = 1. / double(spatial_volume);
	for(int i = 0; i < 2; ++i) poly[i] *= norm;
}





int main(){
	Timer a0;
	a0.start();
	
	//omp_set_num_threads(1);
	
	int numthreads = 0;
	#pragma omp parallel
	numthreads = omp_get_num_threads();
	cout << "Number of threads: " << numthreads << endl;
	
	//create one RNG per thread
	generator = new std::mt19937[numthreads];
	for(int i = 0; i < numthreads; ++i) generator[i].seed(time(NULL)*(i+1));
	
	
	DIRS = 4;
	TDir = DIRS - 1;
	int ls = 8; //The number of points in each direction must be an even number!!!!!!!!!
	int Nx=ls;
	int Ny=ls;
	int Nz=ls;
	int Nt=ls;
	for(int i = 0; i < 4; ++i) Grid[i] = 1;
	Grid[0] = Nx;
	if(DIRS==2) Grid[1] = Nt;
	else if(DIRS > 2) Grid[1] = Ny;
	if(DIRS==3) Grid[2] = Nt;
	else if(DIRS > 3)  Grid[2] = Nz;
	if(DIRS==4) Grid[3] = Nt;
	
	
	double beta = 2.;
	bool hotstart = false;
	int ovrn = 3;
	
	volume = 1;
	for(int i = 0; i < 4; ++i) volume *= Grid[i];	
	half_volume = volume / 2;
	spatial_volume = 1;
	for(int i = 0; i < TDir; ++i) spatial_volume *= Grid[i];
	
	// creates the lattice array and initializes it to 0, cold start
	double *lat = new double[volume*DIRS](); 
	
	
	if(hotstart) {
		//Initializes lattice array with random phase (hot start) between 0-2Pi
		#pragma omp parallel for	
		for(int id = 0; id < volume*DIRS; ++id) 
			lat[id] = M_PI * rand02(generator[omp_get_thread_num()]);
	}
	double plaq[2];
	double poly[2];
	plaquette(lat, plaq);
	for(iter = 1; iter <= 5000; ++iter){
		metropolis(lat, beta);
		for(int ovr = 0; ovr < ovrn; ++ovr)
			overrelaxation(lat, beta);
			
		plaquette(lat, plaq);
		polyakov(lat, poly);
		if( (iter%10)==0){
			cout << "iter: " << iter << " \tplaq: " << 1.-plaq[0] << "\t" << plaq[1] << endl;
			cout << "           " << " \tL: " << poly[0] << "\t" << poly[1] << "\t|L|: " << sqrt(poly[0]*poly[0]+poly[1]*poly[1]) << endl;
		}
	}
	cout << "Acceptation ratio: " << accept_ratio/double(volume*DIRS*iter) << endl;
	
	delete[] lat;
	delete[] generator;	
	a0.stop();
	std::cout << "Time: " << a0.getElapsedTime() << endl;
	return 0;
	
}
