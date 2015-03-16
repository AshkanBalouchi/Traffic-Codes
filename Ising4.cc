/*
 *  ising2.cc
 *  
 *
 *  Created by Ashkan Balouchi on 9/2/14.
 *  Copyright 2014 LSU. All rights reserved.
 *
 */


/*
 *  Ising.cc
 *  
 *
 *  Created by Ashkan Balouchi on 8/25/14.
 *  Copyright 2014 LSU. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>


#include <math.h>
#include <ctime>
#include <cstdlib>

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>		/* time */
#include "gsl_rng.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <map>

#include "Ising4.h"



using namespace std;



int num=10;
int n = 8;  //lattice side size
int sample = 1000000;
int samplefreq = 10;




int n2=n*n;
const int N=n*n*n;  //number of spins

int Nd;  //Ndifferent

//double random;
int randomint;
int* sigma=new int[N]; 

int* east=new int[N];
int* west=new int[N]; 
int* north=new int[N]; 
int* south=new int[N]; 
int* top=new int[N]; 
int* down=new int[N]; 



int i,j;
int timestep = N;
int timet = sample * timestep;
int initialtime = timet/10;
int initialsample = sample/10;
int totsim = sample / samplefreq;

double* delta=new double[13];

int x;

long double m,m2,m4;

int s,s2,s4;

double k;  //Tempreture



int main(int argc, char* argv[]){
	
		
	Address(east,west,north,south,top,down,n);
			
	string filename1;
	
	stringstream ss,runnum;
	
	ss << n;
	runnum << num;

	//---------------------------------timing initialisation
	
	time_t time0, *time_ptr = 0;       
	time0 = time(time_ptr); 
	
	
	
	
	filename1="ising-l"+ss.str()+"-run"+runnum.str()+".txt";
	std::ofstream result(filename1.c_str());   //initialize file
	
	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);	
	
	
	gsl_rng_set ( r, num*1000000);

	
	
	
	for (k=0.21; k < 0.24; k+=0.01) {
		
		
		double auxm=0;
		double auxm2=0;
		double auxm4=0;
		
		/*	double auxs=0;
			double auxs2=0;
			double auxs4=0;*/
		
		energygap ( delta , k );
		
		initializerandom( sigma , r , N);
		
		cout << "ave" << (double)sum(sigma,N)/N <<endl;
		
			
			evolve ( sigma, delta, N, initialsample, r, east, west, north, south, top, down);
	
	
		for (i=0; i < totsim; i++) {
			
			evolve ( sigma, delta, N, samplefreq, r, east, west, north, south, top, down);
			
			m = Ave(sigma,N);
			m2=m*m;
			m4=m2*m2;
			
			auxm += m;
			auxm2 += m2;
			auxm4 += m4;
			
			
		
		/*	 s = sum(sigma,N);
			 s2 = s*s;
			 s4 = s2*s2;
			 
			 auxs += s;
			 auxs2 += s2;
			 auxs4 += s4;*/
			
		}
		
		result << k << "   " << auxm / sample << "   " << auxm2 / sample << "   " << auxm4 / sample << endl;
	// 	result << k << "   " << auxs / (sample*N) << "   " << auxs2 / (sample*N*N) << "   " << auxs4 / (sample*N*N*N*N) << endl;
		
		
		cout << " l = "<< n << " k = " << k << "simulation done!" <<endl;
	
	}

	
	
	//----------------------- output simulation time
	
	double took = difftime(time(time_ptr),time0);
	std::cout << std::endl << "Simulation took " << took << " seconds" << std::endl << std::endl;
	
	
	
	result.close();
	
	return (0);
}

