/*
 *  IsingEnergyMagnetization.cc
 *  
 *
 *  Created by Ashkan Balouchi on 10/5/14.
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



int num=20;
int n = 48;  //lattice side size
int sample = 500000;
int samplefreq = 10;

int maxe=1500000;

double* m=new double[maxe];
double* absm=new double[maxe];
double* m2=new double[maxe];
double* m4=new double[maxe];
double* c=new double[maxe];

double mag,mag2;



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
int initialtime = timet/2;
int initialsample = sample/2;
int totsim = sample / samplefreq;

double* delta=new double[13];

int x;

// long double m,m2,m4;

int s,s2,s4;

double k;  //Tempreture

int E;  //Energy
double M;  //Magentization



int main(int argc, char* argv[]){
	
	
	for (i=0; i<maxe; i++) {
		m[i]=0;
		absm[i]=0;
		m2[i]=0;
		m4[i]=0;
		c[i]=0;
	}
	
	
	Address(east,west,north,south,top,down,n);
	
	string filename1,filename2;
	
	stringstream ss,runnum;
	
	ss << n;
	runnum << num; 
	
	//---------------------------------timing initialisation
	
	time_t time0, *time_ptr = 0;       
	time0 = time(time_ptr); 
	
	
	
	
	/*	filename1="ising-l"+ss.str()+"-run"+runnum.str()+".txt";
	 std::ofstream result(filename1.c_str());   //initialize file*/
	
	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);	
	
	
	gsl_rng_set ( r, num*1000000);
	
	
	for (k=0.221654; k < 0.223; k+=0.004) {
		
		stringstream ssk;
		ssk << k;
		
		
	//	filename1="ising-EHistogram-l"+ss.str()+"-k"+ssk.str()+"-run"+runnum.str()+".txt";
	//	std::ofstream result(filename1.c_str());   //initialize file
		
		
		filename2="ising-EMHL-l"+ss.str()+"-k"+ssk.str()+"-run"+runnum.str()+".txt";
		std::ofstream result2(filename2.c_str());   //initialize file

		
		double auxm=0;
		double auxm2=0;
		double auxm4=0;
		
		/*	double auxs=0;
		 double auxs2=0;
		 double auxs4=0;*/
		
		energygap ( delta , k );
		
		initializerandom( sigma , r , N);
		
		cout << "ave" << (double)sum(sigma,N)/N <<endl;
		
		E = Energy (sigma, k, N, east, west, north, south, top, down);
		M = sum(sigma,N);
		
		cout << "Energy0: " << E << "  M0: " << M << endl;		
		
		evolve (sigma, delta, N, initialsample, r, east, west, north, south, top, down);
		
		E = Energy (sigma, k, N, east, west, north, south, top, down);
		M = sum(sigma,N);
		
		cout << "Energy: " << E << "  M: " << M << endl;	
		
		for (i=0; i < totsim; i++) {
			
			evolveEnergySpin ( sigma, delta, N, E, M, samplefreq, r, east, west, north, south, top, down);
			
			//evolveEnergy (sigma, delta, N, E, samplefreq, r, east, west, north, south, top, down);
			
			mag = (double) M/ (double) N;
			mag2 = mag*mag;
			m[-E] += mag;
			m2[-E] += mag2;
			m4[-E] += mag2*mag2;
			absm[-E] += abs(mag);
			c[-E] += 1;
			
			//result << E << endl;
		}
		
		//result << k << "   " << auxm / sample << "   " << auxm2 / sample << "   " << auxm4 / sample << endl;
		// 	result << k << "   " << auxs / (sample*N) << "   " << auxs2 / (sample*N*N) << "   " << auxs4 / (sample*N*N*N*N) << endl;
		
		
		cout << " l = "<< n << " k = " << k << "simulation done!" <<endl;
	//	result.close();
		
		for (i=0; i<maxe; i++) {
			
			if (c[i]>0) {
				
				result2 << -i << "		" << c[i] << "	 " << m[i]/c[i] << "		"<< absm[i]/c[i] << "	" << m2[i]/c[i] << "		" << m4[i]/c[i] << endl;
			}
		}
	
		
		result2.close();
	}
	
	
	
	
	//----------------------- output simulation time
	
	double took = difftime(time(time_ptr),time0);
	std::cout << std::endl << "Simulation took " << took << " seconds" << std::endl << std::endl;
	
	
	
	//	result.close();
	
	return (0);
}

