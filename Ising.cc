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
#include <stdio.h>
#include <math.h>
#include <complex.h>


using namespace std;

double average(bool* sigma,int N){

	double aux=0;
	
	for (int i=0; i<N; i++) {
		aux+=sigma[i];
	}
	
	aux /= N;
	
	return (aux*2-1);
}


int sum(bool* sigma,int N){
	
	int aux=0;
	
	for (int i=0; i<N; i++) {
		aux+=sigma[i];
	}
	
	return (aux);
}


int main(){

	int n = 48;  //lattice side size
	int n2=n*n;
	int N=n*n*n;  //number of spins
//	int N2= N*N;
//	int N4=N*N*N*N;
	
	int Nd;  //Ndifferent

	//int N2=n*n*n;  //array size
	int random;
	int randomint;
 	bool* sigma=new bool[N]; 
	
	int i,j;
	double time=N*250000.0/2.0;
	int timestep = 1000;
	int sample = time/timestep;
	int delta[7];

	long double m,m2,m4;
	
	int s,s2,s4;
	
	double k;  //Tempreture
	
	int east[N],west[N],south[N],north[N],top[N],down[N];
	
	
	bool timec[timestep];
	for (j=0; j<timestep-1; j++) {
		timec[j] = false;
	}
	timec[timestep] = true;
	
	
	for (i=0; i<N; i++) {
		east[i]=(i+1) % N;
		west[i]=(i-1+N) % N;
		north[i]=(i+n) % N;
		south[i]=(i-n+N) % N;
		top[i]=(i+n*n) % N;
		down[i]=(i-n*n+N) % N;

	}
	
	string filename1;
	
	stringstream ss;
	
	ss << n;
	
	filename1="ising-l"+ss.str()+"-new.txt";
	std::ofstream result(filename1.c_str());   //initialize file

	
	for (k=0.215; k < 0.23; k+=0.0002) {
	
	
	 double auxm=0;
	 double auxm2=0;
	 double auxm4=0;
		
	//	long int auxs=0;
	//	long int auxs2=0;
	//	long int auxs4=0;

	
	for (i=0; i<7; i++) {
		
		if (exp (-2*k*(6-2*i)) > 1) {
			delta[i] = RAND_MAX;
		}else {
			delta[i]= exp (-2*k*(6-2*i))*RAND_MAX ;
		}

	}
		
		
		double interval =((double) RAND_MAX+1) / (double) N;
	
	for (i=0; i<N; i++) {
		random =rand();
		
		int halfmax = RAND_MAX/2;
		
		if (random < halfmax) {
			sigma[i]=false;
		} else {
			sigma[i]=true;
		}

	}	

		
	for (i=1; i < time; i++) {
	//	randomint = rand()% N;
		randomint = rand()/interval;

		Nd = (sigma[randomint] ^ sigma[east[randomint]]) + (sigma[randomint] ^ sigma[west[randomint]]) + (sigma[randomint] ^ sigma[north[randomint]]) + (sigma[randomint] ^ sigma[south[randomint]]) + (sigma[randomint] ^ sigma[top[randomint]]) + (sigma[randomint] ^ sigma[down[randomint]]);
	   
		random =rand();
		if (random < delta[Nd]) {
			sigma[randomint] = !sigma[randomint] ;
		}
	}

	
		j=0;
		for (i=1; i < time; i++) {
	//	randomint = rand()% N;
		randomint = rand()/interval;

	
	Nd = (sigma[randomint] ^ sigma[east[randomint]]) + (sigma[randomint] ^ sigma[west[randomint]]) + (sigma[randomint] ^ sigma[north[randomint]]) + (sigma[randomint] ^ sigma[south[randomint]]) + (sigma[randomint] ^ sigma[top[randomint]]) + (sigma[randomint] ^ sigma[down[randomint]]);

	random = rand();
				
		if (random < delta[Nd]) {
			sigma[randomint] = !sigma[randomint] ;
		} 

		j+=1;
		
		if (timec[j]) {
		
		j=0;
		m = average(sigma,N);
		m2=m*m;
		m4=m2*m2;
		
		auxm += m;
		auxm2 += m2;
		auxm4 += m4;
			
	/*		s = sum(sigma,N)*2-N;
			s2 = s*s;
			s4 = s2*s2;
			
			auxs += s;
			auxs2 += s2;
			auxs4 += s4;*/
			
		}
	
	}
	
	result << k << "   " << auxm / sample << "   " << auxm2 / sample << "   " << auxm4 / sample << endl;
//	result << k << "   " << auxs / (sample*N) << "   " << auxs2 / (sample*N*N) << "   " << auxs4 / (sample*N*N*N*N) << endl;

		
	cout << " l = "<< n << " k = " << k << "simulation done!" <<endl;
	
   // cout << sigma[-3] << "  " << sigma [ N-3] << endl;
	
	}
	
/*	for (k=0.22; k < 0.23; k+=0.0005) {
		
		
		long double auxm=0;
		long double auxm2=0;
		long double auxm4=0;
		
		
		for (i=0; i<7; i++) {
			delta[i]= exp (-2*k*(6-2*i));
		}
		
		for (i=0; i<N; i++) {
			random =(double)rand()/(double)RAND_MAX;
			
			if (random<0.5) {
				sigma[i]=false;
			} else {
				sigma[i]=true;
			}
			
		}	
		
		for (i=1; i < time; i++) {
			randomint = rand()% N;
			Nd = (sigma[randomint] ^ sigma[(randomint+1) % N]) + (sigma[randomint] ^ sigma[(randomint-1+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n) % N]) + (sigma[randomint] ^ sigma[(randomint-n+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n2) % N]) + (sigma[randomint] ^ sigma[(randomint-n2+N) % N]);
			random = (double)rand()/(double)RAND_MAX;
			if (random < delta[Nd]) {
				sigma[randomint] = !sigma[randomint] ;
			}
		}
		
		
		for (i=1; i < time; i++) {
			randomint = rand()% N;
			
			Nd = (sigma[randomint] ^ sigma[(randomint+1) % N]) + (sigma[randomint] ^ sigma[(randomint-1+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n) % N]) + (sigma[randomint] ^ sigma[(randomint-n+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n2) % N]) + (sigma[randomint] ^ sigma[(randomint-n2+N) % N]);
			
			random = (double)rand()/(double)RAND_MAX;
			
			if (random < delta[Nd]) {
				sigma[randomint] = !sigma[randomint] ;
			} 
			
			if (i % timestep ==0) {
				
				m = average(sigma,N);
				
				m2=m*m;
				m4=m2*m2;
				
				auxm += m;
				auxm2 += m2;
				auxm4 += m4;
				
			}
			
		}
		
		result << k << "   " << auxm / sample << "   " << auxm2 / sample << "   " << auxm4 / sample << endl;
		
		cout << " l = "<< n << " k = " << k << "simulation done!" <<endl;
		
		// cout << sigma[-3] << "  " << sigma [ N-3] << endl;
		
	}

	for (k=0.23; k < 0.25; k+=0.001) {
		
		
		long double auxm=0;
		long double auxm2=0;
		long double auxm4=0;
		
		
		for (i=0; i<7; i++) {
			delta[i]= exp (-2*k*(6-2*i));
		}
		
		for (i=0; i<N; i++) {
			random =(double)rand()/(double)RAND_MAX;
			
			if (random<0.5) {
				sigma[i]=false;
			} else {
				sigma[i]=true;
			}
			
		}	
		
		for (i=1; i < time; i++) {
			randomint = rand()% N;
			Nd = (sigma[randomint] ^ sigma[(randomint+1) % N]) + (sigma[randomint] ^ sigma[(randomint-1+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n) % N]) + (sigma[randomint] ^ sigma[(randomint-n+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n2) % N]) + (sigma[randomint] ^ sigma[(randomint-n2+N) % N]);
			random = (double)rand()/(double)RAND_MAX;
			if (random < delta[Nd]) {
				sigma[randomint] = !sigma[randomint] ;
			}
		}
		
		
		for (i=1; i < time; i++) {
			randomint = rand()% N;
			
			Nd = (sigma[randomint] ^ sigma[(randomint+1) % N]) + (sigma[randomint] ^ sigma[(randomint-1+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n) % N]) + (sigma[randomint] ^ sigma[(randomint-n+N) % N]) + (sigma[randomint] ^ sigma[(randomint+n2) % N]) + (sigma[randomint] ^ sigma[(randomint-n2+N) % N]);
			
			random = (double)rand()/(double)RAND_MAX;
			
			if (random < delta[Nd]) {
				sigma[randomint] = !sigma[randomint] ;
			} 
			
			if (i % timestep ==0) {
				
				m = average(sigma,N);
				
				m2=m*m;
				m4=m2*m2;
				
				auxm += m;
				auxm2 += m2;
				auxm4 += m4;
				
			}
			
		}
		
		result << k << "   " << auxm / sample << "   " << auxm2 / sample << "   " << auxm4 / sample << endl;
		
		cout << " l = "<< n << " k = " << k << "simulation done!" <<endl;
		
		// cout << sigma[-3] << "  " << sigma [ N-3] << endl;
		
	}*/
	
	

	
	
	
	result.close();
	
	return (0);
}
