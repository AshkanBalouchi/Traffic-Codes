/*
 *  JamHistFreeCarsDensity.cc
 *  
 *
 *  Created by Ashkan Balouchi on 1/10/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 *  
 */

// This Code if for time evolution, jam hisogram, free cars!



//#include "bit.h"
#include "traffic_obj.h"
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
#include <fftw3.h>

using namespace std;


int main(){
	
	
	
	int num=2; 
	
	int n=30;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di=0.136;    //initial density
	double dm=0.1405;    //max density 
	double dd=0.0005;    //density step
	double p;     //probability of break
	int vmax = 5;    //max velocity
	int m;     //number of cars
	int totaltime=200000;        //number of totaltime steps
	int initime=400000;     //number of initialization time steps
	int bettime=1000000;     //number of initialization time steps
	int* x=new int[N];       //positions
	int* v=new int[N];       //velocities
	double* vd=new double[N];    //gaps aux
	int* g=new int[N];       //gaps
	double* gapd=new double[N];    //gaps aux
	double auxgap[N];
	double auxv[N];
	int jam,jamaux;
	int number,numaux;
	//	int jamnum[N];
	
	double* dden=new double[N]; 
	double Gaptotal;
	int* bit=new int[N]; 
	double* fw=new double[N];
	
	
	double gap0,gap02;
	
	
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=5;
	int count;
	int k;
	double rho;
	
	
	int jamcount,fjamcount;
	int jamcar[1000000];
	int jamlength[1000000];
	
	
	int len;  //Total length of the jams
	int cars;  //Total number of cars in jam
	
	
	
	//Fourier needs:
	
	int time=floor(totaltime/timestep);
	
	
	//File Names:
	//-------------------------
	string filename1,filename2,filename3,filename4,filename5,filename6,filename7,filename8,filename9,filename10;
	//-------------------------
	
	
	stringstream sss,ssl,runnum;
	
	runnum << num;
	sss << vmax;
	ssl << n;
	
	p=0.1;
	
	
	rho=(1-p)/(vmax+1-2*p);
	rho=round(rho,5);
	
	cout << " Transition density = " << rho << endl;
	
	
	cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
	
	
	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);	
	
	int num2 = num*1000000;
	
	gsl_rng_set (r,num2);
	
	
	for (d=di; d<dm; d+=dd) {
		
		
		stringstream ss1;
		ss1 << d;
		
		
		filename4="JamLengthTimeEvo-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+".txt";
		std::ofstream resultjt(filename4.c_str());   //initialize file
		
		filename5="NumCarsInJamTimeEvo-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+".txt";
		std::ofstream resultnt(filename5.c_str());   //initialize file
		
		
		filename7="NumberOfJamsHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resulthist(filename7.c_str());   //initialize file
		
		filename8="JamLengthHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resultl(filename8.c_str());   //initialize file

		filename9="JamCarsHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resultn(filename9.c_str());   //initialize file
		
		
		filename6="FreeCarsDensity-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+".txt";
		std::ofstream resultd(filename6.c_str());   //initialize file

		
		
		m=floor(d*N);
		
		initializeuniform(N,d,x,v);
		
		
		for (int t=0; t<initime; t++) {  
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v,r);
			MoveCars(N,d,vmax,x,v,g);
			
		}  //initializing simulation
		
		jamaux = 0 ;
		numaux = 0;
		
		for (int t=0; t<totaltime; t++) {  //main simulation
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v,r);
			MoveCars(N,d,vmax,x,v,g);
			
			ClusterCount(N, d, vmax, jamcount, fjamcount, jamcar, jamlength, x, v, g);
			
			resulthist << jamcount << endl;
			
			len=0;
			cars=0;
			
			for (int i=1; i<jamcount+fjamcount+1; i++) {
				if (jamlength[i] > 0) {
					resultl << jamlength[i] << endl;
					len+=jamlength[i];
					resultn << jamcar[i] << endl;
					cars+=jamcar[i];
	
				}
			} 
			
			resultd << t << "	" << ((N*d)-cars)/(N-len) << endl;
			
			
			jam=0;
			number=0;
			
			for (int i=0; i<m; i++) {
				if (g[i] < vmax/2+1) {
					jam += g[i];
					number += 1;
				}
				
			}
			
			//			jamaux += jam;
			//			numaux += number;
			
			resultjt << t << "	" << jam << endl;
			resultnt << t << "	" << number << endl;
			
			
			
		}			
		   

		
		resulthist.close();
		resultn.close();
		resultl.close();
		
		resultd.close();
		
		resultnt.close();
		resultjt.close();

		
		cout << "l = " << n  << "  num = " << num << "  d = " << d << "  simulation done!"<<endl;
		
	}
	
	

	
	
	
	
	return (0);
	
	
}










