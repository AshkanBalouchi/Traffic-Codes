/*
 *  FindGap1.cc
 *  
 *
 *  Created by Ashkan Balouchi on 12/18/13.
 *  Copyright 2013 LSU. All rights reserved.
 *
 */


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
//#include <fftw3.h>

using namespace std;



int main(){
	
	
	int N=50000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	double m;     //number of cars
	int totaltime=500000;        //number of totaltime steps
	int* x=new int[N]; 
	int* xb=new int[N]; 
	int* v=new int[N]; 
	int* g=new int[N]; 
	double* gapd=new double[N]; 
	
	int* nb=new int[N]; 
	
	double* nbdis=new double[1000]; 
	//double* gapdf=new double[N]; 
	double gap0,gap02;
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=1;
	int count,countb;
	int i,j,k;
	double rho;
	
	
	
	int time=floor(totaltime/timestep);
	
	stringstream sss;
	sss << vmax;
	
	for (p=0.1; p<0.2; p+=0.1) {
		
		
		string filename1,filename2,filename3,filename4;
		vmax=9;
		rho=(1-p)/(vmax+1-2*p);
		rho=round(rho,3);
		cout << "  rho= " << rho << endl;
		
		di=0.001;
		dd=0.001;
		dm=0.151;
		
		cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
		


		
			for (d=di; d<dm; d+=dd) {
	
				
			stringstream ss;
			ss << d;
		
			filename1="FrequentBrakeDis-v9-d"+ss.str()+"-l50k.txt";
				
			std::ofstream result(filename1.c_str());


			m=floor(d*N);
			
			Zeroint(nb,N);
			Zero(nbdis,1000);

			gap0=0;
			gap02=0;
			
				
			
			
			
			initializeuniform(N,d,x,v);
			
			
			
			
			for (int t=0; t<100000; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
			}  //initializing simulation
			
			
			
			
			
			for (int t=0; t<totaltime; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDownBreakDistribuion(N,d,vmax,countb,x,v,g,xb,nb);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
				
			for (i=0; i<m; i++) {
				nbdis[nb[i]]++;
			}	
				
				
				
			}  //main simulation
			
				for (i=0; i<100; i++) {
					
					result << i << "	" << nbdis[i]/(double(totaltime)*m) << endl;
				}
			
		//	result4.close();

			
			//result2.close();
			//result3.close();
		
			//result2<< d << "    ";

		
				
						
		
			
			result.close();

				
			cout <<"d = " << d << "  simulation done!"<<endl;
			
			
			
		}	
		
		
		
		cout << "Slow Down Probibility "<< p <<" is done!"<<endl;
	}
	
	return (0);
	
}




