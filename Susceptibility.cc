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
	
	
	int n=20;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	int m;     //number of cars
	int totaltime=800000;        //number of totaltime steps
	int* x=new int[N]; 
	int* xb=new int[N]; 
	int* v=new int[N]; 
	int* g=new int[N]; 
	double* gapd=new double[N]; 
	
	int* nb=new int[N]; 
	//double* gapdf=new double[N]; 
	double gap0,gap02;
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=10;
	int count,countb;
	int i,j,k;
	double rho;
	double X41,X42;
	
	
	
	int time=floor(totaltime/timestep);
	
	stringstream sss;
	sss << vmax;
	
	for (p=0.1; p<0.2; p+=0.1) {

		
		string filename1,filename2,filename3,filename4;
		vmax=9;
		rho=(1-p)/(vmax+1-2*p);
		rho=round(rho,3);
		cout << "  rho= " << rho << endl;
	
		
		di=0.005;
		dd=0.001;
		dm=0.15;
		
		cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
		

		stringstream ss;
		ss << n;
		
		/*filename1="v9-l"+ss.str()+"k-susceptibility-test-1.txt";
		filename2="v9-l"+ss.str()+"k-susceptibility-test-2.txt";
		filename3="v9-l"+ss.str()+"k-susceptibility-test-3.txt";*/
		
		filename3="v5-l"+ss.str()+"k-susceptibility-test-1.txt";

		filename4="v5-l"+ss.str()+"k-susceptibility-test-2.txt";



		
	//	std::ofstream result1(filename1.c_str());
	//	std::ofstream result2(filename2.c_str());
		std::ofstream result3(filename3.c_str());
		std::ofstream result4(filename4.c_str());


		
			for (d=di; d<dm; d+=dd) {

					
			m=floor(d*N);
			
			
			
			initializeuniform(N,d,x,v);
			
			
			
			
			for (int t=0; t<1000000; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
			}  //initializing simulation
			
			
	
			
			
				X41=0;
				X42=0;
				
			for (int t=0; t<totaltime; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
				
				
				
				
				if (t % timestep == 0) {
			
					
					
					int ll=0;
					double aux3=0;
					double aux1=0;
					double aux2=0;
					double aux4=0;
					for (j=0; j<m; j++) {
						aux1+=v[j];
						aux2+=v[j]*v[j];
						for (i=0; i<m; i++) {
						//	ll = ((i+j) % m);
							aux3+=v[j]*v[i];
							
						}
					}
				/*	for (j=0; j<m; j++) {
						aux4+=(aux1/double(m)-v[j])*(aux1/double(m)-v[j]);
					}
					
					X41+=((aux3/double(m)-aux2)-aux4)/(aux2-aux1*aux1/double(m));*/
					X42+=aux3/aux2;
					X41+=aux2-aux3/m;

					
					
				}
				
			}  //main simulation
			
			
			
	
				
				
				
				
				
				/*result1<< d << "	" << X41/(double(time)) << endl;
				result2<< d << "	" << X41/(double(time)*double(m)) << endl;
				result3<< d << "	" << X42/(double(time)) << endl;*/
				result3<< d << "	" << X41/(double(time)) << endl;

				result4<< d << "	" << X42/(double(time)*double(m)) << endl;


		
			



				
			cout <<"d = " << d << "  simulation done!"<<endl;
			
			
			
		}	
		
		//result1.close();
		//result2.close();
		result3.close();
		result4.close();


		
		
		cout << "Slow Down Probibility "<< p <<" is done!"<<endl;
	}
	
	return (0);
	
}




