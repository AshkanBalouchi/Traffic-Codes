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

//#include <map>
//#include <random>
//#include <fftw3.h>

using namespace std;



int main(){
	
	
	int n=100;
	int num=10;
	
	
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	
	di=0.088;
	dd=0.0002;
	dm=0.093;
	
	
	
	double p;     //probability of break
	int vmax=8;    //max velocity
	int m;     //number of cars
	int totaltime=200000;        //number of totaltime steps
	int initiationtime=1000000;  
	int* x=new int[N/5]; 
	//int* xb=new int[N]; 
	int* v=new int[N/5]; 
	int* g=new int[N/5]; 
	//double* gapd=new double[N]; 
	
	//int* nb=new int[N]; 
	//double* gapdf=new double[N]; 
	double gap0,gap02;
	double vbar,vbar2,v2bar,vvar;    //averaged velocity
	int timestep=10;  //time step
	int averageingtime=100; 
	int tmod,tnum;
	
	int vel[averageingtime*N/5]; 
	double* vav=new double[N/5]; 


	int count,countb;
	int i,j,k;
	double rho;
	double X41,X42;
	
	
	
	int time=floor(totaltime/timestep);
	
	stringstream sss,runnum;
	
	
	//int num;
	//num=3;
	runnum << num;
	
	
	sss << vmax;
	
	for (p=0.1; p<0.2; p+=0.1) {

		gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);			

		
		
		string filename1,filename2,filename3,filename4;
		rho=(1-p)/(vmax+1-2*p);
		rho=round(rho,3);
		cout << "  rho= " << rho << endl;
	
		
		/*di=0.079;
		dd=0.0001;
		dm=0.095;*/
		
		cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
		

		stringstream ss,sss;
		ss << n;
		sss << vmax;
		
		/*filename1="v9-l"+ss.str()+"k-susceptibility-test-1.txt";
		filename2="v9-l"+ss.str()+"k-susceptibility-test-2.txt";
		filename3="v9-l"+ss.str()+"k-susceptibility-test-3.txt";*/
		
		filename3="v"+sss.str()+"-l"+ss.str()+"k-susceptibility-time-normalized-"+runnum.str()+".txt";

		filename4="v"+sss.str()+"-l"+ss.str()+"k-susceptibility-time-top-"+runnum.str()+".txt";



		
	//	std::ofstream result1(filename1.c_str());
	//	std::ofstream result2(filename2.c_str());
		std::ofstream result3(filename3.c_str());
		std::ofstream result4(filename4.c_str());


		
			for (d=di; d<dm; d+=dd) {

					
			m=floor(d*N);
			
			
			
			initializeuniform(N,d,x,v);
			
			
				srand (num*1000000);

			
			for (int t=0; t<initiationtime; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v,r);
				MoveCars(N,d,vmax,x,v,g);
			}  //initializing simulation
			

				
				for (int t=0; t<averageingtime*timestep; t++) {  //main simulation
					
					Accelerate(N,d,vmax,x,v,g);
					SlowDown(N,d,vmax,x,v,g);
					Randomize(N,d,p,v,r);
					MoveCars(N,d,vmax,x,v,g);
					
					if (t % timestep == 0) {
						
						tmod = t/timestep;
						
						for (i=0; i<m; i++) {
						
						vel[tmod+averageingtime*i]=v[i];
							
						}

					}
					
				} 
	
			
			
				X41=0;
				X42=0;
				
				for (j=0; j<m; j++) {
					vav[j]=0;
					for (i=0; i<averageingtime; i++) {
						vav[j]+=vel[i+averageingtime*j];
				  }
					vav[j] = vav[j]/averageingtime;
				}
		
				
				
			for (int t=0; t<totaltime; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v,r);
				MoveCars(N,d,vmax,x,v,g);
				
				
				
				
				if (t % timestep == 0) {
			
					tnum=t/timestep;
					tmod=tnum % averageingtime;
					
					
					
					
					double aux3=0;
					double aux1=0;
					double aux2=0;
					double aux4=0;
					
					
					for (j=0; j<m; j++) {
						vav[j]=(vav[j]*averageingtime - vel[tmod+averageingtime*j]+v[j])/averageingtime;
						vel[tmod+averageingtime*j]=v[j];
					}
					
					
					vbar=0;
					vvar=0;
					v2bar=0;
					
					for (j=0; j<m; j++) {
						vbar+=v[j];
					}
					
					vbar = vbar/double(m);
					
					for (j=0; j<m; j++) {
						vvar += (v[j]-vbar)*(v[j]-vbar);
						for (i=0; i<m; i++) {
							aux3+=(v[j]-vav[j])*(v[i]-vav[i]);
							
						}
					}
					
				//	vvar=vvar/double(m);
					
					
				/*	for (j=0; j<m; j++) {
						aux4+=(aux1/double(m)-v[j])*(aux1/double(m)-v[j]);
					}
					
					X41+=((aux3/double(m)-aux2)-aux4)/(aux2-aux1*aux1/double(m));*/
					X41+=aux3/vvar;
					X42+=aux3;

				
					
					
					
					
					
				}
				
			}  //main simulation
			
			
			
	
				
				
				
				
				
				/*result1<< d << "	" << X41/(double(time)) << endl;
				result2<< d << "	" << X41/(double(time)*double(m)) << endl;
				result3<< d << "	" << X42/(double(time)) << endl;*/
				result3<< d << "	" << X41/(double(time)) << endl;

				result4<< d << "	" << X42/(double(time)*double(m)) << endl;


		
			



				
			cout <<"v = " << vmax << " d = " << d << "  simulation done!"<<endl;
			
			
			
		}	
		
		//result1.close();
		//result2.close();
		result3.close();
		result4.close();


		
		
		cout << "Slow Down Probibility "<< p <<" is done!"<<endl;
	}
	
	return (0);
	
}




