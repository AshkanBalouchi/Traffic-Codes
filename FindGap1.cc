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
	
	
	int N=100000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	double m;     //number of cars
	int totaltime=1000000;        //number of totaltime steps
	int* x=new int[N]; 
	int* xb=new int[N]; 
	int* v=new int[N]; 
	int* g=new int[N]; 
	double* gapd=new double[N]; 
	
	int* nb=new int[N]; 
	//double* gapdf=new double[N]; 
	double gap0,gap02;
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=1;
	int count,countb;
	int i,j,k;
	double rho;
	
	
	
	int time=floor(totaltime/timestep);
	//cout << "vmax = ";
	//cin	>> vmax;
	//vmax=9;
	/*cout << "p = ";
	 cin >> p ;*/
	/*cout << " initial d = ";
	 cin >> di;
	 cout << " density step = ";
	 cin >> dd;
	 cout << " max density = ";
	 cin >> dm;
	 cout << std::endl;*/
	stringstream sss;
	sss << vmax;
	
	for (p=0.1; p<0.2; p+=0.1) {
		
		
		//string address="/Users/ashkanbalouchi/Desktop/Traffic/TrafficSim/FFT-Analyse/bit/results/vmax"+sss.str()+"/gaptest/";
		//string name1="bit.d";
		//string name2=".p";
		//string name3=".v";
		//string name7=".txt";
		string filename1,filename2,filename3,filename4;
		//stringstream ss2,ss3;
		//ss2 << p;
		//ss3 << vmax;
		//string name8="GapD.d.";
		//string name9="Gap2D.d.";
		//string name10="Gap3D.d.";
		
		
		
//	filename2="/Users/ashkanbalouchi/Desktop/Traffic/TrafficSim/FFT-Analyse/bit/results/vmax11/gaptest/l4e4-all.txt";
	//	filename3="/Users/ashkanbalouchi/Desktop/Traffic/TrafficSim/FFT-Analyse/bit/results/vmax11/gaptest/l4e4-stdv-all.txt";
	//	filename4="/Users/ashkanbalouchi/Desktop/Traffic/TrafficSim/FFT-Analyse/bit/results/vmax11/gap3/timeevolution-l5e4-test2.txt";
		
		vmax=9;
		rho=(1-p)/(vmax+1-2*p);
		rho=round(rho,3);
		cout << "  rho= " << rho << endl;
		/* di=rho*0.1;
		 di=round(di,3);
		 dm=rho*1.2;
		 dm=round(dm,3);
		 dd=(dm-di)/22.0;
		 dd=round(dd,3);*/
		
		//di=2*rho;
		//dd=rho;
		//dm=5*rho;
		
		di=0.001;
		dd=0.001;
		dm=0.151;
		
		cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
		
	//	std::ofstream result2(filename2.c_str());
	//	std::ofstream result3(filename3.c_str());
//		std::ofstream result4(filename4.c_str());
		
		
		filename1="v9-l100k-corr.txt";
		
		std::ofstream result4(filename1.c_str());
		
			for (d=di; d<dm; d+=dd) {
	//	for (d=0.067; d<0.07301; d+=0.0002) {

		//for (d=0.08; d<0.088; d+=0.0001) {
				
				stringstream ss;
				ss << d;
				
				filename2="v9-d"+ss.str()+"-l100k.txt";
				
				std::ofstream result2(filename2.c_str());

			
			vbar=0;
			vbar2=0;
			v2bar=0;
			stringstream ss1;
			ss1 << d;
			
			
			
		//	filename4=address+"timeevolution.d."+ss1.str()+".l4e4-all"+name7;
			//filename1=address+name8+ss1.str()+name2+ss2.str()+name3+ss3.str()+"e5-e6-2e6-e2"+name7;
			//filename2=address+name9+ss1.str()+name2+ss2.str()+name3+ss3.str()+name7;
			//filename3=address+name10+ss1.str()+name2+ss2.str()+name3+ss3.str()+name7;
			
			/*	filename1=address+name8+ss1.str()+name2+ss2.str()+name3+ss3.str()+".random"+name7;
			 filename2=address+name9+ss1.str()+name2+ss2.str()+name3+ss3.str()+".random"+name7;
			 filename3=address+name10+ss1.str()+name2+ss2.str()+name3+ss3.str()+".random"+name7;*/
			
			//std::ofstream result1(filename1.c_str());  
			//std::ofstream result2(filename2.c_str());
			//std::ofstream result3(filename3.c_str());   //initialize file
			
		//	std::ofstream result4(filename4.c_str());

			
			m=floor(d*N);
			
			//	Zero(gapdf,N);
			gap0=0;
			gap02=0;
			//	Zero(gapdf2,N);
			//	Zero(gapdf3,N);
			
			
			
			initializeuniform(N,d,x,v);
			
			
			
			
			for (int t=0; t<100000; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
			}  //initializing simulation
			
			//	FindGapDistribution(N,d,g,gapd);
			
			
			double auxgap[N];
			for (int jj=0; jj<N; jj++) {
				auxgap[jj]=0;
			}
			
			
			for (int t=0; t<totaltime; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDownBreakDistribuion(N,d,vmax,countb,x,v,g,xb,nb);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
				
				
				
				
				if (t % timestep == 0) {
					FindGapDistribution(N,d,g,gapd);
					//	Find2GapDistribution(N,d,g,gapd2);
					//	Find3GapDistribution(N,d,g,gapd3);
					
										
				//	double aux2=0;
					for (int gg=0; gg<N; gg++) {
						auxgap[gg]+=gapd[gg];
		//				aux2+=gapd[gg]*gapd[gg];
					}
					
				//	gap0+=aux/(m);
				//	gap02+=aux2/(m*m);
					
				//	result4 << t/timestep << "    "<< aux/(m) <<endl;
					
				}
				
			}  //main simulation
			
			
			
		//	result4.close();

			
			//result2.close();
			//result3.close();
		
			//result2<< d << "    ";

			
				for (int jj=0; jj<N; jj++) {
					auxgap[jj]=auxgap[jj]/(m*double(totaltime)/timestep);
					result2 << jj << "	" << auxgap[jj] << endl;
			}
			
			result2 << endl;
				
				int ll=0;
				double aux=0;
				double aux2=0;
				for (j=0; j<N; j++) {
					aux2+=auxgap[j]*auxgap[j];
					for (i=0; i<N; i++) {
						ll = ((i+j) % N);
						aux+=auxgap[j]*auxgap[ll];

					}
				}
				result4<< d << "	" << aux/aux2 << endl;
			
			//result2 << d << "    "<< gap0/(double(totaltime)/timestep)<<endl;
			//result3 << d<<  "  " << sqrt(gap02/(double(totaltime)/timestep)-gap0/(double(totaltime)/timestep)*gap0/(double(totaltime)/timestep)) << endl;
			
			//result3 << d<<  "  " << gap02/(double(totaltime)/timestep) << endl;
			
			
			
			
			
			/*		 initializeRandom2(N,d,x,v);
			 
			 FindGaps(N,d,vmax,x,v,g);
			 
			 FindGapDistribution(N,d,g,gapd);
			 Find2GapDistribution(N,d,g,gapd2);
			 Find3GapDistribution(N,d,g,gapd3);
			 
			 double Gaptotal=0;
			 
			 double Gaptotal2=0;
			 
			 double Gaptotal3=0;
			 
			 for (int j=0; j<N; j++) { 
			 Gaptotal+=gapd[j]; 
			 Gaptotal2+=gapd2[j];
			 Gaptotal3+=gapd3[j];
			 
			 }
			 
			 for (int j=0; j<20/d; j++) {
			 result1 << j<<  "  " << gapd[j]/double(Gaptotal) << endl;
			 result2 << j<<  "  " << gapd2[j]/double(Gaptotal2) << endl;
			 result3 << j<<  "  " << gapd3[j]/double(Gaptotal3) << endl;
			 }
			 
			 
			 
			 result1.close();
			 result2.close();
			 result3.close();*/
			
				result2.close();

				
			cout <<"d = " << d << "  simulation done!"<<endl;
			
			
			
		}	
		
//		result2.close();
		//result3.close();
		result4.close();
		
		
		cout << "Slow Down Probibility "<< p <<" is done!"<<endl;
	}
	
	return (0);
	
}




