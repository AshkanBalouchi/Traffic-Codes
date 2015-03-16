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
	
	
	int n=5;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	int m;     //number of cars
	int totaltime=5000000;        //number of totaltime steps
	int* x=new int[N]; 
	//int* xb=new int[N]; 
	int* v=new int[N]; 
	//int* vr=new int[N]; 
	//int* vs=new int[N]; 
	//int* blr=new int[N]; 
	//int* bls=new int[N]; 

	//int* blrdis=new int[N]; 
	//int* blsdis=new int[N]; 

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
	double NNVCorr;
	double NNVV;

	
	double aux1;
	double aux2;
	double aux3;
	double aux4;
	double aux5;
	
	
	int time=floor(totaltime/timestep);
	
	vmax=9;

	stringstream sss;
	sss << vmax;
	
	for (p=0.1; p<0.2; p+=0.1) {

		
		string filename1,filename2,filename3,filename4;
		rho=(1-p)/(vmax+1-2*p);
		rho=round(rho,3);
		cout << "  rho= " << rho << endl;
	
		
		di=0.06;
		dd=0.0005;
		dm=0.09;
		
		cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
		

		stringstream ss;
		ss << n;
		
		filename3="VCorr-v"+sss.str()+"-l"+ss.str()+"k-fine.txt";
	//	filename4="v"+sss.str()+"-l"+ss.str()+"k-NNVVariance2.txt";


		std::ofstream result3(filename3.c_str());
	//	std::ofstream result4(filename4.c_str());



		
			for (d=di; d<dm; d+=dd) {

			
				
				aux1=0;
				aux2=0;
				aux3=0;
				aux4=0;
				aux5=0;
				
				
				//	NNVCorr=0;
			//	NNVV=0;

				
			m=floor(d*N);
			
			//	stringstream ss1;
			//	ss1 << d;

			/*	filename1="v"+sss.str()+"-l"+ss.str()+"k-d"+ss1.str()+"-brakelightgap.txt";
				filename2="v"+sss.str()+"-l"+ss.str()+"k-d"+ss1.str()+"-brakelighton.txt";


				std::ofstream result1(filename1.c_str());
				std::ofstream result2(filename2.c_str());*/


			
			
			initializeuniform(N,d,x,v);
			
			
			
			for (int t=0; t<1000000; t++) {  //main simulation
				
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
			}  //initializing simulation
			
			
	
			
			
				X41=0;
				X42=0;
				
			//	Zeroint(blrdis,N);
			//	Zeroint(blsdis,N);

				
			for (int t=0; t<totaltime; t++) {  //main simulation
				
			/*	for (j=0; j<m; j++) {
					vr[j]=v[j];
				}*/
				Accelerate(N,d,vmax,x,v,g);
				SlowDown(N,d,vmax,x,v,g);
			//	SlowDownBrakeLight(N,d,vmax,x,v,g,bls);
				Randomize(N,d,p,v);
				MoveCars(N,d,vmax,x,v,g);
				
			/*	for (i=0; i<m; i++) {
					if (v[i] < vr[i]){
						blr[i]=1;
					}else {
						blr[i]=0;
					}

				//	cout << v[i]<<"  " << vr[i]<< "  "<<blr[i]<<endl;
				}*/
				
				
				
				
			//	int l=0;
			//	double aux11=0;
			//	double aux1=0;
			//	double aux22=0;
			//	double aux2=0;
			/*	for (i=1; i<20; i++) {
					
					aux2=0;
					aux22=0;

					for (j=0; j<m; j++) {
						
						aux1=1;
						aux11=1;

						for (k=0; k<i; k++) {
							l = (j+k) % m;
							aux1=aux1*blr[l];
							aux11=aux11*bls[l];	

						}
						
						aux2+=aux1;
						aux22+=aux11;

					}
					
					blrdis[i]+=aux2;
					blsdis[i]+=aux22;

					
				}*/
				
				
				
				
				if (t % timestep == 0) {
			
					
					 

				
					for (j=0; j<m-3; j++) {
						aux1+=v[j];
						aux2+=v[j]*v[j];
					    aux3+=v[j]*v[j+1];
						aux4+=v[j]*v[j+2];
					    aux5+=v[j]*v[j+3];

						}
					aux1+=v[m-3]+v[m-2]+v[m-1];
					aux2+=v[m-3]*v[m-3]+v[m-2]*v[m-2]+v[m-1]*v[m-1];
					aux3+=v[m-3]*v[m-2]+v[m-2]*v[m-1]+v[m-1]*v[0];
					aux4+=v[m-3]*v[m-1]+v[m-2]*v[0]+v[m-1]*v[1];
					aux5+=v[m-3]*v[0]+v[m-2]*v[1]+v[m-1]*v[2];

				
				//	NNVCorr+=(aux3/m-(aux1*aux1/(m*m)))/(aux2/m-(aux1*aux1/(m*m)));
				//	NNVV+=(aux3/m-(aux1*aux1/(m*m)));


					
					
				}
				
			}  //main simulation
			
			
			
	
				
				result3<< d << "	" << std::setprecision(10) << aux1/(m*double(time)) << "		" << std::setprecision(10)<< aux2/(m*double(time)) << "		" << std::setprecision(10)<< aux3/(m*double(time)) << "		" << std::setprecision(10)<< aux4/(m*double(time)) << "		" << std::setprecision(10)<< aux5/(m*double(time)) << endl;


				
				
				
			/*	result1<< d << "	" << X41/(double(time)) << endl;
				result2<< d << "	" << X41/(double(time)*double(m)) << endl;
				result3<< d << "	" << X42/(double(time)) << endl;
				result4<< d << "	" << X42/(double(time)*double(m)) << endl;

*/
		
			


		/*		double auxx1,auxx2;
			for (i=1; i<20; i++) {
				auxx1=blsdis[i]/(double(totaltime)*double(m));
				auxx2=blrdis[i]/(double(totaltime)*double(m));
				result1<< i << "	" << auxx1 << endl;
				result2<< i << "	" << auxx2 << endl;

			}	*/
			
				
				
			//	result1.close();
			//	result2.close();
			
				
				cout <<" l = " << n << " v = "<< vmax << "  d = " << d << "  simulation done!"<<endl;

			
		}	
		
		result3.close();
	//	result4.close();


		
		
		cout << "Slow Down Probibility "<< p <<" is done!"<<endl;
	}
	
	return (0);
	
}




