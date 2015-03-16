/*
 *  FullAnalyse.cc
 *  
 *
 *  Created by Ashkan Balouchi on 1/14/14.
 *  Copyright 2014 LSU. All rights reserved.
 *
 */

//#include "FullAnalyse.h"


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
	
	int n=2;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	int m;     //number of cars
	int totaltime=10000000;        //number of totaltime steps
	int initime=1000000;     //number of initialization time steps
	int* x=new int[N];       //positions
	int* v=new int[N];       //velocities
	double* vd=new double[N];    //gaps aux
	int* g=new int[N];       //gaps
	double* gapd=new double[N];    //gaps aux
	double auxgap[N];
	double auxv[N];
	
	double* dden=new double[N]; 
	double Gaptotal;
	int* bit=new int[N]; 
	double* fw=new double[N];
	
	
	double gap0,gap02;

	
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=100;
	int count;
	int k;
	double rho;
	
	
	
	//Fourier needs:
	
	int time=floor(totaltime/timestep);
	int* fourierbit=new int[N*time]; 

	fftw_complex *in,*out,*in2,*out2;
	fftw_plan p1,p2;
	
	in = (fftw_complex*) fftw_malloc(N* sizeof(fftw_complex));
	out = (fftw_complex*) fftw_malloc(N* sizeof(fftw_complex));
	in2 = (fftw_complex*) fftw_malloc(N* sizeof(fftw_complex));
	out2 = (fftw_complex*) fftw_malloc(N* sizeof(fftw_complex));
	p1 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//cout << "vmax = ";
	//cin	>> vmax;
	vmax=8;
	
	
	//File Names:
	//-------------------------
	string filename1,filename2,filename3,filename4,filename5,filename6;
	//-------------------------
	
	
	stringstream sss,ssl;
	sss << vmax;
	ssl << n;
	
	p=0.1;
	
	
	rho=(1-p)/(vmax+1-2*p);
	rho=round(rho,5);
	cout << " Transition density = " << rho << endl;
	
	
	
	di=0.085;
	dd=0.0002;
	dm=0.1;
	
	cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
	
	filename3="AveVel-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-low"+".txt";
	std::ofstream resultv(filename3.c_str());   //initialize file

	filename4="VelD-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-low"+".txt";
	std::ofstream resultvd(filename4.c_str());   //initialize file
	
	filename1="GapD-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-low"+".txt";
	std::ofstream result1(filename1.c_str());   //initialize file

	
	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);			

	
	for (d=di; d<dm; d+=dd) {
	
		vbar=0;
		vbar2=0;
		v2bar=0;
		stringstream ss1;
		ss1 << d;
		
		
		
		m=floor(d*N);
		gap0=0;
		gap02=0;
		
		filename2="GapDFull-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultgdt(filename2.c_str());   //initialize file
		
		filename5="Fourier1D-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultf(filename5.c_str());   //initialize file
		
		filename6="DDC-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddc(filename6.c_str());   //initialize file
		
		initializeuniform(N,d,x,v);
		
		for (int t=0; t<initime; t++) {  
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v,r);
			MoveCars(N,d,vmax,x,v,g);
		}  //initializing simulation
		
		
		for (int jj=1; jj<N; jj++) {
			auxgap[jj]=0;
		}
		for (int jj=0; jj<vmax+1; jj++) {
			auxv[jj]=0;
		}
		
		Zero(fw,N);

		
		for (int t=0; t<totaltime-2000*timestep; t++) {  //main simulation
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v,r);
			MoveCars(N,d,vmax,x,v,g);
		
			
			
			if (t % timestep == 0) {
				
				FindGapDistribution(N,d,g,gapd);
				FindVelDistribution(N,d,vmax,v,vd);

				
				vbar += average(v,m);
				v2bar += average2(v,m);
				vbar2 += average(v,m)*average(v,m);
				
				for (int gg=1; gg<N; gg++) {
					auxgap[gg]+=gapd[gg];
				}
				for (int gg=0; gg<vmax+1; gg++) {
					auxv[gg]+=vd[gg];
				}
			
			}
		
		}			
	
		for (int t=totaltime-2000*timestep; t<totaltime; t++) {  //main simulation
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v,r);
			MoveCars(N,d,vmax,x,v,g);
			
			
			
			if (t % timestep == 0) {
				
				FindGapDistribution(N,d,g,gapd);
				FindVelDistribution(N,d,vmax,v,vd);
				
				
				vbar += average(v,m);
				v2bar += average2(v,m);
				vbar2 += average(v,m)*average(v,m);
				
				for (int gg=1; gg<N; gg++) {
					auxgap[gg]+=gapd[gg];
				}
				for (int gg=0; gg<vmax+1; gg++) {
					auxv[gg]+=vd[gg];
				}
				
			
				Zeroint(bit,N);
				for (int i=0; i<m; i++) {	bit[x[i]]=1;	} //finds site bits
				for (int i=0; i<N; i++) {
					in[i][0]=bit[i];			
					in[i][1]=0;			
				}
				fftw_execute(p1); 
				for (int i=0; i<N; i++) {
					fw[i]+=out[i][0]*out[i][0]+out[i][1]*out[i][1];
				}
			
			}
			
		}			
		
		
		for (int i=1; i<N; i++) {
			resultf << i <<  "   " << fw[i]/2000. << std::endl;		
		}
		
		cout << fw[0]/time << std::endl;
		cout <<"d = " << d << "  1d fourier done"<<std::endl;
		
		resultf.close();
		
		in2[0][0]=m-1;
		in2[0][1]=0;
		for (int i=1; i<N; i++) {
			in2[i][0]=fw[i]/(2000.*(double)m)-1;				
			in2[i][1]=0;			
		}
		
		fftw_execute(p2); 
		
		for (int i=0; i<N; i++) {
			resultddc << i <<  "   " << out2[i][0]/double(N)<< std::endl;		
		}
		
		
		resultddc.close();
		cout << "D.D.Correlation done!"<<endl;

		
		result1<< d << "    ";
		resultvd<< d << "    ";

		
		for (int jj=1; jj<vmax+2; jj++) {
			result1 << auxgap[jj]/(m*double(time)) << "    ";
			resultvd << auxv[jj-1]/(m*double(time)) << "    ";

		}
		result1 << endl;
		resultvd << endl;

		
		for (int jj=1; jj<m+1; jj++) {
			resultgdt << jj<< "    "<<auxgap[jj]/(m*double(time)) <<endl;
		}
		
												  
		vbar /=(double)time;
		v2bar /=(double)time;
		vbar2 /=(double)time;
		resultv << d << "\t" << vbar << "\t" << vbar2 << "\t" << v2bar <<endl;										  
												  
		
		resultgdt.close();
		cout <<"d = " << d << "  simulation done!"<<endl;

		
		
		
		
		
		
		
	
	}
	
	
	
	fftw_destroy_plan(p1);
	
	fftw_free(in); fftw_free(out);	
	
	fftw_destroy_plan(p2);
	
	fftw_free(in2); fftw_free(out2);
	
	result1.close();
    resultv.close();
    resultvd.close();

	
	return (0);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	





}










