
// This!
/*
 *  FullAnalyseBreak.cc
 *  
 *
 *  Created by Ashkan Balouchi on 1/28/14.
 *  Copyright 2014 LSU. All rights reserved.
 *
 */



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
	
	int n=5;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di;    //initial density
	double dm;    //max density 
	double dd;    //density step
	double p;     //probability of break
	int vmax;    //max velocity
	int m;     //number of cars
	int countb;
	int countf;
	double countbbar;
	double countfbar;
	int totaltime=20000000;        //number of totaltime steps
	int initime=1000000;     //number of initialization time steps
	int* car=new int[N];       //car numbers
	int* x=new int[N];       //positions
	int* xb=new int[N];       //breaking positions
	int* xf=new int[N];       //Free positions
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
	
	double* fwb=new double[N];
	double* fwb1=new double[N];
	
	double* ddcb=new double[N];
	double* ddcb1=new double[N];
	
	double* ddcb2=new double[N];
	
	double* ddcb3=new double[N];
	
	double* fwf=new double[N];
	double* ddcf=new double[N];
	
	
	
	double gap0,gap02;
	
	
	double vbar,vbar2,v2bar;    //averaged velocity
	int timestep=100;
	int k;
	double rho;
	
	int timef=5000;
	
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
	vmax=9;
	
	
	//File Names:
	//-------------------------
	string filename1,filename2,filename3,filename4,filename5,filename6,filename7,filename8,filename9,filename10,filename11,filename12,filename13,filename14,filename15,filename16;
	//-------------------------
	
	
	stringstream sss,ssl;
	sss << vmax;
	ssl << n;
	
	p=0.1;
	
	
	rho=(1-p)/(vmax+1-2*p);
	rho=round(rho,5);
	cout << " Transition density = " << rho << endl;
	
	di=0.08;
	dd=0.0002;
	dm=0.09;
	
	/*	di=0.08;
	 dd=0.0001;
	 dm=0.09;*/
	
	cout << "initial density = " << di<< ", Final density = "<<dm<<", density step: "<<dd<<endl;
	
	filename3="AveVel-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
	std::ofstream resultv(filename3.c_str());   //initialize file
	
	filename4="VelD-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
	std::ofstream resultvd(filename4.c_str());   //initialize file
	
	filename1="GapD-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
	std::ofstream result1(filename1.c_str());   //initialize file
	
	filename11="NumCarB-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
	std::ofstream resultcb(filename11.c_str());   //initialize file
	
	filename12="NumCarF-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
	std::ofstream resultcf(filename12.c_str());   //initialize file
	
	
	for (d=di; d<dm; d+=dd) {
		
		vbar=0;
		vbar2=0;
		v2bar=0;
		stringstream ss1;
		ss1 << d;
		
		
		
		m=floor(d*N);
		//cout << m <<endl;
		gap0=0;
		gap02=0;
		
		filename2="GapDFull-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultgdt(filename2.c_str());   //initialize file
		
		filename5="Fourier1D-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultf(filename5.c_str());   //initialize file
		
		filename6="DDC-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddc(filename6.c_str());   //initialize file
		
		filename7="Fourier1DBreak-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultfb(filename7.c_str());   //initialize file
		
		filename16="Fourier1DBreak1-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultfb1(filename16.c_str());   //initialize file
		
		filename8="DDCB-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddcb(filename8.c_str());   //initialize file
		
		filename13="DDCB1-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddcb1(filename13.c_str());   //initialize file
		
		filename14="DDCB2-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddcb2(filename14.c_str());   //initialize file
		
		filename15="DDCB3-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddcb3(filename15.c_str());   //initialize file
		
		filename9="Fourier1DFree-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultff(filename9.c_str());   //initialize file
		
		filename10="DDCF-d"+ss1.str()+"-l"+ssl.str()+"k-v"+sss.str()+"-p0.1"+".txt";
		std::ofstream resultddcf(filename10.c_str());   //initialize file
		
		
		initializeuniform(N,d,x,v);
		
		for (int t=0; t<initime; t++) {  
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDown(N,d,vmax,x,v,g);
			Randomize(N,d,p,v);
			MoveCars(N,d,vmax,x,v,g);
		}  //initializing simulation
		
		
		for (int jj=1; jj<N; jj++) {
			auxgap[jj]=0;
		}
		for (int jj=0; jj<vmax+1; jj++) {
			auxv[jj]=0;
		}
		
		Zero(fw,N);
		Zero(fwb,N);
		Zero(fwb1,N);
		
		Zero(ddcb,N);
		Zero(ddcb1,N);
		Zero(ddcb2,N);
		Zero(ddcb3,N);
		
		Zero(fwf,N);
		Zero(ddcf,N);
		
		countbbar=0;
		countfbar=0;
		
		
		for (int t=0; t<totaltime-timef*timestep; t++) {  //main simulation
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDownBreakJam(N,d,vmax,countb,countf,x,v,g,xb,xf,car);
			Randomize(N,d,p,v);
			MoveCars(N,d,vmax,x,v,g);
			
			countbbar+=countb;
			countfbar+=countf;
			
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
		
		for (int t=totaltime-timef*timestep; t<totaltime; t++) {  //main simulation
			
			
			
			Accelerate(N,d,vmax,x,v,g);
			SlowDownBreakJam(N,d,vmax,countb,countf,x,v,g,xb,xf,car);
			Randomize(N,d,p,v);
			MoveCars(N,d,vmax,x,v,g);
			
			/*
			 for (int kk=0; kk<countb; kk++) {
			 cout << car[kk]<< "   ";
			 }
			 cout << endl;*/
			
			countbbar+=countb;
			countfbar+=countf;
			//cout << countb << endl;
			
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
				
				if (countb>0) {
					
					
					
					Zeroint(bit,N);
					for (int i=0; i<countb; i++) {	bit[xb[i]]=1;	} //finds site bits
					for (int i=0; i<N; i++) {
						in[i][0]=bit[i];			
						in[i][1]=0;			
					}
					fftw_execute(p1); 
					for (int i=0; i<N; i++) {
						fwb[i]+=out[i][0]*out[i][0]+out[i][1]*out[i][1];
						fwb1[i]+=countb*(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
						
					}
					
					
					//in2[0][0]=countb-1;
					in2[0][0]=countb-1;
					in2[0][1]=0;
					for (int i=1; i<N; i++) {
						in2[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])/(double)countb-1;				
						in2[i][1]=0;			
					}
					
					fftw_execute(p2); 
					
					
					for (int i=0; i<N; i++) {
						ddcb[i]+=out2[i][0];
						ddcb1[i]+=countb*out2[i][0];
						
					}
					
					in2[0][0]=countb;
					in2[0][0]=0;
					
					for (int i=1; i<N; i++) {
						in2[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])/(double)countb-1;				
						in2[i][1]=0;			
					}
					
					fftw_execute(p2); 
					
					
					for (int i=0; i<N; i++) {
						ddcb2[i]+=out2[i][0];
						ddcb3[i]+=countb*out2[i][0];
						
					}
					
				}
				
				Zeroint(bit,N);
				for (int i=0; i<countf; i++) {	bit[xf[i]]=1;	} //finds site bits
				for (int i=0; i<N; i++) {
					in[i][0]=bit[i];			
					in[i][1]=0;			
				}
				fftw_execute(p1); 
				for (int i=0; i<N; i++) {
					fwf[i]+=out[i][0]*out[i][0]+out[i][1]*out[i][1];
				}
				
				
				in2[0][0]=countf-1;
				in2[0][1]=0;
				for (int i=1; i<N; i++) {
					in2[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])/(double)countf-1;				
					in2[i][1]=0;			
				}
				
				fftw_execute(p2); 
				
				
				for (int i=0; i<N; i++) {
					ddcf[i]+=out2[i][0];
				}
				
				
				
			}
			
		}			
		
		
		
		
		
		for (int i=1; i<N; i++) {
			resultf << i <<  "   " << fw[i]/double(timef) << std::endl;	
			resultfb << i <<  "   " << fwb[i]/double(timef) << std::endl;
			resultfb << i <<  "   " << fwb[i]/(double(timef)*countbbar) << std::endl;
			
			resultff << i <<  "   " << fwf[i]/double(timef) << std::endl;	
		}
		
		
		in2[0][0]=m-1;
		in2[0][1]=0;
		for (int i=1; i<N; i++) {
			in2[i][0]=fw[i]/(double(timef)*(double)m)-1;				
			in2[i][1]=0;			
		}
		
		fftw_execute(p2); 
		
		cout << fw[0]/time << std::endl;
		cout <<"d = " << d << "  1d fourier done"<<std::endl;
		
		resultf.close();
		resultfb.close();
		resultfb1.close();
		resultff.close();
		
		countbbar/=totaltime;
		countfbar/=totaltime;
		
		
		for (int i=0; i<N; i++) {
			resultddc << i <<  "   " << out2[i][0]/double(N)<< std::endl;	
			resultddcb << i <<  "   " << ddcb[i]/(double(timef)*double(N))<< std::endl;		
			resultddcb1 << i <<  "   " << ddcb1[i]/(countbbar*double(N))<< std::endl;
			resultddcb2 << i <<  "   " << ddcb2[i]/(double(timef)*double(N))<< std::endl;		
			resultddcb3 << i <<  "   " << ddcb3[i]/(countbbar*double(N))<< std::endl;
			resultddcf << i <<  "   " << ddcf[i]/(double(timef)*double(N))<< std::endl;		
			
		}
		
		countbbar/=totaltime;
		countfbar/=totaltime;	
		
		resultddc.close();
		resultddcb.close();
		resultddcb1.close();
		resultddcb2.close();
		resultddcb3.close();
		resultddcf.close();
		
		cout << "D.D.Correlation done!"<<endl;
		
		
		result1<< d << "    ";
		resultvd<< d << "    ";
		
		
		
		for (int jj=1; jj<vmax+2; jj++) {
			result1 << auxgap[jj]/(m*double(time)) << "    ";
			resultvd << auxv[jj-1]/(m*double(time)) << "    ";
			
		}
		result1 << endl;
		resultvd << endl;
		
		
		for (int jj=1; jj<N; jj++) {
			resultgdt << jj<< "    "<<auxgap[jj]/(m*double(time)) <<endl;
		}
		
		
		vbar /=(double)time;
		v2bar /=(double)time;
		vbar2 /=(double)time;
		resultv << d << "\t" << vbar << "\t" << vbar2 << "\t" << v2bar <<endl;										  
		resultcb << d << "\t" << countbbar << endl;
		resultcf << d << "\t" << countfbar << endl;
		
		
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
	resultcb.close();
	resultcf.close();
	
	
	
	
	return (0);
	
	

	
}


















