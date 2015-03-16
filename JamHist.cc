/*
 *  JamHist.cc
 *  
 *
 *  Created by Ashkan Balouchi on 1/10/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */


/*
 *  JamTimeEvo.cc
 *  
 *
 *  Created by Ashkan Balouchi on 1/7/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */


/*
 *  VelGap.cpp
 *  
 *
 *  Created by Ashkan Balouchi on 1/5/15.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */


/*
 *  GapDis.cc
 *  
 *
 *  Created by Ashkan Balouchi on 8/18/14.
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
	
	
	
	int num=1; 
	
	int n=30;
	int N=n*1000;  //number of sites
	double d;     //density of cars
	double di=0.11;    //initial density
	double dm=0.13;    //max density 
	double dd=0.0005;    //density step
	double p;     //probability of break
	int vmax = 5;    //max velocity
	int m;     //number of cars
	int totaltime=100000;        //number of totaltime steps
	int initime=500000;     //number of initialization time steps
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
	int jamcar[1000];
	int jamlength[1000];
	
	
	
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
	
	/*filename9="JamLengthAv-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-2.txt";
	std::ofstream resultjav(filename9.c_str());
	
	filename10="NumCarsInJamAv-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-2.txt";
	std::ofstream resultnav(filename10.c_str());*/
	
	
	for (d=di; d<dm; d+=dd) {
		
		
		stringstream ss1;
		ss1 << d;
		
		
		
		
		/*	filename1="VelLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		 std::ofstream resultv1(filename1.c_str());   //initialize file
		 
		 filename2="GapLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		 std::ofstream resultg1(filename2.c_str());   //initialize file
		 
		 filename3="VelLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-2.txt";
		 std::ofstream resultv2(filename3.c_str());   //initialize file
		 
		 filename4="GapLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-2.txt";
		 std::ofstream resultg2(filename4.c_str());   //initialize file
		 
		 filename5="VelLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-3.txt";
		 std::ofstream resultv3(filename5.c_str());   //initialize file
		 
		 filename6="GapLabelTime-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-3.txt";
		 std::ofstream resultg3(filename6.c_str());   //initialize file*/
		
		
	/*	filename7="JamLengthTimeEvo-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-2.txt";
		std::ofstream resultj(filename7.c_str());   //initialize file
		
		filename8="NumCarsInJamTimeEvo-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-2.txt";
		std::ofstream resultn(filename8.c_str());   //initialize file  */
		
		
		filename7="NumberOfJamsHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resulthist(filename7.c_str());   //initialize file
		
		filename8="JamLengthHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resultl(filename8.c_str());   //initialize file

		filename9="JamCarsHist-l"+ssl.str()+"k-v"+sss.str()+"-p0.1-d"+ss1.str()+"-1.txt";
		std::ofstream resultn(filename9.c_str());   //initialize file

		
		
		m=floor(d*N);
		
		initializeuniform(N,d,x,v);
		
		//	Zeroint(jam,N);
		//	Zeroint(jamnum,N);
		
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
			
			for (int i=1; i<jamcount+fjamcount+1; i++) {
				if (jamlength[i] > 0) {
					resultl << jamlength[i] << endl;
					resultn << jamcar[i] << endl;
				//	resulth << head[i] << endl;
				}
			}
			
			
		}			
		   
	
/*		resultjav << d << "	" << (double)jamaux/ (double) totaltime << endl;
		resultnav << d << "	" << (double)numaux/ (double) totaltime << endl;  */
		
		
		/*	for (int t=0; t<bettime; t++) {  
		 
		 Accelerate(N,d,vmax,x,v,g);
		 SlowDown(N,d,vmax,x,v,g);
		 Randomize(N,d,p,v,r);
		 MoveCars(N,d,vmax,x,v,g);
		 
		 }  //initializing simulation
		  
		
		 for (int t=0; t<totaltime; t++) {  //main simulation
		  
		 Accelerate(N,d,vmax,x,v,g);
		 SlowDown(N,d,vmax,x,v,g);
		 Randomize(N,d,p,v,r);
		 MoveCars(N,d,vmax,x,v,g);
		
		 
		 if (t % timestep == 0) {
		 
		 for (int i=0; i<m; i++) {
		 resultv2<< i << "	" << t << "	" << v[i] << "	" << endl;
		 resultg2<< i << "	" << t << "	" << g[i] << "	" << endl;
		 
		 }
		 
		 }
		 
		 }			
		 
		 
		 
		 for (int t=0; t<bettime; t++) {  
		 
		 Accelerate(N,d,vmax,x,v,g);
		 SlowDown(N,d,vmax,x,v,g);
		 Randomize(N,d,p,v,r);
		 MoveCars(N,d,vmax,x,v,g);
		 
		 }  //initializing simulation
		 
		 
		 for (int t=0; t<totaltime; t++) {  //main simulation
		 
		 Accelerate(N,d,vmax,x,v,g);
		 SlowDown(N,d,vmax,x,v,g);
		 Randomize(N,d,p,v,r);
		 MoveCars(N,d,vmax,x,v,g);
		 
		 
		 if (t % timestep == 0) {
		 
		 for (int i=0; i<m; i++) {
		 resultv3<< i << "	" << t << "	" << v[i] << "	" << endl;
		 resultg3<< i << "	" << t << "	" << g[i] << "	" << endl;
		 
		 }
		 
		 }
		  
		 }		*/	
		
		
		resulthist.close();
		resultn.close();
		resultl.close();
	//	resultj.close();
		/*	resultv1.close();
		 resultg1.close();
		 resultv2.close();
		 resultg2.close();
		 resultv3.close();
		 resultg3.close(); */
		
		
		cout << "l = " << n  << "  num = " << num << "  d = " << d << "  simulation done!"<<endl;
		
	}
	
	
	// resultjav.close();
	// resultnav.close();
	
	
	
	
	return (0);
	
	
}










