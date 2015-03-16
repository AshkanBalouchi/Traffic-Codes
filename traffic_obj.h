/*
 *  initialize.h
 *  
 *
 *  Created by Ashkan Balouchi on 11/26/12.
 *  Copyright 2012 LSU. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <math.h>
#include <ctime>
#include <cstdlib>





#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>  /* time */

#include "gsl_rng.h"

#include <map>
//#include <random>

using namespace std;

void initialize(int N, double d,int *x,int *v){
	
	int m=floor(d*N);   //# of cars
	
	for (int i=0;i<m; i++) {
		x[i]=i;   //initialize position of cars
		v[i]=0;   //initialize velocity of cars
	}
	
	
}



void initializeuniform(int N, double d,int *x,int *v){
	
	int m=floor(d*N);   //# of cars
	double dis=1/d;
	for (int i=0;i<m; i++) {
		x[i]=floor(i*dis);   //initialize position of cars
		v[i]=0;   //initialize velocity of cars
	}
	
	
}


void initializeuniformOpenBoundary(int N, double d,int vmax,int *x,int *v){
	
	int m=floor(d*N);   //# of cars
	double dis=1/d;
	for (int i=1;i<m; i++) {
		x[i]=floor(i*dis);   //initialize position of cars
		v[i]=0;   //initialize velocity of cars
	}
	
	
}


void initializeOpenBoundary(int N, double d,int vmax,int *x,int *v){
	
		x[0]=0;   //initialize position of cars
		v[0]=vmax;   //initialize velocity of cars
	
	
	
}


void initializeRandom(int N, double d,int *x,int *v){
	
	int m=floor(d*N);   //# of cars
	double random;
	int rndposition;
	int count;
	count =0;
	
	while (count < m) {
		
		random =(double)rand()/(double)RAND_MAX;
		rndposition=floor(N*random);
		if (x[rndposition]==0) {
		count++;
		x[rndposition]=rndposition;   //initialize position of cars
		v[rndposition]=0;             //initialize velocity of cars
	
		}
	}
	
	
}

void initializeRandom2(int N, double d,int *x,int *v){
	
	int m=floor(d*N);   //# of cars
	double random;
	int count;
	count=0;
	for (int i=0; i<N; i++) {
		random =(double)rand()/(double)RAND_MAX;
if (random<d) {
	x[count]=i;
	v[count]=0;
	count++;
}
	}	
	
}



void FindGaps(int N, double d,int vmax,int *x, int *v, int *g){
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m-1; i++) {
		g[i]=x[i+1]-x[i]-1;       //finds gaps 
	}
	g[m-1]=x[0]-x[m-1]-1+N;
	
	for (int i=0; i<m; i++) {
		if (g[i]<0) {
			g[i]+=N;
		}
		if (g[i]>N-1) {
			g[i]=g[i]-N;
			
	}
	}
	
}

void FindGapsOpenBoundary(int N, double d,int m,int vmax,int *x, int *v, int *g){
	
	for (int i=0;i<m-1; i++) {
		g[i]=x[i+1]-x[i]-1;       //finds gaps 
	}
	g[m-1]=N;
	
}



void Accelerate(int N, double d,int vmax,int *x, int *v, int *g){

	FindGaps(N,d,vmax,x,v,g);
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		if (v[i] < vmax ) {  // && v[i]<g[i] ??
			v[i]++;
		}
		
	}	
	}


void AccelerateOpenBoundary(int N, double d,int m,int vmax,int *x, int *v, int *g){
	
	FindGapsOpenBoundary(N,d,m,vmax,x,v,g);
		
	for (int i=0;i<m; i++) {
		if (v[i] < vmax ) {  // && v[i]<g[i] ??
			v[i]++;
		}
		
	}	
}

void Acceleratemax(int N, double d,int vmax,int *x, int *v, int *g){
	
	FindGaps(N,d,vmax,x,v,g);
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		v[i] = vmax;
	}	
}


void AcceleratemaxOpenBoundary(int N, double d,int m,int vmax,int *x, int *v, int *g){
	
	FindGapsOpenBoundary(N,d,m,vmax,x,v,g);
		
	for (int i=0;i<m; i++) {
		v[i] = vmax;
	}	
}

void SlowDown(int N, double d, int vmax, int *x, int *v, int *g){
	
	FindGaps(N,d,vmax,x,v,g);

	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		if (g[i] < v[i]) {
			v[i] = g[i]; //(g[i]-1 ??
		}
		
	}	
}


void SlowDownBreak(int N, double d, int vmax, int &countb, int &countf, int *x, int *v, int *g, int *xb, int *xf){
	
	FindGaps(N,d,vmax,x,v,g);
	countb=0;
	countf=0;
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		if (g[i] < v[i]) {
			v[i] = g[i];
			xb[countb]=x[i];
			countb++;
		} else {
			xf[countf]=x[i];
			countf++;
		}

		
	}	
}


void SlowDownBrakeLight(int N, double d, int vmax, int *x, int *v, int *g, int *bls){
	
	FindGaps(N,d,vmax,x,v,g);

	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		if (g[i] < v[i]) {
			v[i] = g[i];
			bls[i]=1;
		} else {
			bls[i]=0;
		}
		
		
	}	
}

void SlowDownBreakDistribuion(int N, double d, int vmax, int &countb, int *x, int *v, int *g, int *xb, int *nb){
	
	FindGaps(N,d,vmax,x,v,g);
	countb=0;
	
	int m=floor(d*N);     //# of cars
	
	int nbaux;
	
	
	
	for (int i=0;i<m; i++) {
		
		nbaux=nb[i];
		
		if (g[i] < v[i]) {
			v[i] = g[i];
			xb[countb]=x[i];
			countb++;
			nbaux++;
		}
		
		if (nbaux>nb[i]) {
			nb[i]++;
		} else {
			nb[i]=0;
		}

		
		}
		
		
	}	


void SlowDownBreakJam(int N, double d, int vmax, int &countb, int &countf, int *x, int *v, int *g, int *xb, int *xf, int *car){
	
	FindGaps(N,d,vmax,x,v,g);
	countb=0;
	countf=0;
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		if (g[i] < v[i]) {
			v[i] = g[i];
			if (g[i]<vmax/2.) {
				car[countb]=i;
				xb[countb]=x[i];
				countb++;
			}else {
				xf[countf]=x[i];
				countf++;
			}
			
		} else {
			xf[countf]=x[i];
			countf++;
		}
		
		
	}	
}



void SlowDownOpenBoundary(int N, double d,int m,int vmax, int *x, int *v, int *g){
	
	FindGapsOpenBoundary(N,d,m,vmax,x,v,g);

	
	for (int i=0;i<m; i++) {
		if (g[i] < v[i]) {
			v[i] = g[i]; //(g[i]-1 ??
		}
		
	}	
}



void Randomize(int N, double d, double p, int *v, gsl_rng *r){
	
	int m=floor(d*N);     //# of cars
	

	// srand (1);

	
	//std::random_device rd;
	//int rd=1;
	
    // Choose a random mean between 1 and 6
    //default_random_engine e1(rd());
    //uniform_int_distribution<int> uniform_dist(1, 100000);
	
//	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);			

	
	for (int i=0;i<m; i++) {
		if (v[i]>0){
//			double random =(double)rand()/(double)RAND_MAX;
	
			
			//	double random =(rand()% 32000)/32000.0;
			
			
			double random = gsl_rng_uniform (r);
			
		//	cout << random << endl;
			
		//	double random = gsl_rng_uniform ;


		//	double random = uniform_dist(e1)/100000.0;

			
			//cout<<random<<endl;
			
			if(p>random){
				v[i]=v[i]-1; 
		}
		
		}
		
		
}
}

void RandomizeOpenBoundary(int N, double d,int m,double p, int *v){
	
	
	for (int i=0;i<m; i++) {
		if (v[i]>0){
			double random =(double)rand()/(double)RAND_MAX;
			//cout<<random<<endl;
			
			if(p>random){
				v[i]=v[i]-1; 
			}
			
		}
		
		
	}
}

void MoveCars(int N, double d, int vmax,int *x, int *v, int *g){
	
	int m=floor(d*N);     //# of cars
	
	for (int i=0;i<m; i++) {
		x[i]=x[i]+v[i];
		if (x[i]>=N) {
			x[i]=x[i]-N;
		}
		
	}	
}




void MoveCarsOpenBoundary(int N, double d,int &m,int vmax,int *x, int *v, int *g){
	
	int xx[m];
	int vv[m];
	
	double random =(double)rand()/(double)RAND_MAX;
	//cout << random<<endl;
    if ((random<d) && (x[0]+v[0]>0)) {
	
		for (int i=0;i<m; i++) {
			xx[i+1]=x[i]+v[i];
			vv[i+1]=v[i];
		}	
		xx[0]=0;
		vv[0]=std::min(xx[1]-1,vmax);
		m++;
	
	} else {
		for (int i=0;i<m; i++) {
			xx[i]=x[i]+v[i];
			vv[i]=v[i];
		}
	}

	if (x[m-1]>N-1) {
		m=m-1;
	}
	for (int i=0; i<m; i++) {
		x[i]=xx[i];
		v[i]=vv[i];
	}
	
	
}








double average(int *a,int m){
	double A=0;
	double B=0;
	for (int i=0; i<m; i++) {
		A+=a[i];
	}
	B=A/m;
	return (B);
	}



double average2(int *a,int m){     //finds average of the quantity squared
	double A=0;
	for (int i=0; i<m; i++) {
		A+=a[i]*a[i];
	}
	A=A/m;
	return (A);
}


void track(double *track,int *a,int m){  //keeps track of ime evolution to find time average
for (int i=0; i<m; i++) {
	track[i]+=a[i];
}
}



void FindMobility(double *vtrack,double *c,int m, double t){     //  c:mobility
	for (int i=0; i<m; i++) {
		c[i]=vtrack[i]/(t+1);
	}
}


void Zero(double *a,int m){
for (int i=0; i<m; i++) {
	a[i]=0;
}
}

void Zeroint(int *a,int m){
	for (int i=0; i<m; i++) {
		a[i]=0;
	}
}

void FindDynamicCorrelation (double *vtrack,double *c,int m,double t,double *G4){
	FindMobility(vtrack,c,m,t);
	double aux1[m],
	aux2=0;
	for (int k=0; k<m; k++) {
		aux2+=c[k];
	}
	//aux2/=m;
	//aux2=average(c,m);
	Zero(aux1,m);
	Zero(G4,m);

	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			
		aux1[i]+=c[j]*c[i];
	
		}
		G4[i]=(aux1[i]-c[i]*c[i])/(m-1)-(aux2-c[i])*(aux2-c[i])/((m-1)*(m-1));
	}
	
} 


double susceptibility(double *vtrack,int *v,double *c,int m,double t,double *G4){  
	FindDynamicCorrelation(vtrack,c,m,t,G4);
	double X=0;
	for (int i=0; i<m; i++) {
		X+=G4[i];
	}
	double vbar=average(v,m);
	double v2bar=average2(v,m);
	double X4=(1/(v2bar-vbar*vbar))*X;
	return (X4);
	
}


int power(int a, int b){
	int aux=1;
	for (int i=0; i<b; i++) {
		aux*=a;
	}
	return(aux);
}

double round(double a,int b){
	double aux1,aux2;
	aux2=a*power(10,b);
	if (aux2-(floor(aux2)) <= 0.5) {
		aux1=floor(aux2)/power(10,b);
	} else {
		aux1=(floor(aux2)+1)/power(10,b);
	}
	return (aux1);


}


void FindVelDistribution(int N, double d,int vmax,int *v,double *vd){
	
	int m=floor(d*N);     //# of cars
	
	Zero(vd,N);
	
	for (int i=0;i<m; i++) {
		vd[v[i]]++;       //finds gap distributions 
	}
}

void FindGapDistribution(int N, double d,int *g,double *gapd){
	
	int m=floor(d*N);     //# of cars
	
	Zero(gapd,N);
	
	for (int i=0;i<m; i++) {
		gapd[g[i]+1]++;       //finds gap distributions 
	}
}

void FindVelDistributionOpenBoundary(int N, double d,int m,int vmax,int *v,double *vd){
		
	Zero(vd,N);
	
	for (int i=0;i<m; i++) {
		vd[v[i]]++;       //finds gap distributions 
	}
}

void FindGapDistributionOpenBoundary(int N, double d,int m,int *g,double *gapd){
		
	Zero(gapd,N);
	
	for (int i=0;i<m; i++) {
		gapd[g[i]+1]++;       //finds gap distributions 
	}
}


void Find2GapDistribution(int N, double d,int *g,double *gapd){
	
	int m=floor(d*N);     //# of cars
	
	Zero(gapd,N);
	
	for (int i=0;i<m-1; i++) {
		gapd[g[i]+g[i+1]+2]++;       //finds gap distributions 
	}
	gapd[g[m-1]+g[0]+2]++;
}

void Find3GapDistribution(int N, double d,int *g,double *gapd){
	
	int m=floor(d*N);     //# of cars
	
	Zero(gapd,N);
	
	for (int i=0;i<m-2; i++) {
		gapd[g[i]+g[i+1]+g[i+2]+3]++;       //finds gap distributions 
	}
	gapd[g[m-2]+g[m-1]+g[0]+2]++;
	gapd[g[m-1]+g[0]+g[1]+2]++;

}

	
void CorrectGapDistribution(int N,int m,int i,long double *pg,long double *pgc){
	

	
	for (int j=1; j<N; j++) {
	if (j==i) {
		if (m*pg[i]-1>0) {
			pgc[j]=(m*pg[i]-1)/(double(m)-1);}
		else {
			pgc[j]=0;
		}

	}else {
		pgc[j]=m*pg[j]/(m-1);}
		
	}
	}
	
	




void ClusterCount(int N, double d, int vmax, int &jamcount, int &fjamcount, int *jamcar,int *jamlength, int *x, int *v, int *g){
	
	jamcount=0;	
	int m=floor(d*N);     //# of cars
	int jid[m];
	/*
	int head[m],tail[m],dis[m];
	
	Zeroint(jamlength , m);
	Zeroint(jamcar , m);
	Zeroint(dis  , m);
	Zeroint(head , m);
	Zeroint(tail , m);*/
	
	int hh = 5000;
	
	int head[hh],tail[hh],dis[hh];
	
	Zeroint(jamlength , hh);
	Zeroint(jamcar , hh);
	Zeroint(dis  , hh);
	Zeroint(head , hh);
	Zeroint(tail , hh);
	
	
	for (int i=0;i<m; i++) {
		if (g[i] < vmax/2.+1) {
			jid[i] = 1;
		} else {
			jid[i]=0;
		}
	}	 // this loops gives jam ID to the cars
	
	
/*	int j=0;
	while (jid[j] == 0) {
		j++;
	}*/
	
	
	if (jid[0] == 0) {
		
		for (int i=0;i<m-1; i++) {
			
			if (jid[i+1] > jid[i]) {
				jamcount ++ ;
				tail[jamcount] = i;
			}else {
				if (jid[i+1] < jid[i]) {
					head[jamcount] = i;
				}
			}
		
		}
		
		
		if (jid[0] > jid[m-1]) {
			jamcount ++ ;
			tail[jamcount] = m-1;
		}else {
			if (jid[0] < jid[m-1]) {
				head[jamcount] = m-1;
			}
		}
		
		
	} else {
		for (int i=0;i<m-1; i++) {
			
			if (jid[i+1] < jid[i]) {
				jamcount ++ ;
				head[jamcount] = i;
			}else {
				if (jid[i+1] > jid[i]) {
					tail[jamcount+1] = i;
				}
			}
				
		}
		
		if (jid[0] < jid[m-1]) {
			jamcount ++ ;
			head[jamcount] = m-1;
		}else {
			if (jid[0] > jid[m-1]) {
				tail[jamcount+1] = m-1;
			}
		}
		
		
		
		tail[1] = tail[jamcount+1];
		tail[jamcount+1]=0;
	}
	
	
	for (int i=1; i < jamcount; i++) {
		jamlength[i] = x[head[i]] - x[tail[i]]+1;
		dis[i] = x[tail[i+1]]-x[head[i]];
		jamcar[i]=head[i] - tail[i]+1;
	}
	
	jamlength[jamcount] = x[head[jamcount]] - x[tail[jamcount]]+1;
	dis[jamcount] = x[tail[1]]-x[head[jamcount]];
	jamcar[jamcount]=head[jamcount]-tail[jamcount]+1;

	for (int i=1; i < jamcount+1; i++) {
		
		if (jamlength[i] < 0) {
			jamlength[i] += N;
		}
		if (dis[i] < 0) {
			dis[i] += N;
		}
		if (jamcar[i] < 0) {
			jamcar[i] += m;
		}
		
	}
	

	
	int k=0;
	int aux[m];

	if (jamcount > 0) {
		
	
		for (int i=1; i < jamcount; i++) {
			if (dis[i] < jamlength[i]) {
				jamlength[i+1] = jamlength[i]+dis[i]+jamlength[i+1];
				jamlength[i]=0;
				dis[i]=0;
				
				jamcar[i+1] = head[i+1]-tail[i]+1;  
				if (jamcar[i+1] < 0) {
					jamcar[i+1] += m;
				}
				jamcar[i] = 0;
				aux[k]=i;
				k++;
			}
		}
		
		if (dis[jamcount] < jamlength[jamcount]) {
			jamlength[1] = jamlength[jamcount]+dis[jamcount]+jamlength[1];
			jamlength[jamcount]=0;
			dis[jamcount]=0;
			
			jamcar[1] = head[1]-tail[jamcount]+1;
			if (jamcar[1] < 0) {
				jamcar[1] += m;
			}
			jamcar[jamcount] = 0;
			aux[k]=jamcount;
			k++;
		}
	
		
		for (int i=1; i<jamcount+1; i++) {
			if (jamcar[i]>0 && jamcar[i] < 3) {
				k++;
				
				jamlength[i] = 0;
				jamcar[i]=0;
				
				
			}
		}
			
	}
	

	jamcount = jamcount - k;
	
	fjamcount = k;
	
	
}

