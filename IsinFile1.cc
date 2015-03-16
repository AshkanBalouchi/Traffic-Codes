/*
 *  IsingFile.cc
 *  
 *
 *  Created by Ashkan Balouchi on 10/6/14.
 *  Copyright 2014 LSU. All rights reserved.
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
#include <time.h>		/* time */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <map>



using namespace std;

int n = 16;  //lattice side size
int num;
int maxnum=2;

int maxe=1500000;

double* m=new double[maxe];
double* absm=new double[maxe];
double* m2=new double[maxe];
double* m4=new double[maxe];
double* c=new double[maxe];



const int N=n*n*n;  //number of spins


int i,j;
int E;  //Energy

double k;


int main(int argc, char* argv[]){
	
	time_t time0, *time_ptr = 0;       
	time0 = time(time_ptr); 
	
	
	for (i=0; i<maxe; i++) {
		m[i]=0;
		absm[i]=0;
		m2[i]=0;
		m4[i]=0;
		c[i]=0;
	}
	
	
	string filename1,filename2;
	
	
	stringstream ss;
	ss << n;
	
	
	for (k=0.221654; k < 0.223; k+=0.004) {
		
		stringstream ssk;
		ssk << k;
		
		for (num=1; num < maxnum+1 ; num++) {
			
			stringstream runnum;
			runnum << num;
			
			filename1="ising-EMHL-l"+ss.str()+"-k"+ssk.str()+"-run"+runnum.str()+".txt";
			std::ifstream data(filename1.c_str());   //opnes file
			
			while (data >> E >> count >> mag >> absmag >> mag2 >> mag4) {
				if (c[-E] = 0) {
					c[-E] = count;
					m[-E] = mag;
					abs[-E] = absmag;
					m2[-E] = mag2;
					m4[-E] = mag4;
				}else {
					m4[-E] = (m4[-E]*c[-E]+mag4*count)/(c[-E]+count);
					m2[-E] = (m2[-E]*c[-E]+mag2*count)/(c[-E]+count);
					absm[-E] = (absm[-E]*c[-E]+absmag*count)/(c[-E]+count);
					m[-E] = (m[-E]*c[-E]+mag*count)/(c[-E]+count);
					c[-E] += count;
					
					}
			}
			
		}
		
		filename2="ising-EMHL-l"+ss.str()+"-k"+ssk.str()+".txt";
		std::ofstream result(filename2.c_str());   //opnes file
		
		for (i=0; i<1500000; i++) {
			if (c[i]>0) {
				result << -i << "		" << c[i]<< "	" << i << "		" << m[i] << "		" << absm[i] << "		" << m2[i] << "		" << m4[i] << endl;
			}	
		}	
		
		result.close();
		
	}
	
	
	//----------------------- output simulation time
	
	double took = difftime(time(time_ptr),time0);
	std::cout << std::endl << "Simulation took " << took << " seconds" << std::endl << std::endl;
	
	
	
	return (0);
}


/*
 *  IsinFile1.cc
 *  
 *
 *  Created by Ashkan Balouchi on 10/12/14.
 *  Copyright 2014 LSU. All rights reserved.
 *
 */


