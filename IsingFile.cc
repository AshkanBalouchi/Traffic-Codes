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


double* c=new double[1500000];

int n = 24;  //lattice side size
const int N=n*n*n;  //number of spins
int num;

int i,j;
int E;  //Energy


int main(int argc, char* argv[]){
	
	time_t time0, *time_ptr = 0;       
	time0 = time(time_ptr); 
	
	for (i=0; i<1500000; i++) {
		c[E]=0;
	}
	
	
	
	string filename1,filename2;
	
	
	stringstream ss;
	ss << n;
	
	
	for (k=0.221654; k < 0.223; k+=0.004) {
		
		stringstream ssk;
		ssk << k;
			
		for (num=1; num<5; num++) {

			stringstream runnum;
			runnum << num;
		
			filename1="ising-EHistogram-l"+ss.str()+"-k"+ssk.str()+"-run"+runnum.str()+".txt";
			std::ifstream data(filename1.c_str());   //opnes file
			
			while (!data.eof) {
				data >> E;
				c[-E] ++ ;
			}
			
		}
		
		filename2="ising-EH-l"+ss.str()+"-k"+ssk.str()+".txt";
		std::ofstream result(filename2.c_str());   //opnes file

		for (i=0; i<1500000; i++) {
			if (c[i]>0) {
				cout << -i << "		" << c[i];
			}	
		}	
		
		result.close();
	
	}
	
	
	//----------------------- output simulation time

	double took = difftime(time(time_ptr),time0);
	std::cout << std::endl << "Simulation took " << took << " seconds" << std::endl << std::endl;
	
	
	
	return (0);
}


