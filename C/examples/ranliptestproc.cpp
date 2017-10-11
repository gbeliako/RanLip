/************ ranliptestproc - an example of usage  ************************
 *           of ranlip library using static linking
 *
 *      begin                : May 9 2004
 *		version				 : 1.0 
 *		copyright            : (C) 2004 by Gleb Beliakov
 *		email                : gleb@deakin.edu.au
 *
 *
 *    This example shows the usage of procedural interface to access
 *    ranlip library. 
 *
 *    For the procedural interface, we implement caculation of the density  
 *    in double MyDist(double* p, int dim) function and pass its address
 *	  in  SetDistFunctionRanLip(&MyDist) procedure. We also declare our own
 *    uniform random number generator MyRand() ad pass its address in
 *    SetUniformGeneratorRanLip(&MyRand) procedure.
 *
 *    Then we call the required functions as described in  ranlip.h
 *
 *
 * Copyright (C) 2004 Gleb Beliakov (gleb@deakin.edu.au)
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;

#include "ranlipproc.h"

// to measure preprocessing and generation time
clock_t clockS,clockF;
double TotalTime;
void   ResetTime() {TotalTime=0; clockS=clock();}
double ElapsedTime()
{
	clockF=clock();
	double duration=(double)(clockF - clockS) / CLOCKS_PER_SEC;
	TotalTime += duration;
	clockS=clockF;
	return TotalTime;
}

// we need to implement calculation of the probability density in this function

double MyDist(double* p, int dim)
{ // example: multivariate normal distribution
	double r;
	r=0;
	for(int j=0;j<dim;j++) 
		r += p[j]*p[j];
	return exp(-r);
}

double MyRand()
{  // define my own random number generator to replace ranlux
	return  (double)rand()/(RAND_MAX+0.001); 
}

int main(int argc, char *argv[])
{
	int dim=2;
	int i,j;

// box constraints live here ------------------
	double *a, *b;
	a=(double*) malloc(dim*sizeof(double)); 
	b=(double*) malloc(dim*sizeof(double));
	for(i=0;i<dim;i++) {a[i]=0; b[i]=2;}

	ResetTime();

	InitRanLip(dim,a,b);
	SetDistFunctionRanLip(&MyDist);

	SetUniformGeneratorRanLip(&MyRand);
	srand(10); // not SeedRanLip(10);

//	PrepareHatFunctionRanLip(20,4,1);
	PrepareHatFunctionAutoRanLip(50,16);
	cout << "preprocesing time "<<ElapsedTime()<<endl;

	cout<<"Lipschitz constant is "<<LipschitzRanLip()<<endl;
	SeedRanLip(10);

	double rr;
	int timesgen=10000;

	double *P;
	P=(double*) malloc(dim*sizeof(double)); 

	ResetTime();
	for(j=0;j<timesgen;j++)
		RandomVecRanLip(P);

	rr=ElapsedTime();
	cout<<" average generation time: "<< rr/timesgen<<endl;

	cout<<"Acceptance ratio: "<<(double)timesgen/(double)Count_totalRanLip()<<endl;
	if(Count_errorRanLip()>0) 
		cout<<"Lipschitz constant too small, number of errors:"<<Count_errorRanLip()<<endl;

	return 0;
}

