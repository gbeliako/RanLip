/************ ranliptest - an example of usage  ************************
 *         of ranlip library using static linking
 *
 *      begin                : May 9 2004
 *		version				 : 1.0 
 *		copyright            : (C) 2004 by Gleb Beliakov
 *		email                : gleb@deakin.edu.au
 *
 *
 *
 *    This example shows the usage of the class interface to access 
 *    ranlip library.  For the class interface, we derive a 
 *    new class from CRanLib (declared in ranlip.h) , called  
 *    MyRnumGen, with the member Distribution(), which implements
 *    the desired distribution density function. The subsequent calls
 *    to the members of CRanLip are as per documentation of this class  
 *
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

#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

#include <math.h>

#include "ranlip.h"


class MyRnumGen:public CRanLip {
public:
	virtual double	Distribution(double* p) ;

};

double	MyRnumGen::Distribution(double* p) 
{ // example: multivariate normal distribution
	double r;
	r=0;
	for(int j=0;j<Dimension;j++) 
		r += p[j]*p[j];
	return exp(-r);
}





int main(int argc, char *argv[])
{
	int dim=3;
	int i,j;
	MyRnumGen MyGen;

	int timesgen=10000;
	double *P=(double*) malloc(dim*sizeof(double)); 

	double *a, *b;
	a=(double*) malloc(dim*sizeof(double)); 
	b=(double*) malloc(dim*sizeof(double));

	for(i=0;i<dim;i++) {a[i]=-1; b[i]=1;}

	MyGen.Init(dim,a,b);
	MyGen.PrepareHatFunction(20,8,2);
	MyGen.Seed(10);
	
// generation step
	for(j=0;j<timesgen;j++) 
		MyGen.RandomVec(P);

	MyGen.SavePartition("partition.txt");

	MyGen.FreeMem();

	MyGen.Init(dim,a,b);
	MyGen.PrepareHatFunctionAuto(20,8);
	MyGen.Seed(10);

// generation step
	for(j=0;j<timesgen;j++) 
		MyGen.RandomVec(P);

	cout <<"Computed Lipschitz const: "<< MyGen.Lipschitz << endl;

	cout<<"Acceptance ratio: "<<  (double)timesgen/(double)MyGen.count_total <<endl;
	cout<<"errors due to low Lipschitz constant: "<<MyGen.count_errors<< endl;

	MyGen.LoadPartition("partition.txt"); 


// generation step
	for(j=0;j<timesgen;j++) 
		MyGen.RandomVec(P);
	cout<< "Loading saved hat function..."<<endl;
	cout<<"Acceptance ratio: "<<  (double)timesgen/(double)MyGen.count_total <<endl;
	cout<<"errors due to low Lipschitz constant: "<<MyGen.count_errors<< endl;

	return 0;
}
