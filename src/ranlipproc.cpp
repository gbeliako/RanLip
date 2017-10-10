/************ RanLip - universal multivariate random  ************************
 *      variate generatior based on acceptance/rejection
 *           Procedural interface implementation
 *
 *      begin                : May 10 2004
 *		version				 : 1.0 
 *		copyright            : (C) 2004 by Gleb Beliakov
 *		email                : gleb@deakin.edu.au
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



#include "ranlip.h"
#include "ranlipproc.h"


class ProcRanLip:public CRanLip {
	 public: virtual double	Distribution(double* p) ;
			 DistFun	TheDistribution;
			 ProcRanLip() {TheDistribution=NULL;}

			 void SetDistFunction(DistFun dist);
};


double	ProcRanLip::Distribution(double* p) 
{ // pass execution to the supplied function
	if(TheDistribution!=NULL)
	return TheDistribution(p,Dimension);
	else return 1;
}

void ProcRanLip::SetDistFunction(DistFun dist) {TheDistribution=dist;}


ProcRanLip RanLipGenerator;


// initialises the tables and arrays for the generator
// should be the first method to call
void InitRanLip(int dim, double* left, double* right)
{
	RanLipGenerator.Init(dim,left,right);
}


void SetDistFunctionRanLip(DistFun dist)
{
	RanLipGenerator.SetDistFunction(dist);
}

// computes the hat function given the Lipschitz constant
// num refers to the number of cells the domain is subdivided in each direction
// the total number of cells will be num^dim
// numfine is the number of cells in the fine partition. Should be a power of 2
// to speed up calculations
void PrepareHatFunctionRanLip(int num, int numfine, double Lip)
{
	RanLipGenerator.PrepareHatFunction(num,numfine,Lip);
}

// computes the hat function, and automatically computes the Lipschitz constant
void PrepareHatFunctionAutoRanLip(int num, int numfine, double minLip)
{
	RanLipGenerator.PrepareHatFunctionAuto(num,numfine, minLip);
}

// generates a random variate with the required density
void RandomVecRanLip(double* p)
{
	RanLipGenerator.RandomVec(p);
}

// saves the computed hat function into a file
int	SavePartitionRanlip(char * fname)
{
	return RanLipGenerator.SavePartition(fname);
}

// loads the hat function from the file
int	LoadPartitionRanLip(char * fname)
{
	return RanLipGenerator.LoadPartition(fname);
}

// to free memory occupied by the generator
void FreeMemRanLip()
{
	RanLipGenerator.FreeMem();
}

// sets the pointer to the uniform random number generator. The default is
// M. Luescher's ranlux generator (see gsl_randist.h)
void SetUniformGeneratorRanLip( UFunction gen )
{
	RanLipGenerator.SetUniformGenerator(gen);
}

// sets the seed for the uniform random number generator
void SeedRanLip(int seed)
{
	RanLipGenerator.Seed(seed);
}

// gets the seed used for the uniform random number generator
int GetSeedRanLip()
{
	return RanLipGenerator.GetSeed();
}

// returns the total number of generated random variates (accepted and rejected)
long int Count_totalRanLip()
{
	return  RanLipGenerator.count_total;
}

// returns the number of points at which the hat function is smaller than the density
// if nonzero, the random sequence is incorrect (Lipschitz constant was too low).
long int Count_errorRanLip()
{
	return  RanLipGenerator.count_errors;
}

// returns Lipschitz constant computed by PrepareHatFunctionAutoRanLip
double LipschitzRanLip()
{
	return RanLipGenerator.Lipschitz;
}

int	GetDim()
{
	return RanLipGenerator.Dimension;
}