/************ RanLip - universal multivariate random  ************************
 *      variate generatior based on acceptance/rejection
 *                     Procedural interface
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



typedef   double (*UFunction)() ;

// computes the probability density function. Should be implemented by the user.
// p is a vector of size dim
typedef   double (*DistFun)(double* p, int dim) ;


#ifdef __cplusplus
extern "C" {
#endif

// initialises the tables and arrays for the generator
// should be the first method to call
	void InitRanLip(int dim, double* left, double* right);

// sets the pointer to the function which computed the density
	void SetDistFunctionRanLip(DistFun dist);

// computes the hat function given the Lipschitz constant
// num refers to the number of cells the domain is subdivided in each direction
// the total number of cells will be num^dim
// numfine is the number of cells in the fine partition. Should be a power of 2
// to speed up calculations
	void PrepareHatFunctionRanLip(int num, int numfine, double Lip);

// computes the hat function, and automatically computes the Lipschitz constant
	void PrepareHatFunctionAutoRanLip(int num, int numfine, double minLip);

// generates a random variate with the required density
	void RandomVecRanLip(double* p);

// saves the computed hat function into a file
	int	SavePartitionRanlip(char * fname);

// loads the hat function from the file
	int	LoadPartitionRanLip(char * fname);

// to free memory occupied by the generator
	void FreeMemRanLip();

// sets the pointer to the uniform random number generator. The default is
// M. Luescher's ranlux generator (see gsl_randist.h)
	void SetUniformGeneratorRanLip( UFunction gen  );

// sets the seed for the uniform random number generator
    void SeedRanLip(int seed);

// gets the seed used for the uniform random number generator
    int GetSeedRanLip();

// returns Lipschitz constant computed by PrepareHatFunctionAutoRanLip
    double LipschitzRanLip();
 
// returns the total number of generated random variates (accepted and rejected)
    long int Count_totalRanLip();

// returns the number of points at which the hat function is smaller than the density
// if nonzero, the random sequence is incorrect (Lipschitz constant was too low).
    long int Count_errorRanLip();

	int	GetDim();

#ifdef __cplusplus
}
#endif
