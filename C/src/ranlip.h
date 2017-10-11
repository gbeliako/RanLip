/************ CRanLip - universal multivariate random  ************************
 *      variate generatior based on acceptance/rejection
 *
 *      begin                : April 30 2004
 *		version				 : 1.0 
 *		copyright            : (C) 2004 by Gleb Beliakov
 *		email                : gleb@deakin.edu.au
 *
 *
 *  Declaration of CRanLip class 
 * 
 *  Description:  CRanLip is a class which implements a universal nonuniform 
 *     multivariate random variate generator based on acceptance/rejection.
 *     Random varites are generated with the help of a hat function - an 
 *     approximation to the given density function from above. RanLip uses
 *     a piecewise constant hat function, calculated by subdividing the
 *     domain into hypercubes, and by using a large number of sample points
 *     and the Lipschitz constant of the probability density. It is assumed
 *     than the density function is Lipschitz continuous.
 *
 *     This method subdivides the rectangular domain into a number of hyper-
 *     rectangles, computes on each element of this partition a large number
 *     of values of the density function, and estimates its absolute maximum
 *     on this hyperrectangle using the Lipschitz constant of the density.
 *     This maximum is taken as the value of the hat function.
 *
 *  Functionality: CRanLip performs two main functions: compute the hat
 *    function (using PrepareHatFunction() method), and then generate a random
 *    variate (using RandomVec() method) with the required density.
 *
 *  Typical usage: The user needs to derive a new class from CRanLip, in order
 *    to provide the code for evaluation of the required density. The user needs
 *    to override Distribution(double* p) method, as shown below.
 *
 *    Then the user calls: Init(dimension, left_boundary, right_boundary)
 *                         PrepareHatFunction(num_elements, fine_partition, Lipschitz)
 *                         RandomVec(random_vector)
 *
 *    The first call specifies the domain and the dimension
 *    The second call prepares the hat function
 *    The third call returns random vectors with the desired density
 *
 *	Example:

	class MyRnumGen:public CRanLip {
		 public: virtual double	Distribution(double* p) ;
	};

	double	MyRnumGen::Distribution(double* p) 
	{ // multivriate normal distribution
		double r;
			for(int j=0;j<Dimension;j++) {
				r+=sqr_(p[j]);
			} 
		return exp(-r);;
	}
	// declare an instance of this class
	MyRnumGen MyGen;
	#define  dim 3
	double left[dim],right[dim]; // boundaries of the domain
	double R[dim];
	for(i=0;i<dim;i++) {left[i]=-2.0; right[i]=3.0;}  // domain is [-2,3]^dim

	MyGen.Init(dim,a,b);
	MyGen.PrepareHatFunctionAuto(10,32); // will be in total 10^dim elements in the partition, and (32*10)^dim
							// evaluations of the density function (calls to Distribution() method)
	for(i=1;i<= TotalRandomVariates;i++)
	   MyGen.RandomVec(R); // generate random variates

 *
 *   In the above example, Lipschitz constant is unknown and is estimated automatically. If it is known
 *   at least approximately (must be an overestimate), then use MyGen.PrepareHatFunction(10,32,LipConst)
 *
 *
 *  Computational complexity: The user has total control of the number of evaluations of the density
 *     function, when s/he supplies arguments num and numfine in PrepareHatFunction(num,numfine, LipConst)
 *     The domain (hyperrectangle) is subdivided into num^dim hyperrectangles, each of which is further
 *     subdivided into numfine^dim rectangles. The fine partition is used to compute the absolute maximum
 *     of the density on the elements of the rough partition. The total number of computations is
 *     (num*numfine)^dim. These values are compared to each other in (num*numfine*dim)^dim comparisons.
 *
 *     The memory required by the algorithm is numfine^dim + 2*num^dim
 *
 *     numfine should be a power of 2 to speed up computations (CRanLip will enforce it). It can be >=2.
 *     However, high values of num translate into less efficient generation, because of very long
 *     tables needed for the discrete random number generator.
 *     
 *     It is important to realise that the quality of the hat function directly depends on the number of     
 *     elements in the partitions, and the user should find a right balance between the quality and the  
 *     time/memory required to build the hat function. Better hat functions require more computing time
 *     but then the speed of random variate generation is improved because of smaller rejection constant.
 *     A more accurate estimate of the Lipschitz constant may also reduce the rejection constant.
 *
 *     For dimensions >3 the preprocessing time may be significant. The tables in the accompanying
 *     paper (ranlip.ps) provide some indications.
 *
 *   Other functionality: The user can store the computed hat function in a file using SavePartition()  
 *     method. LoadPartition() will load the hat function. In this case no preprocessing is required, ie 
 *     the sequence of commands is MyGen.LoadPartition("partitionfile");
 *	                               MyGen.RandomVec(R); // generate random variates
 *
 *     Statistics for generation step can be computed from MyGen.count_total (this is the total number
 *     of calls to the uniform random vector generator (ie dim calls to uniform random number generator))
 *     total number of generated variates / count_total gives the acceptance ratio. Higher ratio indicates
 *     higher efficiency of the generation step.
 *
 *     MyGen.count_errors contains the number of errors due to incorrect specification of the Lipschitz
 *     constant (and hence incorrect hat function). count_errors is incremented when the value of the hat 
 *     function is below the density at this point. Non-zero value indicates that generation was performed
 *     incorrectly, and should be re-done after computing the hat function with a higher Lipschitz
 *     constant. Although some user may choose to ignore low error count.
 *
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

#if !defined(RLIPNUMGEN_H)
#define RLIPNUMGEN_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlipdist.h"

// returns a uniform random number in (0,1)
typedef   double (*UFunction)() ;

#define	sqr_(a) ((a)*(a))
#define	max_(a,b) ((a>b)?(a):(b))

class CRanLip  
{
public:
	CRanLip();
	virtual ~CRanLip();

	UFunction   UGen;		// pointer to the uniform random number generator
	int		Dimension;	// the number of variables
			// Lipschitz constant, volume of each hypercube and an array of probabilities
	        // for the discrete generator
	double  Lipschitz, Volume, *Probabilities;
	int		TotalElements; // the number of hypercubes in the subdivision ofthe domain
	long int count_total, count_errors; // to provide some statistics

	gsl_ran_discrete_t *Dist;  // points to the discrete random variate generator
	size_t		m_chosenElement;         // chosen element of the partition

// these are just global vectors to save time on memory allocation
	double	*V;
	double	*m_boundLeft, *m_boundRight; // left and right bounds of the domain
	double	*m_tempLeft, *m_tempRight; 


// initialises the tables and arrays for the generator
// should be the first method to call
	void Init(int dim, double* left, double* right);

// to free memory occupied by the partition and tables for the alias method
// called from the destructor, but can also be called before
	void FreeMem();

// sets the pointer to the uniform random number generator. The default is
// M. Luescher's ranlux generator (see gsl_randist.h)
	void SetUniformGenerator(UFunction gen) {UGen=gen;};

// computes the hat function given the Lipschitz constant
// num refers to the number of cells the domain is subdivided in each direction
// the total number of cells will be num^dim
// numfine is the number of cells in the fine partition. Should be a power of 2
// to speed up calculations
	void PrepareHatFunction(int num, int numfine, double Lip);
// computes the hat function, and automatically computes the Lipschitz constant, no smaller than minLip
	void PrepareHatFunctionAuto(int num, int numfine, double minLip=0);

// generates a random variate with the required density
	void RandomVec(double* p);

// generates a random vector uniform on (0,1)^dim
	void RandomVecUniform(double* p);

// sets the seed for the uniform random number generator
	void Seed(int seed);
// returns the value of the seed
	int	 GetSeed() {return TheSeed;};

// returns a uniform random number on (0,1) (the default generator is ranlux )
	inline double	UniformRNumber() { return UGen(); };
//  returns a uniform random number on (a,b) 
	inline double	UniformRNumber(double a, double b){ return UGen()*(b-a)+a;};

// saves the computed hat function into a file
	void	SavePartition(char * fname);

// loads the hat function from the file
	void	LoadPartition(char * fname);

// computes the probability density function. Should be provided by the user. The default method returns 1.
	virtual double	Distribution(double* p) ;

private:

// the following methods are to compute the indices in a multidimensional array
	void	GetIJK(int b);
	void	GetIJKfineBin(int b);
	int		GetIndexfromIJK( int* IJK);

// to compute the value of the hat function on an element of the partition
	double	ComputeMaxBin();
// to estimate the Lipschitz constant of the density on an element of the partition
	double	ComputeLipschitzBin();

// compute the values of the density function in the nodes of a finemesh
	void	ComputeArray();
// compute the values of the density function in the nodes of the rough mesh
	void	ComputeArrayCache(int J);


	unsigned int mask1, mask2, bits;
	double	*h, *hfine;
	int		*m_tempint,*m_tempintfine;
	int		*m_delta;
	int		Computed;  // 1 if computed, 0 otherwise
	int		num_partition, num_small_partition,num_small_partition_p1;
	int		TheSeed;
	double	*vals;
	int		totvals;
	double*	LipschitzH;
	double* cache;

};

#endif // !defined(RLIPNUMGEN_H)
