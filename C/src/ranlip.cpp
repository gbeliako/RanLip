/************ CRanLip - universal multivariate random  ************************
 *      variate generatior based on acceptance/rejection
 *
 *      begin                : April 30 2004
 *		version				 : 1.0 
 *		copyright            : (C) 2004 by Gleb Beliakov
 *		email                : gleb@deakin.edu.au
 *
 *
 *  implementation of CRanLip class 
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


extern "C" int TheSeed;
extern 	ranlux_state_t RLSTATE;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CRanLip::CRanLip()
{
	Computed=0;
	SetUniformGenerator(ranlux_get_double_V);

	Probabilities=NULL;
	Dist=NULL;
	m_boundLeft=NULL;
	m_boundRight=NULL;
	m_tempLeft=NULL;
	m_tempRight=NULL;
	m_tempint=NULL;
	h=hfine=NULL;
	Lipschitz=0;
	TheSeed=10;
}

CRanLip::~CRanLip() {FreeMem();}

void CRanLip::FreeMem()
{
	if(Dist !=NULL)			gsl_ran_discrete_free(Dist);
	if(Probabilities!=NULL)	free(Probabilities);
	if(m_boundLeft!=NULL)	free(m_boundLeft);
	if(m_boundRight!=NULL)	free(m_boundRight);
	if(m_tempLeft!=NULL)	free(m_tempLeft);
	if(m_tempRight!=NULL)	free(m_tempRight);
	if(m_tempint!=NULL)		free(m_tempint);
	if(m_tempintfine!=NULL)	free(m_tempintfine);
	if(h!=NULL)				free(h);
	if(hfine!=NULL)			free(hfine);
	if(V!=NULL)				free(V);

	Computed=0;
	Probabilities=NULL;
	Dist=NULL;
	m_boundLeft=NULL;
	m_boundRight=NULL;
	m_tempLeft=NULL;
	m_tempRight=NULL;
	m_tempint=m_tempintfine=NULL;
	h=hfine=V=NULL;
}



double	CRanLip::Distribution(double* p) // dummy method
{	return 1; }

void CRanLip::Seed(int seed)
{
	TheSeed=seed;
	ranlux_set_seed(seed);
	count_total=	count_errors=0;
}

void CRanLip::Init(int dim, double* left, double* right)
{
	m_boundLeft=(double*) malloc(dim*sizeof(double)); 
	m_boundRight=(double*) malloc(dim*sizeof(double));
	m_tempLeft=(double*) malloc(dim*sizeof(double));
	m_tempRight=(double*) malloc(dim*sizeof(double));
	V =(double*) malloc(dim*sizeof(double));
	h =(double*) malloc(dim*sizeof(double));
	hfine =(double*) malloc(dim*sizeof(double));

	m_tempint=(int*) malloc(dim*sizeof(int));
	m_tempintfine=(int*) malloc(dim*sizeof(int));

	Dimension=dim;
	int i;
	for(i=0;i<Dimension;i++) 
	{
		m_boundLeft[i]=left[i];
		m_boundRight[i]=right[i];
	}
	Computed=0;
}

void CRanLip::RandomVecUniform(double* p)
{
	size_t j = m_chosenElement = gsl_ran_discrete(Dist);

	int i;
	for(i =0;i<Dimension;i++) V[i] = UniformRNumber();

	// calculate the position i,j,k,... in the array, starting at 0
	GetIJK(j);

//	m_tempint contains the indices
	for(i =0;i<Dimension;i++) p[i] = m_boundLeft[i] + h[i]*m_tempint[i];
	for(i =0;i<Dimension;i++) p[i] += V[i]*h[i];

	count_total++;
}

void CRanLip::RandomVec(double* p)
{
	while(Computed) {  // must be 1 to generate random variates
		RandomVecUniform(p);
		double Y = UniformRNumber();
		Y *= Probabilities[m_chosenElement];

		if(Distribution(p) > Probabilities[m_chosenElement]) //return;
		{  count_errors++; 
//		    printf("fail %f %f\n",Distribution(p),Probabilities[m_chosenElement]); 
		return;}

		if(Y <= Distribution(p)) {/*count_success++;*/ return;}
	}
}


void CRanLip::PrepareHatFunction(int num, int numfine, double Lip)
{
	int i,j;
	double d;

	num=max_(num,1);
	numfine=max_(numfine-1,1);

	for(bits=1;bits<32;bits++) if((1<<bits) >= numfine+1) break;
	numfine=(1<<bits) -1;
	mask1= (1<<(bits)) - 1;

	Lipschitz=max_(Lip,1e-10);


	TotalElements=1;
	num_partition=num;   num_small_partition=numfine;  num_small_partition_p1=num_small_partition+1;
	for(i=0;i<Dimension;i++) TotalElements*=num; 

	totvals=1;
	for(i=0;i<Dimension;i++) totvals*=num_small_partition_p1; 

	Probabilities=(double*) malloc(TotalElements*sizeof(double));

// working arrays
	LipschitzH=(double*)malloc(Dimension*sizeof(double));
	m_delta=(int*)malloc(Dimension*sizeof(int));
    vals=(double*)malloc(totvals*sizeof(double));


	Volume=1;
	m_delta[Dimension-1]=1;

	for(i=0;i<Dimension;i++)
	{
		h[i]=(m_boundRight[i]-m_boundLeft[i])/num;
		hfine[i]=h[i]/numfine;
		Volume *= h[i];

		LipschitzH[i]=0.5*Lipschitz*hfine[i] *2; //will divide by 2 later

		if(i>0) m_delta[Dimension-i-1] = m_delta[Dimension-i] * num_small_partition_p1;
	}

// make a special case if numfine==2, not to re-compute many values, store them in a table
// modify ComputeArray by using cache
	if(numfine<=1) {
		cache=(double*)malloc(TotalElements*sizeof(double));
		for(j=0;j<TotalElements;j++) {
			GetIJK(j);
			for(i=0;i<Dimension;i++) 
				m_tempLeft[i] = m_boundLeft[i] + h[i]*m_tempint[i];  // left top corner
			cache[j]=Distribution(m_tempLeft);
		}

		for(j=0;j<TotalElements;j++) {
			ComputeArrayCache(j);
			d=ComputeMaxBin();
			Probabilities[j]=d;
		}
		free(cache);
	} // special case
	else { // general case
		for(j=0;j<TotalElements;j++) {
			GetIJK(j);
			for(i=0;i<Dimension;i++) {
				m_tempLeft[i] = m_boundLeft[i] + h[i]*m_tempint[i];  // left top corner
	//			m_tempRight[i]=m_tempLeft[i]+h[i];  // bottom right corner
			}
			ComputeArray();
			d=ComputeMaxBin();
			Probabilities[j]=d;
		}
	}

	free(vals);
	free(LipschitzH);
	free(m_delta);

	for(j=0;j<TotalElements;j++) Probabilities[j] *= Volume;

	Dist = gsl_ran_discrete_preproc(TotalElements,Probabilities);

	// now get back to function values in probabilities
	for(j=0;j<TotalElements;j++) Probabilities[j] /= Volume;

	count_total=	count_errors=0;
	Computed=1;
}


void CRanLip::PrepareHatFunctionAuto(int num, int numfine, double minLip)
{
	int i,j;

	num=max_(num,1);
	numfine=max_(numfine-1,1);

	for(bits=1;bits<32;bits++) if((1<<bits) >= numfine+1) break;
	numfine=(1<<bits) -1;
	mask1= (1<<(bits)) - 1;


	TotalElements=1;
	num_partition=num;   num_small_partition=numfine;  num_small_partition_p1=num_small_partition+1;
	for(i=0;i<Dimension;i++) TotalElements*=num; 

	totvals=1;
	for(i=0;i<Dimension;i++) totvals*=num_small_partition_p1; 

	Probabilities=(double*) malloc(TotalElements*sizeof(double));

// working arrays
	LipschitzH=(double*)malloc(Dimension*sizeof(double));
	m_delta=(int*)malloc(Dimension*sizeof(int));
    vals=(double*)malloc(totvals*sizeof(double));


	Volume=1;
	m_delta[Dimension-1]=1;

	for(i=0;i<Dimension;i++)
	{
		h[i]=(m_boundRight[i]-m_boundLeft[i])/num;
		hfine[i]=h[i]/numfine;
		Volume *= h[i];

		if(i>0) m_delta[Dimension-i-1] = m_delta[Dimension-i] * num_small_partition_p1;
	}

	Lipschitz=0;
	double d;

// make a special case if numfine==2, not to re-compute many values, store them in a table
// modify ComputeArray by using cache
	if(numfine<=1) {
		cache=(double*)malloc(TotalElements*sizeof(double));
		for(j=0;j<TotalElements;j++) {
			GetIJK(j);
			for(i=0;i<Dimension;i++) 
				m_tempLeft[i] = m_boundLeft[i] + h[i]*m_tempint[i];  // left top corner
			cache[j]=Distribution(m_tempLeft);
		}

		for(j=0;j<TotalElements;j++) {
			ComputeArrayCache(j);
			d= Dimension* max_( ComputeLipschitzBin(), minLip); // an extra factor
//			printf("%f\n",d);
			if(d>Lipschitz) Lipschitz=d;
			for(i=0;i<Dimension;i++) LipschitzH[i]=0.5 * d * hfine[i] *2 ; //will divide by 2 later

			d=ComputeMaxBin();
			Probabilities[j]=d;
		}
		free(cache);
	} // special case
	else { // general case

		for(j=0;j<TotalElements;j++) {
			GetIJK(j);
			for(i=0;i<Dimension;i++) {
				m_tempLeft[i] = m_boundLeft[i] + h[i]*m_tempint[i];  // left top corner
			}
			ComputeArray();
			d =Dimension * max_(ComputeLipschitzBin(),minLip);
			if(d>Lipschitz) Lipschitz=d;
			for(i=0;i<Dimension;i++) LipschitzH[i]=0.5 * d * hfine[i] * 2; //will divide by 2 later

			d=ComputeMaxBin();
			Probabilities[j]=d;
		}
	}
	Lipschitz=max_(Lipschitz, 1e-10);

	free(vals);
	free(LipschitzH);
	free(m_delta);

	for(j=0;j<TotalElements;j++) Probabilities[j] *= Volume;

	Dist = gsl_ran_discrete_preproc(TotalElements,Probabilities);

	// now get back to function values in probabilities
	for(j=0;j<TotalElements;j++) Probabilities[j] /= Volume;

	count_total=	count_errors=0;
	Computed=1;
}

void CRanLip::ComputeArray()
{
	int i,j,k,d1=Dimension-1;
	double v=hfine[d1];
	for(j=0;j<totvals;j++)
	{
		GetIJKfineBin(j);
		for(i=0;i<Dimension;i++) 
			V[i] = m_tempLeft[i] + hfine[i]*m_tempintfine[i];  

		vals[j]=Distribution(V);
		for(k=1;k<num_small_partition_p1;k++)
		{ //innermost index, no need for GetIJKfineBin
			j++;
			V[d1] += v;
			vals[j]=Distribution(V);
		}
	
	}
}


void CRanLip::ComputeArrayCache(int J)
{ // use previously computed values
	int j,i,k;
	GetIJK(J);
	vals[0]=cache[J];

	for(i=0;i<Dimension;i++) 
			m_tempLeft[i] = m_boundLeft[i] + h[i]*m_tempint[i];  // left top corner

	for(j=1;j<totvals;j++)
	{
		GetIJKfineBin(j);
		k=GetIndexfromIJK(m_tempintfine);

		if(k<TotalElements)
			vals[j]=cache[k]; 
		else { // was not computed
			for(i=0;i<Dimension;i++) 
				V[i] = m_tempLeft[i] + hfine[i]*m_tempintfine[i];  
			vals[j]=Distribution(V);
		}
	}
}

int CRanLip::GetIndexfromIJK( int* IJK)
{
	// I have m_tempint and m_tempintfine=IJK, so that the needed element has indices
	// m_tempint[i]+IJK[i]
	int k,i,j,n;
	k=0;
	n=1;
	for(i=Dimension-1;i>=0;i--) {
		j=m_tempint[i]+IJK[i];
		k += n * j;
		if(j>=num_partition) return TotalElements+1; // outside the computed array
		n *= num_partition;
	}
	return k;
}


double CRanLip::ComputeMaxBin()
{
	double d=-10e20;
	double u,v, D;
	int i,j;
	for(j=0;j<totvals;j++)
	{
		u=vals[j];
		GetIJKfineBin(j);

		for(i=0;i<Dimension;i++) {
			if(m_tempintfine[i] < num_small_partition)
			{
				v=vals[j+m_delta[i]];
				D=(u+v)+LipschitzH[i]; 
				d=max_(d,D);
			}
		}
	}

	return d/2;
}


double CRanLip::ComputeLipschitzBin()
{
	double d=-10e20;
	double u,v;
	int i,j;
	for(j=0;j<totvals;j++)
	{
		u=vals[j];
		GetIJKfineBin(j);

		for(i=0;i<Dimension;i++) {
			if(m_tempintfine[i] < num_small_partition)  
			 {
				v=vals[j+m_delta[i]];
				v=fabs(u-v)/hfine[i];
				d=max_(d,v);
			}
		}
	}

	return d;
}

void	CRanLip::GetIJK(int b)
{
	int i;
	div_t t;
	for(i=1;i<Dimension;i++) {
		t=div(b,num_partition);
		m_tempint[Dimension-i]=t.rem;
		b=t.quot;
	}
	m_tempint[0]=b;
}


void	CRanLip::GetIJKfineBin(int b)
{
	int i;
	for(i=Dimension-1;i>0;i--) {
		m_tempintfine[i] = b & mask1;//rem;
		b=b>>bits;//quot;
	}
	m_tempintfine[0]=b;
}



void	CRanLip::SavePartition(char * fname)
{
	if(!Computed) return;

	int j;
	FILE * fp=fopen(fname,"w");
	fprintf(fp,"Dim,Elements,Volume %d %d %d %f\n",Dimension,TotalElements,num_partition, Volume);

	for(j=0;j<Dimension;j++) {
		fprintf(fp,"%f %f\n",m_boundLeft[j],m_boundRight[j]);
	}

	for(j=0;j<TotalElements;j++) {
		fprintf(fp,"%f\n",Probabilities[j]);
	}
	fclose(fp);
}

void	CRanLip::LoadPartition(char * fname)
{
	FreeMem();


	int i,j;
	FILE * fp=fopen(fname,"r");
	fscanf(fp,"Dim,Elements,Volume %d %d %d %lf\n",&Dimension,&TotalElements, &num_partition, &Volume);

	double* t_boundLeft=(double*) malloc(Dimension*sizeof(double)); 
	double* t_boundRight=(double*) malloc(Dimension*sizeof(double));

	for(j=0;j<Dimension;j++) {
		fscanf(fp,"%lf %lf\n",&t_boundLeft[j], &t_boundRight[j]);
	}

	Init(Dimension,t_boundLeft,t_boundRight);
	free(t_boundLeft);
	free(t_boundRight);

	double r;
	Probabilities=(double*) calloc(TotalElements,sizeof(double));
	if(Probabilities==NULL) return;

	for(j=0;j<TotalElements;j++) {
		fscanf(fp,"%lf\n",&r);
		Probabilities[j] = r;
	}
	fclose(fp);

	for(i=0;i<Dimension;i++)
	{
		h[i]=(m_boundRight[i]-m_boundLeft[i])/num_partition;
	}

// for the discrete generator
	Dist = gsl_ran_discrete_preproc(TotalElements,Probabilities);

	count_total=	count_errors=0;
	Computed=1;
}



