/* randist/gsl_randist.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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

#ifndef __GSL_RANDIST_H__
#define __GSL_RANDIST_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS

#ifdef __cplusplus1
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct {                /* struct for Walker algorithm */
    size_t K;
    size_t *A;
    double *F;
} gsl_ran_discrete_t;

gsl_ran_discrete_t * gsl_ran_discrete_preproc (size_t K, const double *P);
void gsl_ran_discrete_free(gsl_ran_discrete_t *g);
size_t gsl_ran_discrete (/*const gsl_rng *r,*/ const gsl_ran_discrete_t *g);
double gsl_ran_discrete_pdf (size_t k, const gsl_ran_discrete_t *g);



static  unsigned long int ranlux_get (void *vstate);
static  double ranlux_get_double (void *vstate);
static  void ranlux_set_lux (void *state, unsigned long int s, unsigned int luxury);
static  void ranlux_set (void *state, unsigned long int s);
static  void ranlux389_set (void *state, unsigned long int s);

static const unsigned long int mask_lo = 0x00ffffffUL;	/* 2^24 - 1 */
static const unsigned long int mask_hi = ~0x00ffffffUL;
static const unsigned long int two24 = 16777216;	/* 2^24 */

typedef struct
  {
    unsigned int i;
    unsigned int j;
    unsigned int n;
    unsigned int skip;
    unsigned int carry;
    unsigned long int u[24];
  }
ranlux_state_t;

void ranlux_set_seed (unsigned long int s);
static double  UnifRand() { return (double)rand()/(RAND_MAX+0.001); };
double ranlux_get_double_V () ;

__END_DECLS

#endif /* __GSL_RANDIST_H__ */
