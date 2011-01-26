/*************************************************************************
Copyright (c)
    2007, Sergey Bochkanov (ALGLIB project).
    1988, Pierre L'Ecuyer

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/


#include <stdafx.h>
#include "hqrnd.h"


/*$ Declarations $*/
static ae_int_t hqrnd_hqrndmax = 2147483563;
static ae_int_t hqrnd_hqrndm1 = 2147483563;
static ae_int_t hqrnd_hqrndm2 = 2147483399;
static ae_int_t hqrnd_hqrndmagic = 1634357784;
static ae_int_t hqrnd_hqrndintegerbase(hqrndstate* state,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndrandomize(hqrndstate* state, ae_state *_state)
{

    _hqrndstate_clear(state);

    hqrndseed(ae_randominteger(hqrnd_hqrndm1, _state), ae_randominteger(hqrnd_hqrndm2, _state), state, _state);
}


/*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndseed(ae_int_t s1,
     ae_int_t s2,
     hqrndstate* state,
     ae_state *_state)
{

    _hqrndstate_clear(state);

    state->s1 = s1%(hqrnd_hqrndm1-1)+1;
    state->s2 = s2%(hqrnd_hqrndm2-1)+1;
    state->v = (double)1/(double)hqrnd_hqrndmax;
    state->magicv = hqrnd_hqrndmagic;
}


/*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrnduniformr(hqrndstate* state, ae_state *_state)
{
    double result;


    result = state->v*hqrnd_hqrndintegerbase(state, _state);
    return result;
}


/*************************************************************************
This function generates random integer number in [0, N)

1. N must be less than HQRNDMax-1.
2. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t hqrnduniformi(hqrndstate* state, ae_int_t n, ae_state *_state)
{
    ae_int_t mx;
    ae_int_t result;


    
    /*
     * Correct handling of N's close to RNDBaseMax
     * (avoiding skewed distributions for RNDBaseMax<>K*N)
     */
    ae_assert(n>0, "HQRNDUniformI: N<=0!", _state);
    ae_assert(n<hqrnd_hqrndmax-1, "HQRNDUniformI: N>=RNDBaseMax-1!", _state);
    mx = hqrnd_hqrndmax-1-(hqrnd_hqrndmax-1)%n;
    do
    {
        result = hqrnd_hqrndintegerbase(state, _state)-1;
    }
    while(result>=mx);
    result = result%n;
    return result;
}


/*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrndnormal(hqrndstate* state, ae_state *_state)
{
    double v1;
    double v2;
    double result;


    hqrndnormal2(state, &v1, &v2, _state);
    result = v1;
    return result;
}


/*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndunit2(hqrndstate* state, double* x, double* y, ae_state *_state)
{
    double v;
    double mx;
    double mn;

    *x = 0;
    *y = 0;

    do
    {
        hqrndnormal2(state, x, y, _state);
    }
    while(!(ae_fp_neq(*x,0)||ae_fp_neq(*y,0)));
    mx = ae_maxreal(ae_fabs(*x, _state), ae_fabs(*y, _state), _state);
    mn = ae_minreal(ae_fabs(*x, _state), ae_fabs(*y, _state), _state);
    v = mx*ae_sqrt(1+ae_sqr(mn/mx, _state), _state);
    *x = *x/v;
    *y = *y/v;
}


/*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndnormal2(hqrndstate* state,
     double* x1,
     double* x2,
     ae_state *_state)
{
    double u;
    double v;
    double s;

    *x1 = 0;
    *x2 = 0;

    for(;;)
    {
        u = 2*hqrnduniformr(state, _state)-1;
        v = 2*hqrnduniformr(state, _state)-1;
        s = ae_sqr(u, _state)+ae_sqr(v, _state);
        if( ae_fp_greater(s,0)&&ae_fp_less(s,1) )
        {
            
            /*
             * two Sqrt's instead of one to
             * avoid overflow when S is too small
             */
            s = ae_sqrt(-2*ae_log(s, _state), _state)/ae_sqrt(s, _state);
            *x1 = u*s;
            *x2 = v*s;
            return;
        }
    }
}


/*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************/
double hqrndexponential(hqrndstate* state,
     double lambdav,
     ae_state *_state)
{
    double result;


    ae_assert(ae_fp_greater(lambdav,0), "HQRNDExponential: LambdaV<=0!", _state);
    result = -ae_log(hqrnduniformr(state, _state), _state)/lambdav;
    return result;
}


/*************************************************************************

L'Ecuyer, Efficient and portable combined random number generators
*************************************************************************/
static ae_int_t hqrnd_hqrndintegerbase(hqrndstate* state,
     ae_state *_state)
{
    ae_int_t k;
    ae_int_t result;


    ae_assert(state->magicv==hqrnd_hqrndmagic, "HQRNDIntegerBase: State is not correctly initialized!", _state);
    k = state->s1/53668;
    state->s1 = 40014*(state->s1-k*53668)-k*12211;
    if( state->s1<0 )
    {
        state->s1 = state->s1+2147483563;
    }
    k = state->s2/52774;
    state->s2 = 40692*(state->s2-k*52774)-k*3791;
    if( state->s2<0 )
    {
        state->s2 = state->s2+2147483399;
    }
    
    /*
     * Result
     */
    result = state->s1-state->s2;
    if( result<1 )
    {
        result = result+2147483562;
    }
    return result;
}


ae_bool _hqrndstate_init(hqrndstate* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _hqrndstate_init_copy(hqrndstate* dst, hqrndstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->s1 = src->s1;
    dst->s2 = src->s2;
    dst->v = src->v;
    dst->magicv = src->magicv;
    return ae_true;
}


void _hqrndstate_clear(hqrndstate* p)
{
}


/*$ End $*/
