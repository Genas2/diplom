/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
#include "apserv.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
This  function  generates  1-dimensional  general  interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1d(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            if( i!=n-1 )
            {
                x->ptr.p_double[i] = a+(i+0.2*(2*ae_randomreal(_state)-1))*h;
            }
            else
            {
                x->ptr.p_double[i] = b;
            }
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional equidistant interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dequidist(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            x->ptr.p_double[i] = a+i*h;
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*h;
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-1 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb1(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb1: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*(2*i+1)/(2*n), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-2 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb2(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb2: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*i/(n-1), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function checks that all values from X[] are distinct. It does more
than just usual floating point comparison:
* first, it calculates max(X) and min(X)
* second, it maps X[] from [min,max] to [1,2]
* only at this stage actual comparison is done

The meaning of such check is to ensure that all values are "distinct enough"
and will not cause interpolation subroutine to fail.

NOTE:
    X[] must be sorted by ascending (subroutine ASSERT's it)

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool aredistinct(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    double a;
    double b;
    ae_int_t i;
    ae_bool nonsorted;
    ae_bool result;


    ae_assert(n>=1, "APSERVAreDistinct: internal error (N<1)", _state);
    if( n==1 )
    {
        
        /*
         * everything is alright, it is up to caller to decide whether it
         * can interpolate something with just one point
         */
        result = ae_true;
        return result;
    }
    a = x->ptr.p_double[0];
    b = x->ptr.p_double[0];
    nonsorted = ae_false;
    for(i=1; i<=n-1; i++)
    {
        a = ae_minreal(a, x->ptr.p_double[i], _state);
        b = ae_maxreal(b, x->ptr.p_double[i], _state);
        nonsorted = nonsorted||ae_fp_greater_eq(x->ptr.p_double[i-1],x->ptr.p_double[i]);
    }
    ae_assert(!nonsorted, "APSERVAreDistinct: internal error (not sorted)", _state);
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_eq((x->ptr.p_double[i]-a)/(b-a)+1,(x->ptr.p_double[i-1]-a)/(b-a)+1) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void bvectorsetlengthatleast(/* Boolean */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rvectorsetlengthatleast(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Cols(X)<N or Rows(X)<M, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsetlengthatleast(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{


    if( x->rows<m||x->cols<n )
    {
        ae_matrix_set_length(x, m, n, _state);
    }
}


/*************************************************************************
This function checks that all values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitevector(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteVector: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(x->ptr.p_double[i], _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitecvector(/* Complex */ ae_vector* z,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCVector: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(z->ptr.p_complex[i].x, _state)||!ae_isfinite(z->ptr.p_complex[i].y, _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitematrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitecmatrix(/* Complex */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteCMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitertrmatrix(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteRTRMatrix: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitectrmatrix(/* Complex */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCTRMatrix: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are  finite  or
NaN's.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfiniteornanmatrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteOrNaNMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteOrNaNMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !(ae_isfinite(x->ptr.pp_double[i][j], _state)||ae_isnan(x->ptr.pp_double[i][j], _state)) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag2(double x, double y, ae_state *_state)
{
    double w;
    double xabs;
    double yabs;
    double z;
    double result;


    xabs = ae_fabs(x, _state);
    yabs = ae_fabs(y, _state);
    w = ae_maxreal(xabs, yabs, _state);
    z = ae_minreal(xabs, yabs, _state);
    if( ae_fp_eq(z,0) )
    {
        result = w;
    }
    else
    {
        result = w*ae_sqrt(1+ae_sqr(z/w, _state), _state);
    }
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag3(double x, double y, double z, ae_state *_state)
{
    double w;
    double result;


    w = ae_maxreal(ae_fabs(x, _state), ae_maxreal(ae_fabs(y, _state), ae_fabs(z, _state), _state), _state);
    if( ae_fp_eq(w,0) )
    {
        result = 0;
        return result;
    }
    x = x/w;
    y = y/w;
    z = z/w;
    result = w*ae_sqrt(ae_sqr(x, _state)+ae_sqr(y, _state)+ae_sqr(z, _state), _state);
    return result;
}


/*************************************************************************
Safe division.

This function attempts to calculate R=X/Y without overflow.

It returns:
* +1, if abs(X/Y)>=MaxRealNumber or undefined - overflow-like situation
      (no overlfow is generated, R is either NAN, PosINF, NegINF)
*  0, if MinRealNumber<abs(X/Y)<MaxRealNumber or X=0, Y<>0
      (R contains result, may be zero)
* -1, if 0<abs(X/Y)<MinRealNumber - underflow-like situation
      (R contains zero; it corresponds to underflow)

No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
ae_int_t saferdiv(double x, double y, double* r, ae_state *_state)
{
    ae_int_t result;

    *r = 0;

    
    /*
     * Two special cases:
     * * Y=0
     * * X=0 and Y<>0
     */
    if( ae_fp_eq(y,0) )
    {
        result = 1;
        if( ae_fp_eq(x,0) )
        {
            *r = _state->v_nan;
        }
        if( ae_fp_greater(x,0) )
        {
            *r = _state->v_posinf;
        }
        if( ae_fp_less(x,0) )
        {
            *r = _state->v_neginf;
        }
        return result;
    }
    if( ae_fp_eq(x,0) )
    {
        *r = 0;
        result = 0;
        return result;
    }
    
    /*
     * make Y>0
     */
    if( ae_fp_less(y,0) )
    {
        x = -x;
        y = -y;
    }
    
    /*
     *
     */
    if( ae_fp_greater_eq(y,1) )
    {
        *r = x/y;
        if( ae_fp_less_eq(ae_fabs(*r, _state),ae_minrealnumber) )
        {
            result = -1;
            *r = 0;
        }
        else
        {
            result = 0;
        }
    }
    else
    {
        if( ae_fp_greater_eq(ae_fabs(x, _state),ae_maxrealnumber*y) )
        {
            if( ae_fp_greater(x,0) )
            {
                *r = _state->v_posinf;
            }
            else
            {
                *r = _state->v_neginf;
            }
            result = 1;
        }
        else
        {
            *r = x/y;
            result = 0;
        }
    }
    return result;
}


/*************************************************************************
This function calculates "safe" min(X/Y,V) for positive finite X, Y, V.
No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safeminposrv(double x, double y, double v, ae_state *_state)
{
    double r;
    double result;


    if( ae_fp_greater_eq(y,1) )
    {
        
        /*
         * Y>=1, we can safely divide by Y
         */
        r = x/y;
        result = v;
        if( ae_fp_greater(v,r) )
        {
            result = r;
        }
        else
        {
            result = v;
        }
    }
    else
    {
        
        /*
         * Y<1, we can safely multiply by Y
         */
        if( ae_fp_less(x,v*y) )
        {
            result = x/y;
        }
        else
        {
            result = v;
        }
    }
    return result;
}


/*************************************************************************
This function makes periodic mapping of X to [A,B].

It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
such that X = T + K*(B-A).

NOTES:
* K is represented as real value, although actually it is integer
* T is guaranteed to be in [A,B]
* T replaces X

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
void apperiodicmap(double* x,
     double a,
     double b,
     double* k,
     ae_state *_state)
{

    *k = 0;

    ae_assert(ae_fp_less(a,b), "APPeriodicMap: internal error!", _state);
    *k = ae_ifloor((*x-a)/(b-a), _state);
    *x = *x-*k*(b-a);
    while(ae_fp_less(*x,a))
    {
        *x = *x+(b-a);
        *k = *k-1;
    }
    while(ae_fp_greater(*x,b))
    {
        *x = *x-(b-a);
        *k = *k+1;
    }
    *x = ae_maxreal(*x, a, _state);
    *x = ae_minreal(*x, b, _state);
}


/*************************************************************************
'bounds' value: maps X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
double boundval(double x, double b1, double b2, ae_state *_state)
{
    double result;


    if( ae_fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ae_fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


ae_bool _apbuffers_init(apbuffers* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->ia1, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ia2, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _apbuffers_init_copy(apbuffers* dst, apbuffers* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->ia1, &src->ia1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ia2, &src->ia2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra1, &src->ra1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra2, &src->ra2, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _apbuffers_clear(apbuffers* p)
{
    ae_vector_clear(&p->ia1);
    ae_vector_clear(&p->ia2);
    ae_vector_clear(&p->ra1);
    ae_vector_clear(&p->ra2);
}


/*$ End $*/
