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

#ifndef _apserv_h
#define _apserv_h

#include "aenv.h"
#include "ialglib.h"


/*$ Declarations $*/


/*************************************************************************
Buffers for internal functions which need buffers:
* check for size of the buffer you want to use.
* if buffer is too small, resize it; leave unchanged, if it is larger than
  needed.
* use it.

We can pass this structure to multiple functions;  after first run through
functions buffer sizes will be finally determined,  and  on  a next run no
allocation will be required.
*************************************************************************/
typedef struct
{
    ae_vector ia1;
    ae_vector ia2;
    ae_vector ra1;
    ae_vector ra2;
} apbuffers;


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void bvectorsetlengthatleast(/* Boolean */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rvectorsetlengthatleast(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
If Cols(X)<N or Rows(X)<M, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsetlengthatleast(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
This function checks that all values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitevector(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
This function checks that all values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitecvector(/* Complex */ ae_vector* z,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitematrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitecmatrix(/* Complex */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitertrmatrix(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state);


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitectrmatrix(/* Complex */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state);


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are  finite  or
NaN's.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfiniteornanmatrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag2(double x, double y, ae_state *_state);


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag3(double x, double y, double z, ae_state *_state);


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
ae_int_t saferdiv(double x, double y, double* r, ae_state *_state);


/*************************************************************************
This function calculates "safe" min(X/Y,V) for positive finite X, Y, V.
No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safeminposrv(double x, double y, double v, ae_state *_state);


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
     ae_state *_state);


/*************************************************************************
'bounds' value: maps X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
double boundval(double x, double b1, double b2, ae_state *_state);
ae_bool _apbuffers_init(apbuffers* p, ae_state *_state, ae_bool make_automatic);
ae_bool _apbuffers_init_copy(apbuffers* dst, apbuffers* src, ae_state *_state, ae_bool make_automatic);
void _apbuffers_clear(apbuffers* p);


/*$ End $*/
#endif

