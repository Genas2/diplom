/*************************************************************************
Copyright (c) 2007-2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _idwint_h
#define _idwint_h

#include "aenv.h"
#include "ialglib.h"
#include "tsort.h"
#include "apserv.h"
#include "nearestneighbor.h"
#include "reflections.h"
#include "hblas.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"


/*$ Declarations $*/


/*************************************************************************
IDW interpolant.
*************************************************************************/
typedef struct
{
    ae_int_t n;
    ae_int_t nx;
    ae_int_t d;
    double r;
    ae_int_t nw;
    kdtree tree;
    ae_int_t modeltype;
    ae_matrix q;
    ae_vector xbuf;
    ae_vector tbuf;
    ae_vector rbuf;
    ae_matrix xybuf;
    ae_int_t debugsolverfailures;
    double debugworstrcond;
    double debugbestrcond;
} idwinterpolant;


/*$ Body $*/


/*************************************************************************
IDW interpolation

INPUT PARAMETERS:
    Z   -   IDW interpolant built with one of model building
            subroutines.
    X   -   array[0..NX-1], interpolation point

Result:
    IDW interpolant Z(X)

  -- ALGLIB --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************/
double idwcalc(idwinterpolant* z,
     /* Real    */ ae_vector* x,
     ae_state *_state);


/*************************************************************************
IDW interpolant using modified Shepard method for uniform point
distributions.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    D   -   nodal function type, either:
            * 0     constant  model.  Just  for  demonstration only, worst
                    model ever.
            * 1     linear model, least squares fitting. Simpe  model  for
                    datasets too small for quadratic models
            * 2     quadratic  model,  least  squares  fitting. Best model
                    available (if your dataset is large enough).
            * -1    "fast"  linear  model,  use  with  caution!!!   It  is
                    significantly  faster than linear/quadratic and better
                    than constant model. But it is less robust (especially
                    in the presence of noise).
    NQ  -   number of points used to calculate  nodal  functions  (ignored
            for constant models). NQ should be LARGER than:
            * max(1.5*(1+NX),2^NX+1) for linear model,
            * max(3/4*(NX+2)*(NX+1),2^NX+1) for quadratic model.
            Values less than this threshold will be silently increased.
    NW  -   number of points used to calculate weights and to interpolate.
            Required: >=2^NX+1, values less than this  threshold  will  be
            silently increased.
            Recommended value: about 2*NQ

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.
    
NOTES:
  * best results are obtained with quadratic models, worst - with constant
    models
  * when N is large, NQ and NW must be significantly smaller than  N  both
    to obtain optimal performance and to obtain optimal accuracy. In 2  or
    3-dimensional tasks NQ=15 and NW=25 are good values to start with.
  * NQ  and  NW  may  be  greater  than  N.  In  such  cases  they will be
    automatically decreased.
  * this subroutine is always succeeds (as long as correct parameters  are
    passed).
  * see  'Multivariate  Interpolation  of Large Sets of Scattered Data' by
    Robert J. Renka for more information on this algorithm.
  * this subroutine assumes that point distribution is uniform at the small
    scales.  If  it  isn't  -  for  example,  points are concentrated along
    "lines", but "lines" distribution is uniform at the larger scale - then
    you should use IDWBuildModifiedShepardR()


  -- ALGLIB PROJECT --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************/
void idwbuildmodifiedshepard(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t d,
     ae_int_t nq,
     ae_int_t nw,
     idwinterpolant* z,
     ae_state *_state);


/*************************************************************************
IDW interpolant using modified Shepard method for non-uniform datasets.

This type of model uses  constant  nodal  functions and interpolates using
all nodes which are closer than user-specified radius R. It  may  be  used
when points distribution is non-uniform at the small scale, but it  is  at
the distances as large as R.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    R   -   radius, R>0

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.

NOTES:
* if there is less than IDWKMin points within  R-ball,  algorithm  selects
  IDWKMin closest ones, so that continuity properties of  interpolant  are
  preserved even far from points.

  -- ALGLIB PROJECT --
     Copyright 11.04.2010 by Bochkanov Sergey
*************************************************************************/
void idwbuildmodifiedshepardr(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     double r,
     idwinterpolant* z,
     ae_state *_state);


/*************************************************************************
IDW model for noisy data.

This subroutine may be used to handle noisy data, i.e. data with noise  in
OUTPUT values.  It differs from IDWBuildModifiedShepard() in the following
aspects:
* nodal functions are not constrained to pass through  nodes:  Qi(xi)<>yi,
  i.e. we have fitting  instead  of  interpolation.
* weights which are used during least  squares fitting stage are all equal
  to 1.0 (independently of distance)
* "fast"-linear or constant nodal functions are not supported (either  not
  robust enough or too rigid)

This problem require far more complex tuning than interpolation  problems.
Below you can find some recommendations regarding this problem:
* focus on tuning NQ; it controls noise reduction. As for NW, you can just
  make it equal to 2*NQ.
* you can use cross-validation to determine optimal NQ.
* optimal NQ is a result of complex tradeoff  between  noise  level  (more
  noise = larger NQ required) and underlying  function  complexity  (given
  fixed N, larger NQ means smoothing of compex features in the data).  For
  example, NQ=N will reduce noise to the minimum level possible,  but  you
  will end up with just constant/linear/quadratic (depending on  D)  least
  squares model for the whole dataset.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    D   -   nodal function degree, either:
            * 1     linear model, least squares fitting. Simpe  model  for
                    datasets too small for quadratic models (or  for  very
                    noisy problems).
            * 2     quadratic  model,  least  squares  fitting. Best model
                    available (if your dataset is large enough).
    NQ  -   number of points used to calculate nodal functions.  NQ should
            be  significantly   larger   than  1.5  times  the  number  of
            coefficients in a nodal function to overcome effects of noise:
            * larger than 1.5*(1+NX) for linear model,
            * larger than 3/4*(NX+2)*(NX+1) for quadratic model.
            Values less than this threshold will be silently increased.
    NW  -   number of points used to calculate weights and to interpolate.
            Required: >=2^NX+1, values less than this  threshold  will  be
            silently increased.
            Recommended value: about 2*NQ or larger

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.

NOTES:
  * best results are obtained with quadratic models, linear models are not
    recommended to use unless you are pretty sure that it is what you want
  * this subroutine is always succeeds (as long as correct parameters  are
    passed).
  * see  'Multivariate  Interpolation  of Large Sets of Scattered Data' by
    Robert J. Renka for more information on this algorithm.


  -- ALGLIB PROJECT --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************/
void idwbuildnoisy(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t d,
     ae_int_t nq,
     ae_int_t nw,
     idwinterpolant* z,
     ae_state *_state);
ae_bool _idwinterpolant_init(idwinterpolant* p, ae_state *_state, ae_bool make_automatic);
ae_bool _idwinterpolant_init_copy(idwinterpolant* dst, idwinterpolant* src, ae_state *_state, ae_bool make_automatic);
void _idwinterpolant_clear(idwinterpolant* p);


/*$ End $*/
#endif

