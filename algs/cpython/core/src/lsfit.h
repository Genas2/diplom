/*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _lsfit_h
#define _lsfit_h

#include "aenv.h"
#include "ialglib.h"
#include "tsort.h"
#include "ratint.h"
#include "apserv.h"
#include "polint.h"
#include "spline1d.h"
#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "linmin.h"
#include "fbls.h"
#include "minlbfgs.h"
#include "minlm.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"


/*$ Declarations $*/


/*************************************************************************
Polynomial fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************/
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} polynomialfitreport;


/*************************************************************************
Barycentric fitting report:
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
    TaskRCond       reciprocal of task's condition number
*************************************************************************/
typedef struct
{
    double taskrcond;
    ae_int_t dbest;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} barycentricfitreport;


/*************************************************************************
Spline fitting report:
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
    
Fields  below are  filled  by   obsolete    functions   (Spline1DFitCubic,
Spline1DFitHermite). Modern fitting functions do NOT fill these fields:
    TaskRCond       reciprocal of task's condition number
*************************************************************************/
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} spline1dfitreport;


/*************************************************************************
Least squares fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************/
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} lsfitreport;


/*************************************************************************
Nonlinear fitter.

You should use ALGLIB functions to work with fitter.
Never try to access its fields directly!
*************************************************************************/
typedef struct
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t k;
    double epsf;
    double epsx;
    ae_int_t maxits;
    double stpmax;
    ae_bool xrep;
    ae_matrix taskx;
    ae_vector tasky;
    ae_vector w;
    ae_bool xupdated;
    ae_bool needf;
    ae_bool needfg;
    ae_bool needfgh;
    ae_int_t pointindex;
    ae_vector x;
    ae_vector c;
    double f;
    ae_vector g;
    ae_matrix h;
    ae_int_t repterminationtype;
    double reprmserror;
    double repavgerror;
    double repavgrelerror;
    double repmaxerror;
    minlmstate optstate;
    minlmreport optrep;
    rcommstate rstate;
} lsfitstate;


/*$ Body $*/


/*************************************************************************
Fitting by polynomials in barycentric form. This function provides  simple
unterface for unconstrained unweighted fitting. See  PolynomialFitWC()  if
you need constrained fitting.

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO:
    PolynomialFitWC()

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0
            * if given, only leading N elements of X/Y are used
            * if not given, automatically determined from sizes of X/Y
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
    P   -   interpolant in barycentric form.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

NOTES:
    you can convert P from barycentric form  to  the  power  or  Chebyshev
    basis with PolynomialBar2Pow() or PolynomialBar2Cheb() functions  from
    POLINT subpackage.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************/
void polynomialfit(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* p,
     polynomialfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted  fitting by polynomials in barycentric form, with constraints  on
function values or first derivatives.

Small regularizing term is used when solving constrained tasks (to improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO:
    PolynomialFit()

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
            * if given, only leading N elements of X/Y/W are used
            * if not given, automatically determined from sizes of X/Y/W
    XC  -   points where polynomial values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that P(XC[i])=YC[i]
            * DC[i]=1   means that P'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    P   -   interpolant in barycentric form.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

NOTES:
    you can convert P from barycentric form  to  the  power  or  Chebyshev
    basis with PolynomialBar2Pow() or PolynomialBar2Cheb() functions  from
    POLINT subpackage.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* even simple constraints can be inconsistent, see  Wikipedia  article  on
  this subject: http://en.wikipedia.org/wiki/Birkhoff_interpolation
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the one special cases, however, we can  guarantee  consistency.  This
  case  is:  M>1  and constraints on the function values (NOT DERIVATIVES)

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************/
void polynomialfitwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* p,
     polynomialfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weghted rational least  squares  fitting  using  Floater-Hormann  rational
functions  with  optimal  D  chosen  from  [0,9],  with  constraints   and
individual weights.

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal D (least WEIGHTED root
mean square error) is chosen.  Task  is  linear,  so  linear least squares
solver  is  used.  Complexity  of  this  computational  scheme is O(N*M^2)
(mostly dominated by the least squares solver).

SEE ALSO
* BarycentricFitFloaterHormann(), "lightweight" fitting without invididual
  weights and constraints.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where function values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroutine doesn't calculate task's condition number for K<>0.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained barycentric interpolants:
* excessive  constraints  can  be  inconsistent.   Floater-Hormann   basis
  functions aren't as flexible as splines (although they are very smooth).
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints on the function  VALUES at the interval
  boundaries. Note that consustency of the  constraints  on  the  function
  DERIVATIVES is NOT guaranteed (you can use in such cases  cubic  splines
  which are more flexible).
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormannwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state);


/*************************************************************************
Rational least squares fitting using  Floater-Hormann  rational  functions
with optimal D chosen from [0,9].

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal  D  (least  root  mean
square error) is chosen.  Task  is  linear, so linear least squares solver
is used. Complexity  of  this  computational  scheme is  O(N*M^2)  (mostly
dominated by the least squares solver).

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0.
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormann(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state);


/*************************************************************************
Rational least squares fitting using  Floater-Hormann  rational  functions
with optimal D chosen from [0,9].

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal  D  (least  root  mean
square error) is chosen.  Task  is  linear, so linear least squares solver
is used. Complexity  of  this  computational  scheme is  O(N*M^2)  (mostly
dominated by the least squares solver).

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0.
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitpenalized(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     double rho,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted fitting by penalized cubic spline.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are cubic splines with  natural  boundary
conditions. Problem is regularized by  adding non-linearity penalty to the
usual least squares penalty function:

    S(x) = arg min { LS + P }, where
    LS   = SUM { w[i]^2*(y[i] - S(x[i]))^2 } - least squares penalty
    P    = C*10^rho*integral{ S''(x)^2*dx } - non-linearity penalty
    rho  - tunable constant given by user
    C    - automatically determined scale parameter,
           makes penalty invariant with respect to scaling of X, Y, W.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            problem.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    M   -   number of basis functions ( = number_of_nodes), M>=4.
    Rho -   regularization  constant  passed   by   user.   It   penalizes
            nonlinearity in the regression spline. It  is  logarithmically
            scaled,  i.e.  actual  value  of  regularization  constant  is
            calculated as 10^Rho. It is automatically scaled so that:
            * Rho=2.0 corresponds to moderate amount of nonlinearity
            * generally, it should be somewhere in the [-8.0,+8.0]
            If you do not want to penalize nonlineary,
            pass small Rho. Values as low as -15 should work.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD or
                           Cholesky decomposition; problem may be
                           too ill-conditioned (very rare)
    S   -   spline interpolant.
    Rep -   Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

NOTE 1: additional nodes are added to the spline outside  of  the  fitting
interval to force linearity when x<min(x,xc) or x>max(x,xc).  It  is  done
for consistency - we penalize non-linearity  at [min(x,xc),max(x,xc)],  so
it is natural to force linearity outside of this interval.

NOTE 2: function automatically sorts points,  so  caller may pass unsorted
array.

  -- ALGLIB PROJECT --
     Copyright 19.10.2010 by Bochkanov Sergey
*************************************************************************/
void spline1dfitpenalizedw(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     ae_int_t m,
     double rho,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted fitting by cubic  spline,  with constraints on function values or
derivatives.

Equidistant grid with M-2 nodes on [min(x,xc),max(x,xc)] is  used to build
basis functions. Basis functions are cubic splines with continuous  second
derivatives  and  non-fixed first  derivatives  at  interval  ends.  Small
regularizing term is used  when  solving  constrained  tasks  (to  improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitHermiteWC()  -   fitting by Hermite splines (more flexible,
                                less smooth)
    Spline1DFitCubic()      -   "lightweight" fitting  by  cubic  splines,
                                without invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints (optional):
            * 0<=K<M.
            * K=0 means no constraints (XC/YC/DC are not used)
            * if given, only first K elements of XC/YC/DC are used
            * if not given, automatically determined from XC/YC/DC
    M   -   number of basis functions ( = number_of_nodes+2), M>=4.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints  on  the  function  values  AND/OR  its
  derivatives at the interval boundaries.
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.


  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubicwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted  fitting  by Hermite spline,  with constraints on function values
or first derivatives.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are Hermite splines.  Small  regularizing
term is used when solving constrained tasks (to improve stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitCubicWC()    -   fitting by Cubic splines (less flexible,
                                more smooth)
    Spline1DFitHermite()    -   "lightweight" Hermite fitting, without
                                invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints (optional):
            * 0<=K<M.
            * K=0 means no constraints (XC/YC/DC are not used)
            * if given, only first K elements of XC/YC/DC are used
            * if not given, automatically determined from XC/YC/DC
    M   -   number of basis functions (= 2 * number of nodes),
            M>=4,
            M IS EVEN!

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -2 means odd M was passed (which is not supported)
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

IMPORTANT:
    this subroitine supports only even M's


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the several special cases, however, we can guarantee consistency.
* one of this cases is  M>=4  and   constraints  on   the  function  value
  (AND/OR its derivative) at the interval boundaries.
* another special case is M>=4  and  ONE  constraint on the function value
  (OR, BUT NOT AND, derivative) anywhere in [min(x),max(x)]

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermitewc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Least squares fitting by cubic spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitCubicWC().  See  Spline1DFitCubicWC() for more information
about subroutine parameters (we don't duplicate it here because of length)

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubic(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Least squares fitting by Hermite spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitHermiteWC().  See Spline1DFitHermiteWC()  description  for
more information about subroutine parameters (we don't duplicate  it  here
because of length).

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermite(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -1    incorrect N/M were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearw(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted constained linear least squares fitting.

This  is  variation  of LSFitLinearW(), which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then LSFitLinearW()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearwc(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);


/*************************************************************************
Linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinear(/* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);


/*************************************************************************
Constained linear least squares fitting.

This  is  variation  of LSFitLinear(),  which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then  LSFitLinear()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearc(/* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);


/*************************************************************************
Weighted nonlinear least squares fitting using function values only.

Combination of numerical differentiation and secant updates is used to
obtain function Jacobian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]).

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    DiffStep-   numerical differentiation step;
                should not be very small or large;
                large = loss of accuracy
                small = growth of round-off errors

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 18.10.2008 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewf(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     double diffstep,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Nonlinear least squares fitting using function values only.

Combination of numerical differentiation and secant updates is used to
obtain function Jacobian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (f(c,x[0])-y[0])^2 + ... + (f(c,x[n-1])-y[n-1])^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]).

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    DiffStep-   numerical differentiation step;
                should not be very small or large;
                large = loss of accuracy
                small = growth of round-off errors

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 18.10.2008 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatef(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     double diffstep,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Weighted nonlinear least squares fitting using gradient only.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,
    
    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted
    
This subroutine uses only f(c,x[i]) and its gradient.
    
INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also:
    LSFitResults
    LSFitCreateFG (fitting without weights)
    LSFitCreateWFGH (fitting using Hessian)
    LSFitCreateFGH (fitting using Hessian, without weights)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewfg(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_bool cheapfg,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Nonlinear least squares fitting using gradient only, without individual
weights.

Nonlinear task min(F(c)) is solved, where

    F(c) = ((f(c,x[0])-y[0]))^2 + ... + ((f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]) and its gradient.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatefg(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_bool cheapfg,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Weighted nonlinear least squares fitting using gradient/Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(c,x[i]), its gradient and its Hessian.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewfgh(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Nonlinear least squares fitting using gradient/Hessian, without individial
weights.

Nonlinear task min(F(c)) is solved, where

    F(c) = ((f(c,x[0])-y[0]))^2 + ... + ((f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(c,x[i]), its gradient and its Hessian.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatefgh(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     lsfitstate* state,
     ae_state *_state);


/*************************************************************************
Stopping conditions for nonlinear least squares fitting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (according to the scheme used by MINLM unit).


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitsetcond(lsfitstate* state,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void lsfitsetstpmax(lsfitstate* state, double stpmax, ae_state *_state);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not
    
When reports are needed, State.C (current parameters) and State.F (current
value of fitting function) are reported.


  -- ALGLIB --
     Copyright 15.08.2010 by Bochkanov Sergey
*************************************************************************/
void lsfitsetxrep(lsfitstate* state, ae_bool needxrep, ae_state *_state);


/*************************************************************************
NOTES:

1. this algorithm is somewhat unusual because it works with  parameterized
   function f(C,X), where X is a function argument (we  have  many  points
   which are characterized by different  argument  values),  and  C  is  a
   parameter to fit.

   For example, if we want to do linear fit by f(c0,c1,x) = c0*x+c1,  then
   x will be argument, and {c0,c1} will be parameters.
   
   It is important to understand that this algorithm finds minimum in  the
   space of function PARAMETERS (not arguments), so it  needs  derivatives
   of f() with respect to C, not X.
   
   In the example above it will need f=c0*x+c1 and {df/dc0,df/dc1} = {x,1}
   instead of {df/dx} = {c0}.

2. Callback functions accept C as the first parameter, and X as the second

3. If  state  was  created  with  LSFitCreateFG(),  algorithm  needs  just
   function   and   its   gradient,   but   if   state   was  created with
   LSFitCreateFGH(), algorithm will need function, gradient and Hessian.
   
   According  to  the  said  above,  there  ase  several  versions of this
   function, which accept different sets of callbacks.
   
   This flexibility opens way to subtle errors - you may create state with
   LSFitCreateFGH() (optimization using Hessian), but call function  which
   does not accept Hessian. So when algorithm will request Hessian,  there
   will be no callback to call. In this case exception will be thrown.
   
   Be careful to avoid such errors because there is no way to find them at
   compile time - you can see them at runtime only.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool lsfititeration(lsfitstate* state, ae_state *_state);


/*************************************************************************
Nonlinear least squares fitting results.

Called after return from LSFitFit().

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    Info    -   completetion code:
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
    C       -   array[0..K-1], solution
    Rep     -   optimization report. Following fields are set:
                * Rep.TerminationType completetion code:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitresults(lsfitstate* state,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);


/*************************************************************************
Internal subroutine: automatic scaling for LLS tasks.
NEVER CALL IT DIRECTLY!

Maps abscissas to [-1,1], standartizes ordinates and correspondingly scales
constraints. It also scales weights so that max(W[i])=1

Transformations performed:
* X, XC         [XA,XB] => [-1,+1]
                transformation makes min(X)=-1, max(X)=+1

* Y             [SA,SB] => [0,1]
                transformation makes mean(Y)=0, stddev(Y)=1
                
* YC            transformed accordingly to SA, SB, DC[I]

  -- ALGLIB PROJECT --
     Copyright 08.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitscalexy(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     double* xa,
     double* xb,
     double* sa,
     double* sb,
     /* Real    */ ae_vector* xoriginal,
     /* Real    */ ae_vector* yoriginal,
     ae_state *_state);
ae_bool _polynomialfitreport_init(polynomialfitreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _polynomialfitreport_init_copy(polynomialfitreport* dst, polynomialfitreport* src, ae_state *_state, ae_bool make_automatic);
void _polynomialfitreport_clear(polynomialfitreport* p);
ae_bool _barycentricfitreport_init(barycentricfitreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _barycentricfitreport_init_copy(barycentricfitreport* dst, barycentricfitreport* src, ae_state *_state, ae_bool make_automatic);
void _barycentricfitreport_clear(barycentricfitreport* p);
ae_bool _spline1dfitreport_init(spline1dfitreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _spline1dfitreport_init_copy(spline1dfitreport* dst, spline1dfitreport* src, ae_state *_state, ae_bool make_automatic);
void _spline1dfitreport_clear(spline1dfitreport* p);
ae_bool _lsfitreport_init(lsfitreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _lsfitreport_init_copy(lsfitreport* dst, lsfitreport* src, ae_state *_state, ae_bool make_automatic);
void _lsfitreport_clear(lsfitreport* p);
ae_bool _lsfitstate_init(lsfitstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _lsfitstate_init_copy(lsfitstate* dst, lsfitstate* src, ae_state *_state, ae_bool make_automatic);
void _lsfitstate_clear(lsfitstate* p);


/*$ End $*/
#endif

