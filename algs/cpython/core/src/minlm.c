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
#include "minlm.h"


/*$ Declarations $*/
static ae_int_t minlm_lmmodefj = 0;
static ae_int_t minlm_lmmodefgj = 1;
static ae_int_t minlm_lmmodefgh = 2;
static ae_int_t minlm_lmflagnoprelbfgs = 1;
static ae_int_t minlm_lmflagnointlbfgs = 2;
static ae_int_t minlm_lmprelbfgsm = 5;
static ae_int_t minlm_lmintlbfgsits = 5;
static ae_int_t minlm_lbfgsnorealloc = 1;
static double minlm_lambdaup = 2.0;
static double minlm_lambdadown = 0.33;
static double minlm_suspiciousnu = 16;
static ae_int_t minlm_smallmodelage = 3;
static ae_int_t minlm_additers = 5;
static void minlm_lmprepare(ae_int_t n,
     ae_int_t m,
     ae_bool havegrad,
     minlmstate* state,
     ae_state *_state);
static void minlm_clearrequestfields(minlmstate* state, ae_state *_state);
static ae_bool minlm_increaselambda(double* lambdav,
     double* nu,
     ae_state *_state);
static void minlm_decreaselambda(double* lambdav,
     double* nu,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] and Jacobian of f[].


REQUIREMENTS:
This algorithm will request following information during its operation:

* function vector f[] at given point X
* function vector f[] and Jacobian of f[] (simultaneously) at given point

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec()  and jac() callbacks.
First  one  is used to calculate f[] at given point, second one calculates
f[] and Jacobian df[i]/dx[j].

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works  with  general  form function and does not provide Jacobian), but it
will  lead  to  exception  being  thrown  after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateVJ() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateVJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateVJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateVJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateVJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize, check parameters
     */
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_false;
    state->hasfi = ae_true;
    state->hasg = ae_false;
    
    /*
     * second stage of initialization
     */
    minlm_lmprepare(n, m, ae_false, state, _state);
    minlmsetacctype(state, 0, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] only. Finite differences  are  used  to
calculate Jacobian.


REQUIREMENTS:
This algorithm will request following information during its operation:
* function vector f[] at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec() callback.

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works with general form function and does not accept function vector), but
it will  lead  to  exception being thrown after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateV() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also MinLMIteration, MinLMResults.

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatev(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     double diffstep,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(ae_isfinite(diffstep, _state), "MinLMCreateV: DiffStep is not finite!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "MinLMCreateV: DiffStep<=0!", _state);
    ae_assert(n>=1, "MinLMCreateV: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateV: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateV: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateV: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->n = n;
    state->m = m;
    state->algomode = 0;
    state->hasf = ae_false;
    state->hasfi = ae_true;
    state->hasg = ae_false;
    state->diffstep = diffstep;
    
    /*
     * second stage of initialization
     */
    minlm_lmprepare(n, m, ae_false, state, _state);
    minlmsetacctype(state, 1, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

DESCRIPTION:
This  function  is  used  to  find  minimum  of general form (not "sum-of-
-squares") function
    F = F(x[0], ..., x[n-1])
using  its  gradient  and  Hessian.  Levenberg-Marquardt modification with
L-BFGS pre-optimization and internal pre-conditioned  L-BFGS  optimization
after each Levenberg-Marquardt step is used.


REQUIREMENTS:
This algorithm will request following information during its operation:

* function value F at given point X
* F and gradient G (simultaneously) at given point X
* F, G and Hessian H (simultaneously) at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts func(),  grad()  and  hess()
function pointers. First pointer is used to calculate F  at  given  point,
second  one  calculates  F(x)  and  grad F(x),  third one calculates F(x),
grad F(x), hess F(x).

You can try to initialize MinLMState structure with FGH-function and  then
use incorrect version of MinLMOptimize() (for example, version which  does
not provide Hessian matrix), but it will lead to  exception  being  thrown
after first attempt to calculate Hessian.


USAGE:
1. User initializes algorithm state with MinLMCreateFGH() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   pointers (delegates, etc.) to callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgh(ae_int_t n,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateFGH: N<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateFGH: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateFGH: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->n = n;
    state->m = 0;
    state->algomode = 2;
    state->hasf = ae_true;
    state->hasfi = ae_false;
    state->hasg = ae_true;
    
    /*
     * init2
     */
    minlm_lmprepare(n, 0, ae_true, state, _state);
    minlmsetacctype(state, 2, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using:
* value of function vector f[]
* value of Jacobian of f[]
* gradient of merit function F(x)

This function creates optimizer which uses acceleration strategy 2.  Cheap
gradient of merit function (which is twice the product of function  vector
and Jacobian) is used for accelerated iterations (see User Guide for  more
info on this subject).

REQUIREMENTS:
This algorithm will request following information during its operation:

* function vector f[] at given point X
* function vector f[] and Jacobian of f[] (simultaneously) at given point
* gradient of

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts  fvec(),  jac()  and  grad()
callbacks. First one is used to calculate f[] at given point,  second  one
calculates f[] and Jacobian df[i]/dx[j], last one calculates  gradient  of
merit function F(x).

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works  with  general  form function and does not provide Jacobian), but it
will  lead  to  exception  being  thrown  after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateVGJ() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevgj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateVGJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateVGJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateVGJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateVGJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize, check parameters
     */
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_false;
    state->hasfi = ae_true;
    state->hasg = ae_false;
    
    /*
     * second stage of initialization
     */
    minlm_lmprepare(n, m, ae_false, state, _state);
    minlmsetacctype(state, 2, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
                   LEVENBERG-MARQUARDT-LIKE METHOD FOR
                  NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:

This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of F(), gradient of F(), function vector f[]  and  Jacobian of
f[].

This function is considered obsolete since ALGLIB 3.1.0 and is present for
backward  compatibility  only.  We  recommend to use MinLMCreateVGJ, which
provides similar, but more consistent interface.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateFGJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateFGJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateFGJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateFGJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_true;
    state->hasfi = ae_false;
    state->hasg = ae_true;
    
    /*
     * init2
     */
    minlm_lmprepare(n, m, ae_true, state, _state);
    minlmsetacctype(state, 2, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
    CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using  value  of  F(),  function  vector  f[] and Jacobian of f[]. Classic
Levenberg-Marquardt method is used.

This function is considered obsolete since ALGLIB 3.1.0 and is present for
backward  compatibility  only.  We  recommend  to use MinLMCreateVJ, which
provides similar, but more consistent and feature-rich interface.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateFJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateFJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateFJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateFJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_true;
    state->hasfi = ae_false;
    state->hasg = ae_false;
    
    /*
     * init 2
     */
    minlm_lmprepare(n, m, ae_true, state, _state);
    minlmsetacctype(state, 0, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
This function sets stopping conditions for Levenberg-Marquardt optimization
algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                ||G||<EpsG is satisfied, where ||.|| means Euclidian norm,
                G - gradient.
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations   is    unlimited.   Only   Levenberg-Marquardt
                iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                counted because their cost is very low compared to that of
                LM).

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetcond(minlmstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinLMSetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinLMSetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinLMSetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinLMSetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinLMSetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinLMSetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinLMSetCond: negative MaxIts!", _state);
    if( ((ae_fp_eq(epsg,0)&&ae_fp_eq(epsf,0))&&ae_fp_eq(epsx,0))&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state->epsg = epsg;
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinLMOptimize(). Both Levenberg-Marquardt and internal  L-BFGS
iterations are reported.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetxrep(minlmstate* state, ae_bool needxrep, ae_state *_state)
{


    state->xrep = needxrep;
}


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
void minlmsetstpmax(minlmstate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinLMSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinLMSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
This function is used to change acceleration settings

You can choose between three acceleration strategies:
* AccType=0, no acceleration.
* AccType=1, secant updates are used to update quadratic model after  each
  iteration. After fixed number of iterations (or after  model  breakdown)
  we  recalculate  quadratic  model  using  analytic  Jacobian  or  finite
  differences. Number of secant-based iterations depends  on  optimization
  settings: about 3 iterations - when we have analytic Jacobian, up to 2*N
  iterations - when we use finite differences to calculate Jacobian.
* AccType=2, after quadratic model is built and LM step is made, we use it
  as preconditioner for several (5-10) iterations of L-BFGS algorithm.

AccType=1 is recommended when Jacobian  calculation  cost  is  prohibitive
high (several Mx1 function vector calculations  followed  by  several  NxN
Cholesky factorizations are faster than calculation of one M*N  Jacobian).
It should also be used when we have no Jacobian, because finite difference
approximation takes too much time to compute.

AccType=2 is recommended when Jacobian is cheap - much more  cheaper  than
one  Cholesky  factorization.   We   can   reduce   number   of   Cholesky
factorizations at the cost of increased number of  Jacobian  calculations.
Sometimes it helps.

Table below list  optimization  protocols  (XYZ  protocol  corresponds  to
MinLMCreateXYZ) and acceleration types they support (and use by  default).

ACCELERATION TYPES SUPPORTED BY OPTIMIZATION PROTOCOLS:

protocol    0   1   2   comment
V           +   +
VJ          +   +   +
FGH         +       +
VGJ         +   +   +   special protocol, not for widespread use
FJ          +       +   obsolete protocol, not recommended
FGJ         +       +   obsolete protocol, not recommended

DAFAULT VALUES:

protocol    0   1   2   comment
V               x       without acceleration it is so slooooooooow
VJ          x
FGH         x
VGJ                 x   we've implicitly turned (2) by passing gradient
FJ          x           obsolete protocol, not recommended
FGJ                 x   obsolete protocol, not recommended

NOTE: this  function should be called before optimization. Attempt to call
it during algorithm iterations may result in unexpected behavior.

NOTE: attempt to call this function with unsupported protocol/acceleration
combination will result in exception being thrown.

  -- ALGLIB --
     Copyright 14.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetacctype(minlmstate* state,
     ae_int_t acctype,
     ae_state *_state)
{


    ae_assert((acctype==0||acctype==1)||acctype==2, "MinLMSetAccType: incorrect AccType!", _state);
    if( acctype==0 )
    {
        state->maxmodelage = 0;
        state->makeadditers = ae_false;
        return;
    }
    if( acctype==1 )
    {
        ae_assert(state->hasfi, "MinLMSetAccType: AccType=1 is incompatible with current protocol!", _state);
        if( state->algomode==0 )
        {
            state->maxmodelage = 2*state->n;
        }
        else
        {
            state->maxmodelage = minlm_smallmodelage;
        }
        state->makeadditers = ae_false;
        return;
    }
    if( acctype==2 )
    {
        ae_assert(state->algomode==1||state->algomode==2, "MinLMSetAccType: AccType=2 is incompatible with current protocol!", _state);
        state->maxmodelage = 0;
        state->makeadditers = ae_true;
        return;
    }
}


/*************************************************************************
NOTES:

1. Depending on function used to create state  structure,  this  algorithm
   may accept Jacobian and/or Hessian and/or gradient.  According  to  the
   said above, there ase several versions of this function,  which  accept
   different sets of callbacks.

   This flexibility opens way to subtle errors - you may create state with
   MinLMCreateFGH() (optimization using Hessian), but call function  which
   does not accept Hessian. So when algorithm will request Hessian,  there
   will be no callback to call. In this case exception will be thrown.

   Be careful to avoid such errors because there is no way to find them at
   compile time - you can see them at runtime only.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool minlmiteration(minlmstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_bool bflag;
    ae_int_t iflag;
    double v;
    double s;
    double t;
    ae_int_t i;
    ae_int_t k;
    ae_bool result;


    
    /*
     * Reverse communication preparations
     * I know it looks ugly, but it works the same way
     * anywhere from C++ to Python.
     *
     * This code initializes locals by:
     * * random values determined during code
     *   generation - on first subroutine call
     * * values from previous call - on subsequent calls
     */
    if( state->rstate.stage>=0 )
    {
        n = state->rstate.ia.ptr.p_int[0];
        m = state->rstate.ia.ptr.p_int[1];
        iflag = state->rstate.ia.ptr.p_int[2];
        i = state->rstate.ia.ptr.p_int[3];
        k = state->rstate.ia.ptr.p_int[4];
        bflag = state->rstate.ba.ptr.p_bool[0];
        v = state->rstate.ra.ptr.p_double[0];
        s = state->rstate.ra.ptr.p_double[1];
        t = state->rstate.ra.ptr.p_double[2];
    }
    else
    {
        n = -983;
        m = -989;
        iflag = -834;
        i = 900;
        k = -287;
        bflag = ae_false;
        v = 214;
        s = -338;
        t = -686;
    }
    if( state->rstate.stage==0 )
    {
        goto lbl_0;
    }
    if( state->rstate.stage==1 )
    {
        goto lbl_1;
    }
    if( state->rstate.stage==2 )
    {
        goto lbl_2;
    }
    if( state->rstate.stage==3 )
    {
        goto lbl_3;
    }
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    if( state->rstate.stage==15 )
    {
        goto lbl_15;
    }
    if( state->rstate.stage==16 )
    {
        goto lbl_16;
    }
    if( state->rstate.stage==17 )
    {
        goto lbl_17;
    }
    
    /*
     * Routine body
     */
    
    /*
     * prepare
     */
    n = state->n;
    m = state->m;
    state->repiterationscount = 0;
    state->repterminationtype = 0;
    state->repnfunc = 0;
    state->repnjac = 0;
    state->repngrad = 0;
    state->repnhess = 0;
    state->repncholesky = 0;
    
    /*
     * Initial report of current point
     *
     * Note 1: we rewrite State.X twice because
     * user may accidentally change it after first call.
     *
     * Note 2: we set NeedF or NeedFI depending on what
     * information about function we have.
     */
    if( !state->xrep )
    {
        goto lbl_18;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    if( !state->hasf )
    {
        goto lbl_20;
    }
    state->needf = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needf = ae_false;
    goto lbl_21;
lbl_20:
    ae_assert(state->hasfi, "MinLM: internal error 2!", _state);
    state->needfi = ae_true;
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->needfi = ae_false;
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->f = v;
lbl_21:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->xupdated = ae_false;
lbl_18:
    
    /*
     * Prepare control variables
     */
    state->nu = 1;
    state->lambdav = -ae_maxrealnumber;
    state->modelage = state->maxmodelage+1;
    state->deltaxready = ae_false;
    state->deltafready = ae_false;
    
    /*
     * Main cycle.
     *
     * We move through it until either:
     * * one of the stopping conditions is met
     * * we decide that stopping conditions are too stringent
     *   and break from cycle
     *
     */
lbl_22:
    if( ae_false )
    {
        goto lbl_23;
    }
    
    /*
     * First, we have to prepare quadratic model for our function.
     * We use BFlag to ensure that model is prepared;
     * if it is false at the end of this block, something went wrong.
     *
     * We may either calculate brand new model or update old one.
     *
     * Before this block we have:
     * * State.XBase            - current position.
     * * State.DeltaX           - if DeltaXReady is True
     * * State.DeltaF           - if DeltaFReady is True
     *
     * After this block is over, we will have:
     * * State.XBase            - base point (unchanged)
     * * State.FBase            - F(XBase)
     * * State.GBase            - linear term
     * * State.QuadraticModel   - quadratic term
     * * State.LambdaV          - current estimate for lambda
     *
     * We also clear DeltaXReady/DeltaFReady flags
     * after initialization is done.
     */
    bflag = ae_false;
    if( !(state->algomode==0||state->algomode==1) )
    {
        goto lbl_24;
    }
    
    /*
     * Calculate f[] and Jacobian
     */
    if( !(state->modelage>state->maxmodelage||!(state->deltaxready&&state->deltafready)) )
    {
        goto lbl_26;
    }
    
    /*
     * Refresh model (using either finite differences or analytic Jacobian)
     */
    if( state->algomode!=0 )
    {
        goto lbl_28;
    }
    
    /*
     * Optimization using F values only.
     * Use finite differences to estimate Jacobian.
     */
    ae_assert(state->hasfi, "MinLMIteration: internal error when estimating Jacobian (no f[])", _state);
    k = 0;
lbl_30:
    if( k>n-1 )
    {
        goto lbl_32;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->x.ptr.p_double[k] = state->x.ptr.p_double[k]-state->diffstep;
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->fm1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->x.ptr.p_double[k] = state->x.ptr.p_double[k]+state->diffstep;
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->fp1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    v = 1/(2*state->diffstep);
    ae_v_moved(&state->j.ptr.pp_double[0][k], state->j.stride, &state->fp1.ptr.p_double[0], 1, ae_v_len(0,m-1), v);
    ae_v_subd(&state->j.ptr.pp_double[0][k], state->j.stride, &state->fm1.ptr.p_double[0], 1, ae_v_len(0,m-1), v);
    k = k+1;
    goto lbl_30;
lbl_32:
    
    /*
     * Calculate F(XBase)
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->needfi = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repnjac = state->repnjac+1;
    
    /*
     * New model
     */
    state->modelage = 0;
    goto lbl_29;
lbl_28:
    
    /*
     * Obtain f[] and Jacobian
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfij = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->needfij = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repnjac = state->repnjac+1;
    
    /*
     * New model
     */
    state->modelage = 0;
lbl_29:
    goto lbl_27;
lbl_26:
    
    /*
     * State.J contains Jacobian or its current approximation;
     * refresh it using secant updates:
     *
     * f(x0+dx) = f(x0) + J*dx,
     * J_new = J_old + u*h'
     * h = x_new-x_old
     * u = (f_new - f_old - J_old*h)/(h'h)
     *
     * We can explicitly generate h and u, but it is
     * preferential to do in-place calculations. Only
     * I-th row of J_old is needed to calculate u[I],
     * so we can update J row by row in one pass.
     *
     * NOTE: we expect that State.XBase contains new point,
     * State.FBase contains old point, State.DeltaX and
     * State.DeltaY contain updates from last step.
     */
    ae_assert(state->deltaxready&&state->deltafready, "MinLMIteration: uninitialized DeltaX/DeltaF", _state);
    t = ae_v_dotproduct(&state->deltax.ptr.p_double[0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_assert(ae_fp_neq(t,0), "MinLM: internal error (T=0)", _state);
    for(i=0; i<=m-1; i++)
    {
        v = ae_v_dotproduct(&state->j.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = (state->deltaf.ptr.p_double[i]-v)/t;
        ae_v_addd(&state->j.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    }
    ae_v_move(&state->fi.ptr.p_double[0], 1, &state->fibase.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_add(&state->fi.ptr.p_double[0], 1, &state->deltaf.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * Increase model age
     */
    state->modelage = state->modelage+1;
lbl_27:
    
    /*
     * Generate quadratic model:
     *     f(xbase+dx) =
     *       = (f0 + J*dx)'(f0 + J*dx)
     *       = f0^2 + dx'J'f0 + f0*J*dx + dx'J'J*dx
     *       = f0^2 + 2*f0*J*dx + dx'J'J*dx
     *
     * Note that we calculate 2*(J'J) instead of J'J because
     * our quadratic model is based on Tailor decomposition,
     * i.e. it has 0.5 before quadratic term.
     */
    rmatrixgemm(n, n, m, 2.0, &state->j, 0, 0, 1, &state->j, 0, 0, 0, 0.0, &state->quadraticmodel, 0, 0, _state);
    rmatrixmv(n, m, &state->j, 0, 0, 1, &state->fi, 0, &state->gbase, 0, _state);
    ae_v_muld(&state->gbase.ptr.p_double[0], 1, ae_v_len(0,n-1), 2);
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->fbase = v;
    ae_v_move(&state->fibase.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * set control variables
     */
    bflag = ae_true;
lbl_24:
    if( state->algomode!=2 )
    {
        goto lbl_33;
    }
    ae_assert(!state->hasfi, "MinLMIteration: internal error (HasFI is True in Hessian-based mode)", _state);
    
    /*
     * Obtain F, G, H
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfgh = ae_true;
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->needfgh = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repngrad = state->repngrad+1;
    state->repnhess = state->repnhess+1;
    rmatrixcopy(n, n, &state->h, 0, 0, &state->quadraticmodel, 0, 0, _state);
    ae_v_move(&state->gbase.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->fbase = state->f;
    
    /*
     * set control variables
     */
    bflag = ae_true;
    state->modelage = 0;
lbl_33:
    ae_assert(bflag, "MinLM: internal integrity check failed!", _state);
    state->deltaxready = ae_false;
    state->deltafready = ae_false;
    
    /*
     * If Lambda is not initialized, initialize it using quadratic model
     */
    if( ae_fp_less(state->lambdav,0) )
    {
        state->lambdav = 0;
        for(i=0; i<=n-1; i++)
        {
            state->lambdav = ae_maxreal(state->lambdav, ae_fabs(state->quadraticmodel.ptr.pp_double[i][i], _state), _state);
        }
        state->lambdav = 0.001*state->lambdav;
        if( ae_fp_eq(state->lambdav,0) )
        {
            state->lambdav = 1;
        }
    }
    
    /*
     * Test stopping conditions for function gradient
     */
    v = ae_v_dotproduct(&state->gbase.ptr.p_double[0], 1, &state->gbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    v = ae_sqrt(v, _state);
    if( ae_fp_greater(v,state->epsg) )
    {
        goto lbl_35;
    }
    if( state->modelage!=0 )
    {
        goto lbl_37;
    }
    
    /*
     * Model is fresh, we can rely on it and terminate algorithm
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_39;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->xupdated = ae_false;
lbl_39:
    result = ae_false;
    return result;
    goto lbl_38;
lbl_37:
    
    /*
     * Model is not fresh, we should refresh it and test
     * conditions once more
     */
    state->modelage = state->maxmodelage+1;
    goto lbl_22;
lbl_38:
lbl_35:
    
    /*
     * Find value of Levenberg-Marquardt damping parameter which:
     * * leads to positive definite damped model
     * * within bounds specified by StpMax
     * * generates step which decreases function value
     *
     * After this block IFlag is set to:
     * * -2, if model update is needed (either Lambda growth is too large
     *       or step is too short, but we can't rely on model and stop iterations)
     * * -1, if model is fresh, Lambda have grown too large, termination is needed
     * *  0, if everything is OK, continue iterations
     *
     * State.Nu can have any value on enter, but after exit it is set to 1.0
     */
    iflag = -99;
lbl_41:
    if( ae_false )
    {
        goto lbl_42;
    }
    
    /*
     * Do we need model update?
     */
    if( state->modelage>0&&ae_fp_greater_eq(state->nu,minlm_suspiciousnu) )
    {
        iflag = -2;
        goto lbl_42;
    }
    
    /*
     * DampedModel = QuadraticModel+lambda*I
     */
    rmatrixcopy(n, n, &state->quadraticmodel, 0, 0, &state->dampedmodel, 0, 0, _state);
    for(i=0; i<=n-1; i++)
    {
        state->dampedmodel.ptr.pp_double[i][i] = state->dampedmodel.ptr.pp_double[i][i]+state->lambdav;
    }
    
    /*
     * 1. try to solve (RawModel+Lambda*I)*dx = -g.
     *    increase lambda if left part is not positive definite.
     * 2. bound step by StpMax
     *    increase lambda if step is larger than StpMax
     *
     * We use BFlag variable to indicate that we have to increase Lambda.
     * If it is False, we will try to increase Lambda and move to new iteration.
     */
    bflag = ae_true;
    state->repncholesky = state->repncholesky+1;
    if( spdmatrixcholeskyrec(&state->dampedmodel, 0, n, ae_true, &state->choleskybuf, _state) )
    {
        ae_v_moveneg(&state->xdir.ptr.p_double[0], 1, &state->gbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
        fblscholeskysolve(&state->dampedmodel, 1.0, n, ae_true, &state->xdir, &state->choleskybuf, _state);
        v = ae_v_dotproduct(&state->xdir.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
        if( ae_isfinite(v, _state) )
        {
            v = ae_sqrt(v, _state);
            if( ae_fp_greater(state->stpmax,0)&&ae_fp_greater(v,state->stpmax) )
            {
                bflag = ae_false;
            }
        }
        else
        {
            bflag = ae_false;
        }
    }
    else
    {
        bflag = ae_false;
    }
    if( !bflag )
    {
        
        /*
         * Solution failed:
         * try to increase lambda to make matrix positive definite and continue.
         */
        if( !minlm_increaselambda(&state->lambdav, &state->nu, _state) )
        {
            iflag = -1;
            goto lbl_42;
        }
        goto lbl_41;
    }
    
    /*
     * Step in State.XDir and it is bounded by StpMax.
     *
     * We should check stopping conditions on step size here.
     * DeltaX, which is used for secant updates, is initialized here.
     *
     * This code is a bit tricky because sometimes XDir<>0, but
     * it is so small that XDir+XBase==XBase (in finite precision
     * arithmetics). So we set DeltaX to XBase, then
     * add XDir, and then subtract XBase to get exact value of
     * DeltaX.
     *
     * Step length is estimated using DeltaX.
     *
     * NOTE: stopping conditions are tested
     * for fresh models only (ModelAge=0)
     */
    ae_v_move(&state->deltax.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_add(&state->deltax.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_sub(&state->deltax.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->deltaxready = ae_true;
    v = ae_v_dotproduct(&state->deltax.ptr.p_double[0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
    v = ae_sqrt(v, _state);
    if( ae_fp_greater(v,state->epsx) )
    {
        goto lbl_43;
    }
    if( state->modelage!=0 )
    {
        goto lbl_45;
    }
    
    /*
     * Step is too short, model is fresh and we can rely on it.
     * Terminating.
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_47;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->xupdated = ae_false;
lbl_47:
    result = ae_false;
    return result;
    goto lbl_46;
lbl_45:
    
    /*
     * Step is suspiciously short, but model is not fresh
     * and we can't rely on it.
     */
    iflag = -2;
    goto lbl_42;
lbl_46:
lbl_43:
    
    /*
     * Let's evaluate new step:
     * a) if we have Fi vector, we evaluate it using rcomm, and
     *    then we manually calculate State.F as sum of squares of Fi[]
     * b) if we have F value, we just evaluate it through rcomm interface
     *
     * We prefer (a) because we may need Fi vector for additional
     * iterations
     */
    ae_assert(state->hasfi||state->hasf, "MinLM: internal error 2!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_add(&state->x.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    if( !state->hasfi )
    {
        goto lbl_49;
    }
    state->needfi = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->needfi = ae_false;
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->f = v;
    ae_v_move(&state->deltaf.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_sub(&state->deltaf.ptr.p_double[0], 1, &state->fibase.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->deltafready = ae_true;
    goto lbl_50;
lbl_49:
    state->needf = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->needf = ae_false;
lbl_50:
    state->repnfunc = state->repnfunc+1;
    if( ae_fp_greater_eq(state->f,state->fbase) )
    {
        
        /*
         * Increase lambda and continue
         */
        if( !minlm_increaselambda(&state->lambdav, &state->nu, _state) )
        {
            iflag = -1;
            goto lbl_42;
        }
        goto lbl_41;
    }
    
    /*
     * We've found our step!
     */
    iflag = 0;
    goto lbl_42;
    goto lbl_41;
lbl_42:
    state->nu = 1;
    ae_assert(iflag>=-2&&iflag<=0, "MinLM: internal integrity check failed!", _state);
    if( iflag==-2 )
    {
        state->modelage = state->maxmodelage+1;
        goto lbl_22;
    }
    if( iflag==-1 )
    {
        goto lbl_23;
    }
    
    /*
     * Levenberg-Marquardt step is ready.
     * Compare predicted vs. actual decrease and decide what to do with lambda.
     *
     * NOTE: we expect that State.DeltaX contains direction of step,
     * State.F contains function value at new point.
     */
    ae_assert(state->deltaxready, "MinLM: deltaX is not ready", _state);
    t = 0;
    for(i=0; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&state->quadraticmodel.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
        t = t+state->deltax.ptr.p_double[i]*state->gbase.ptr.p_double[i]+0.5*state->deltax.ptr.p_double[i]*v;
    }
    state->predicteddecrease = -t;
    state->actualdecrease = -(state->f-state->fbase);
    if( ae_fp_less_eq(state->predicteddecrease,0) )
    {
        goto lbl_23;
    }
    v = state->actualdecrease/state->predicteddecrease;
    if( ae_fp_greater_eq(v,0.1) )
    {
        goto lbl_51;
    }
    if( minlm_increaselambda(&state->lambdav, &state->nu, _state) )
    {
        goto lbl_53;
    }
    
    /*
     * Lambda is too large, we have to break iterations.
     */
    state->repterminationtype = 7;
    if( !state->xrep )
    {
        goto lbl_55;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->xupdated = ae_false;
lbl_55:
    result = ae_false;
    return result;
lbl_53:
lbl_51:
    if( ae_fp_greater(v,0.5) )
    {
        minlm_decreaselambda(&state->lambdav, &state->nu, _state);
    }
    
    /*
     * Accept step, report it and
     * test stopping conditions on iterations count and function decrease.
     *
     * NOTE: we expect that State.DeltaX contains direction of step,
     * State.F contains function value at new point.
     *
     * NOTE2: we should update XBase ONLY. In the beginning of the next
     * iteration we expect that State.FIBase is NOT updated and
     * contains old value of a function vector.
     */
    ae_v_add(&state->xbase.ptr.p_double[0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( !state->xrep )
    {
        goto lbl_57;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->xupdated = ae_false;
lbl_57:
    state->repiterationscount = state->repiterationscount+1;
    if( state->repiterationscount>=state->maxits&&state->maxits>0 )
    {
        state->repterminationtype = 5;
    }
    if( state->modelage==0 )
    {
        if( ae_fp_less_eq(ae_fabs(state->f-state->fbase, _state),state->epsf*ae_maxreal(1, ae_maxreal(ae_fabs(state->f, _state), ae_fabs(state->fbase, _state), _state), _state)) )
        {
            state->repterminationtype = 1;
        }
    }
    if( state->repterminationtype<=0 )
    {
        goto lbl_59;
    }
    if( !state->xrep )
    {
        goto lbl_61;
    }
    
    /*
     * Report: XBase contains new point, F contains function value at new point
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->xupdated = ae_false;
lbl_61:
    result = ae_false;
    return result;
lbl_59:
    state->modelage = state->modelage+1;
    
    /*
     * Additional iterations for unconstrained problems:
     * preconditioned L-BFGS is used.
     *
     * NOTE: additional iterations are incompatible with secant updates
     * because they invalidate
     */
    if( !(ae_fp_eq(state->stpmax,0)&&state->makeadditers) )
    {
        goto lbl_63;
    }
    ae_assert(state->hasg||state->m!=0, "MinLM: no grad or Jacobian for additional iterations", _state);
    
    /*
     * Make preconditioned iterations
     */
    minlbfgssetcholeskypreconditioner(&state->internalstate, &state->dampedmodel, ae_true, _state);
    minlbfgsrestartfrom(&state->internalstate, &state->xbase, _state);
lbl_65:
    if( !minlbfgsiteration(&state->internalstate, _state) )
    {
        goto lbl_66;
    }
    if( !state->internalstate.needfg )
    {
        goto lbl_67;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->internalstate.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    if( !state->hasg )
    {
        goto lbl_69;
    }
    state->needfg = ae_true;
    state->rstate.stage = 15;
    goto lbl_rcomm;
lbl_15:
    state->needfg = ae_false;
    state->repngrad = state->repngrad+1;
    ae_v_move(&state->internalstate.g.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->internalstate.f = state->f;
    goto lbl_70;
lbl_69:
    state->needfij = ae_true;
    state->rstate.stage = 16;
    goto lbl_rcomm;
lbl_16:
    state->needfij = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repnjac = state->repnjac+1;
    for(i=0; i<=n-1; i++)
    {
        state->internalstate.g.ptr.p_double[i] = 0;
    }
    for(i=0; i<=m-1; i++)
    {
        v = 2*state->fi.ptr.p_double[i];
        ae_v_addd(&state->internalstate.g.ptr.p_double[0], 1, &state->j.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        state->internalstate.f = state->internalstate.f+ae_sqr(state->fi.ptr.p_double[i], _state);
    }
lbl_70:
lbl_67:
    goto lbl_65;
lbl_66:
    minlbfgsresultsbuf(&state->internalstate, &state->xbase, &state->internalrep, _state);
    
    /*
     * Invalidate DeltaX/DeltaF (control variables used for integrity checks)
     */
    state->deltaxready = ae_false;
    state->deltafready = ae_false;
lbl_63:
    goto lbl_22;
lbl_23:
    
    /*
     * Lambda is too large, we have to break iterations.
     */
    state->repterminationtype = 7;
    if( !state->xrep )
    {
        goto lbl_71;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 17;
    goto lbl_rcomm;
lbl_17:
    state->xupdated = ae_false;
lbl_71:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = iflag;
    state->rstate.ia.ptr.p_int[3] = i;
    state->rstate.ia.ptr.p_int[4] = k;
    state->rstate.ba.ptr.p_bool[0] = bflag;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = s;
    state->rstate.ra.ptr.p_double[2] = t;
    return result;
}


/*************************************************************************
Levenberg-Marquardt algorithm results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report;
                see comments for this structure for more info.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresults(minlmstate* state,
     /* Real    */ ae_vector* x,
     minlmreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minlmreport_clear(rep);

    minlmresultsbuf(state, x, rep, _state);
}


/*************************************************************************
Levenberg-Marquardt algorithm results

Buffered implementation of MinLMResults(), which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresultsbuf(minlmstate* state,
     /* Real    */ ae_vector* x,
     minlmreport* rep,
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->terminationtype = state->repterminationtype;
    rep->nfunc = state->repnfunc;
    rep->njac = state->repnjac;
    rep->ngrad = state->repngrad;
    rep->nhess = state->repnhess;
    rep->ncholesky = state->repncholesky;
}


/*************************************************************************
This  subroutine  restarts  LM  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used for reverse communication previously
                allocated with MinLMCreateXXX call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minlmrestartfrom(minlmstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinLMRestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinLMRestartFrom: X contains infinite or NaN values!", _state);
    ae_v_move(&state->xbase.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ba, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 2+1, _state);
    state->rstate.stage = -1;
    minlm_clearrequestfields(state, _state);
}


/*************************************************************************
Prepare internal structures (except for RComm).

Note: M must be zero for FGH mode, non-zero for V/VJ/FJ/FGJ mode.
*************************************************************************/
static void minlm_lmprepare(ae_int_t n,
     ae_int_t m,
     ae_bool havegrad,
     minlmstate* state,
     ae_state *_state)
{
    ae_int_t i;


    if( n<=0||m<0 )
    {
        return;
    }
    if( havegrad )
    {
        ae_vector_set_length(&state->g, n, _state);
    }
    if( m!=0 )
    {
        ae_matrix_set_length(&state->j, m, n, _state);
        ae_vector_set_length(&state->fi, m, _state);
        ae_vector_set_length(&state->fibase, m, _state);
        ae_vector_set_length(&state->deltaf, m, _state);
        ae_vector_set_length(&state->fm2, m, _state);
        ae_vector_set_length(&state->fm1, m, _state);
        ae_vector_set_length(&state->fp2, m, _state);
        ae_vector_set_length(&state->fp1, m, _state);
    }
    else
    {
        ae_matrix_set_length(&state->h, n, n, _state);
    }
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->deltax, n, _state);
    ae_matrix_set_length(&state->quadraticmodel, n, n, _state);
    ae_matrix_set_length(&state->dampedmodel, n, n, _state);
    ae_vector_set_length(&state->xbase, n, _state);
    ae_vector_set_length(&state->gbase, n, _state);
    ae_vector_set_length(&state->xdir, n, _state);
    
    /*
     * prepare internal L-BFGS
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = 0;
    }
    minlbfgscreate(n, ae_minint(minlm_additers, n, _state), &state->x, &state->internalstate, _state);
    minlbfgssetcond(&state->internalstate, 0.0, 0.0, 0.0, ae_minint(minlm_additers, n, _state), _state);
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void minlm_clearrequestfields(minlmstate* state, ae_state *_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->needfgh = ae_false;
    state->needfij = ae_false;
    state->needfi = ae_false;
    state->xupdated = ae_false;
}


/*************************************************************************
Increases lambda, returns False when there is a danger of overflow
*************************************************************************/
static ae_bool minlm_increaselambda(double* lambdav,
     double* nu,
     ae_state *_state)
{
    double lnlambda;
    double lnnu;
    double lnlambdaup;
    double lnmax;
    ae_bool result;


    result = ae_false;
    lnlambda = ae_log(*lambdav, _state);
    lnlambdaup = ae_log(minlm_lambdaup, _state);
    lnnu = ae_log(*nu, _state);
    lnmax = ae_log(ae_maxrealnumber, _state);
    if( ae_fp_greater(lnlambda+lnlambdaup+lnnu,0.25*lnmax) )
    {
        return result;
    }
    if( ae_fp_greater(lnnu+ae_log(2, _state),lnmax) )
    {
        return result;
    }
    *lambdav = *lambdav*minlm_lambdaup*(*nu);
    *nu = *nu*2;
    result = ae_true;
    return result;
}


/*************************************************************************
Decreases lambda, but leaves it unchanged when there is danger of underflow.
*************************************************************************/
static void minlm_decreaselambda(double* lambdav,
     double* nu,
     ae_state *_state)
{


    *nu = 1;
    if( ae_fp_less(ae_log(*lambdav, _state)+ae_log(minlm_lambdadown, _state),ae_log(ae_minrealnumber, _state)) )
    {
        *lambdav = ae_minrealnumber;
    }
    else
    {
        *lambdav = *lambdav*minlm_lambdadown;
    }
}


ae_bool _minlmstate_init(minlmstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fi, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->j, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->h, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xbase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fibase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gbase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->quadraticmodel, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->dampedmodel, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xdir, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->deltax, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->deltaf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->choleskybuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fm2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fm1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fp2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsstate_init(&p->internalstate, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsreport_init(&p->internalrep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minlmstate_init_copy(minlmstate* dst, minlmstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->m = src->m;
    dst->diffstep = src->diffstep;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->maxmodelage = src->maxmodelage;
    dst->makeadditers = src->makeadditers;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->fi, &src->fi, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->j, &src->j, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->h, &src->h, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->needfgh = src->needfgh;
    dst->needfij = src->needfij;
    dst->needfi = src->needfi;
    dst->xupdated = src->xupdated;
    dst->algomode = src->algomode;
    dst->hasf = src->hasf;
    dst->hasfi = src->hasfi;
    dst->hasg = src->hasg;
    if( !ae_vector_init_copy(&dst->xbase, &src->xbase, _state, make_automatic) )
        return ae_false;
    dst->fbase = src->fbase;
    if( !ae_vector_init_copy(&dst->fibase, &src->fibase, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gbase, &src->gbase, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->quadraticmodel, &src->quadraticmodel, _state, make_automatic) )
        return ae_false;
    dst->lambdav = src->lambdav;
    dst->nu = src->nu;
    if( !ae_matrix_init_copy(&dst->dampedmodel, &src->dampedmodel, _state, make_automatic) )
        return ae_false;
    dst->modelage = src->modelage;
    if( !ae_vector_init_copy(&dst->xdir, &src->xdir, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->deltax, &src->deltax, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->deltaf, &src->deltaf, _state, make_automatic) )
        return ae_false;
    dst->deltaxready = src->deltaxready;
    dst->deltafready = src->deltafready;
    dst->repiterationscount = src->repiterationscount;
    dst->repterminationtype = src->repterminationtype;
    dst->repnfunc = src->repnfunc;
    dst->repnjac = src->repnjac;
    dst->repngrad = src->repngrad;
    dst->repnhess = src->repnhess;
    dst->repncholesky = src->repncholesky;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->choleskybuf, &src->choleskybuf, _state, make_automatic) )
        return ae_false;
    dst->actualdecrease = src->actualdecrease;
    dst->predicteddecrease = src->predicteddecrease;
    if( !ae_vector_init_copy(&dst->fm2, &src->fm2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->fm1, &src->fm1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->fp2, &src->fp2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->fp1, &src->fp1, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsstate_init_copy(&dst->internalstate, &src->internalstate, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsreport_init_copy(&dst->internalrep, &src->internalrep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _minlmstate_clear(minlmstate* p)
{
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->fi);
    ae_matrix_clear(&p->j);
    ae_matrix_clear(&p->h);
    ae_vector_clear(&p->g);
    ae_vector_clear(&p->xbase);
    ae_vector_clear(&p->fibase);
    ae_vector_clear(&p->gbase);
    ae_matrix_clear(&p->quadraticmodel);
    ae_matrix_clear(&p->dampedmodel);
    ae_vector_clear(&p->xdir);
    ae_vector_clear(&p->deltax);
    ae_vector_clear(&p->deltaf);
    _rcommstate_clear(&p->rstate);
    ae_vector_clear(&p->choleskybuf);
    ae_vector_clear(&p->fm2);
    ae_vector_clear(&p->fm1);
    ae_vector_clear(&p->fp2);
    ae_vector_clear(&p->fp1);
    _minlbfgsstate_clear(&p->internalstate);
    _minlbfgsreport_clear(&p->internalrep);
}


ae_bool _minlmreport_init(minlmreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _minlmreport_init_copy(minlmreport* dst, minlmreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->terminationtype = src->terminationtype;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
    return ae_true;
}


void _minlmreport_clear(minlmreport* p)
{
}


/*$ End $*/
