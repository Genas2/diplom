/*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _minbleic_h
#define _minbleic_h

#include "aenv.h"
#include "ialglib.h"
#include "linmin.h"
#include "apserv.h"
#include "mincg.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"


/*$ Declarations $*/


/*************************************************************************
This object stores nonlinear optimizer state.
You should use functions provided by MinBLEIC subpackage to work with this
object
*************************************************************************/
typedef struct
{
    ae_int_t n;
    double innerepsg;
    double innerepsf;
    double innerepsx;
    double outerepsx;
    double outerepsi;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_int_t cgtype;
    double mustart;
    double mudecay;
    ae_vector x;
    double f;
    ae_vector g;
    ae_bool needfg;
    ae_bool xupdated;
    rcommstate rstate;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repterminationtype;
    double repdebugeqerr;
    double repdebugfs;
    double repdebugff;
    double repdebugdx;
    ae_vector xcur;
    ae_vector xprev;
    ae_vector xstart;
    ae_int_t itsleft;
    ae_int_t mucounter;
    ae_matrix ce;
    ae_matrix ci;
    ae_matrix cebasis;
    ae_matrix cesvl;
    ae_vector xe;
    ae_int_t cecnt;
    ae_int_t cicnt;
    ae_int_t cedim;
    ae_vector bndl;
    ae_vector hasbndl;
    ae_vector bndu;
    ae_vector hasbndu;
    ae_vector lm;
    ae_int_t lmcnt;
    double mu;
    ae_vector w;
    ae_vector tmp0;
    ae_vector tmp1;
    ae_vector r;
    double v0;
    double v1;
    double v2;
    double t;
    double errfeas;
    double errslack;
    double gnorm;
    double mpgnorm;
    double lmdif;
    double lmnorm;
    double lmgrowth;
    double mba;
    double boundary;
    ae_bool closetobarrier;
    double bndmax;
    linminstate lstate;
    mincgstate cgstate;
    mincgreport cgrep;
} minbleicstate;


/*************************************************************************
This structure stores optimization report:
* InnerIterationsCount      number of inner iterations
* OuterIterationsCount      number of outer iterations
* NFEV                      number of gradient evaluations

There are additional fields which can be used for debugging:
* DebugEqErr                error in the equality constraints (2-norm)
* DebugFS                   f, calculated at projection of initial point
                            to the feasible set
* DebugFF                   f, calculated at the final point
* DebugDX                   |X_start-X_final|
*************************************************************************/
typedef struct
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
    double debugeqerr;
    double debugfs;
    double debugff;
    double debugdx;
} minbleicreport;


/*$ Body $*/


/*************************************************************************
                     BOUND CONSTRAINED OPTIMIZATION
       WITH ADDITIONAL LINEAR EQUALITY AND INEQUALITY CONSTRAINTS

DESCRIPTION:
The  subroutine  minimizes  function   F(x)  of N arguments subject to any
combination of:
* bound constraints
* linear inequality constraints
* linear equality constraints

REQUIREMENTS:
* function value and gradient
* grad(f) must be Lipschitz continuous on a level set: L = { x : f(x)<=f(x0) }
* function must be defined even in the infeasible points (algorithm make take
  steps in the infeasible area before converging to the feasible point)
* starting point X0 must be feasible or not too far away from the feasible set
* problem must satisfy strict complementary conditions

USAGE:

Constrained optimization if far more complex than the unconstrained one.
Here we give very brief outline of the BLEIC optimizer. We strongly recommend
you to read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on optimization, which is available at http://www.alglib.net/optimization/

1. User initializes algorithm state with MinBLEICCreate() call

2. USer adds boundary and/or linear constraints by calling
   MinBLEICSetBC() and MinBLEICSetLC() functions.

3. User sets stopping conditions for underlying unconstrained solver
   with MinBLEICSetInnerCond() call.
   This function controls accuracy of underlying optimization algorithm.

4. User sets stopping conditions for outer iteration by calling
   MinBLEICSetOuterCond() function.
   This function controls handling of boundary and inequality constraints.

5. User tunes barrier parameters:
   * barrier width with MinBLEICSetBarrierWidth() call
   * (optionally) dynamics of the barrier width with MinBLEICSetBarrierDecay() call
   These functions control handling of boundary and inequality constraints.

6. Additionally, user may set limit on number of internal iterations
   by MinBLEICSetMaxIts() call.
   This function allows to prevent algorithm from looping forever.

7. User calls MinBLEICOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

8. User calls MinBLEICResults() to get solution

9. Optionally user may call MinBLEICRestartFrom() to solve another problem
   with same N but another starting point.
   MinBLEICRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size ofX
    X       -   starting point, array[N]:
                * it is better to set X to a feasible point
                * but X can be infeasible, in which case algorithm will try
                  to find feasible point first, using X as initial
                  approximation.

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleiccreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     minbleicstate* state,
     ae_state *_state);


/*************************************************************************
This function sets boundary constraints for BLEIC optimizer.

Boundary constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF.
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbc(minbleicstate* state,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state);


/*************************************************************************
This function sets linear constraints for BLEIC optimizer.

Linear constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetlc(minbleicstate* state,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state);


/*************************************************************************
This function sets stopping conditions for the underlying nonlinear CG
optimizer. It controls overall accuracy of solution. These conditions
should be strict enough in order for algorithm to converge.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                Algorithm finishes its work if 2-norm of the Lagrangian
                gradient is less than or equal to EpsG.
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |X(k+1)-X(k)| <= EpsX is fulfilled.

Passing EpsG=0, EpsF=0 and EpsX=0 (simultaneously) will lead to
automatic stopping criterion selection.

These conditions are used to terminate inner iterations. However, you
need to tune termination conditions for outer iterations too.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetinnercond(minbleicstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_state *_state);


/*************************************************************************
This function sets stopping conditions for outer iteration of BLEIC algo.

These conditions control accuracy of constraint handling and amount of
infeasibility allowed in the solution.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsX    -   >0, stopping condition on outer iteration step length
    EpsI    -   >0, stopping condition on infeasibility
    
Both EpsX and EpsI must be non-zero.

MEANING OF EpsX

EpsX  is  a  stopping  condition for outer iterations. Algorithm will stop
when  solution  of  the  current  modified  subproblem will be within EpsX
(using 2-norm) of the previous solution.

MEANING OF EpsI

EpsI controls feasibility properties -  algorithm  won't  stop  until  all
inequality constraints will be satisfied with error (distance from current
point to the feasible area) at most EpsI.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetoutercond(minbleicstate* state,
     double epsx,
     double epsi,
     ae_state *_state);


/*************************************************************************
This function sets initial barrier width.

BLEIC optimizer uses  modified  barrier  functions  to  handle  inequality
constraints. These functions are almost constant in the inner parts of the
feasible  area,  but  grow rapidly to the infinity OUTSIDE of the feasible
area. Barrier width is a distance from feasible area to  the  point  where
modified barrier function becomes infinite.

Barrier width must be:
* small enough (below some problem-dependent value) in order for algorithm
  to  converge.  Necessary  condition  is that the target function must be
  well described by linear model in the areas as small as barrier width.
* not VERY small (in order to avoid  difficulties  associated  with  rapid
  changes in the modified function, ill-conditioning, round-off issues).

Choosing  appropriate  barrier  width  is  very  important  for  efficient
optimization, and it often requires error  and  trial.  You  can  use  two
strategies when choosing barrier width:
* set barrier width with MinBLEICSetBarrierWidth() call. In this case  you
  should try different barrier widths and examine results.
* set decreasing barrier width by combining  MinBLEICSetBarrierWidth() and
  MinBLEICSetBarrierDecay()  calls.  In  this case algorithm will decrease
  barrier  width  after  each  outer iteration until it encounters optimal
  barrier width.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    Mu      -   >0, initial barrier width

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierwidth(minbleicstate* state,
     double mu,
     ae_state *_state);


/*************************************************************************
This function sets decay coefficient for barrier width.

By default, no barrier decay is used (Decay=1.0).

BLEIC optimizer uses  modified  barrier  functions  to  handle  inequality
constraints. These functions are almost constant in the inner parts of the
feasible  area,  but  grow rapidly to the infinity OUTSIDE of the feasible
area. Barrier width is a distance from feasible area to  the  point  where
modified barrier function becomes infinite. Decay coefficient allows us to
decrease  barrier  width  from  the  initial  (suboptimial)  value   until
optimal value will be met.

We recommend you either to set MuDecay=1.0 (no decay) or use some moderate
value like 0.5-0.7

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    MuDecay -   0<MuDecay<=1, decay coefficient

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierdecay(minbleicstate* state,
     double mudecay,
     ae_state *_state);


/*************************************************************************
This function allows to stop algorithm after specified number of inner
iterations.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    MaxIts  -   maximum number of inner iterations.
                If MaxIts=0, the number of iterations is unlimited.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetmaxits(minbleicstate* state,
     ae_int_t maxits,
     ae_state *_state);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinBLEICOptimize().

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetxrep(minbleicstate* state,
     ae_bool needxrep,
     ae_state *_state);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  lead   to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetstpmax(minbleicstate* state,
     double stpmax,
     ae_state *_state);


/*************************************************************************

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool minbleiciteration(minbleicstate* state, ae_state *_state);


/*************************************************************************
BLEIC results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -3    inconsistent constraints. Feasible point is
                            either nonexistent or too hard to find. Try to
                            restart optimizer with better initial
                            approximation
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    *  4    conditions on constraints are fulfilled
                            with error less than or equal to EpsC
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible,
                            X contains best point found so far.
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresults(minbleicstate* state,
     /* Real    */ ae_vector* x,
     minbleicreport* rep,
     ae_state *_state);


/*************************************************************************
BLEIC results

Buffered implementation of MinBLEICResults() which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresultsbuf(minbleicstate* state,
     /* Real    */ ae_vector* x,
     minbleicreport* rep,
     ae_state *_state);


/*************************************************************************
This subroutine restarts algorithm from new point.
All optimization parameters (including constraints) are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have  same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicrestartfrom(minbleicstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);


/*************************************************************************
Modified barrier function calculated at point X with respect to the
barrier parameter Mu.

Functions, its first and second derivatives are calculated.
*************************************************************************/
void barrierfunc(double x,
     double mu,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state);
ae_bool _minbleicstate_init(minbleicstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _minbleicstate_init_copy(minbleicstate* dst, minbleicstate* src, ae_state *_state, ae_bool make_automatic);
void _minbleicstate_clear(minbleicstate* p);
ae_bool _minbleicreport_init(minbleicreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _minbleicreport_init_copy(minbleicreport* dst, minbleicreport* src, ae_state *_state, ae_bool make_automatic);
void _minbleicreport_clear(minbleicreport* p);


/*$ End $*/
#endif

