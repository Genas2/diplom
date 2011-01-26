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


#include <stdafx.h>
#include "minbleic.h"


/*$ Declarations $*/
static double minbleic_svdtol = 100;
static double minbleic_lmtol = 100;
static double minbleic_maxlmgrowth = 1000;
static double minbleic_minlagrangemul = 1.0E-50;
static double minbleic_maxlagrangemul = 1.0E+50;
static double minbleic_maxouterits = 20;
static ae_int_t minbleic_mucountdownlen = 15;
static ae_int_t minbleic_cscountdownlen = 5;
static void minbleic_clearrequestfields(minbleicstate* state,
     ae_state *_state);
static void minbleic_makeprojection(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* r,
     double* rnorm2,
     ae_state *_state);
static void minbleic_modifytargetfunction(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* r,
     double rnorm2,
     double* f,
     /* Real    */ ae_vector* g,
     double* gnorm,
     double* mpgnorm,
     double* mba,
     double* fierr,
     double* cserr,
     ae_state *_state);
static void minbleic_penaltyfunction(minbleicstate* state,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* r,
     double* mba,
     double* fierr,
     ae_state *_state);


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
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_matrix c;
    ae_vector ct;

    ae_frame_make(_state, &_frame_block);
    _minbleicstate_clear(state);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);

    ae_assert(n>=1, "MinBLEICCreate: N<1", _state);
    ae_assert(x->cnt>=n, "MinBLEICCreate: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinBLEICCreate: X contains infinite or NaN values!", _state);
    
    /*
     * Initialize.
     */
    state->n = n;
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->g, n, _state);
    ae_vector_set_length(&state->bndl, n, _state);
    ae_vector_set_length(&state->hasbndl, n, _state);
    ae_vector_set_length(&state->bndu, n, _state);
    ae_vector_set_length(&state->hasbndu, n, _state);
    ae_vector_set_length(&state->xcur, n, _state);
    ae_vector_set_length(&state->xprev, n, _state);
    ae_vector_set_length(&state->xstart, n, _state);
    ae_vector_set_length(&state->xe, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->bndl.ptr.p_double[i] = _state->v_neginf;
        state->hasbndl.ptr.p_bool[i] = ae_false;
        state->bndu.ptr.p_double[i] = _state->v_posinf;
        state->hasbndu.ptr.p_bool[i] = ae_false;
    }
    minbleicsetlc(state, &c, &ct, 0, _state);
    minbleicsetinnercond(state, 0.0, 0.0, 0.0, _state);
    minbleicsetoutercond(state, 1.0E-6, 1.0E-6, _state);
    minbleicsetbarrierwidth(state, 1.0E-3, _state);
    minbleicsetbarrierdecay(state, 1.0, _state);
    minbleicsetmaxits(state, 0, _state);
    minbleicsetxrep(state, ae_false, _state);
    minbleicsetstpmax(state, 0.0, _state);
    minbleicrestartfrom(state, x, _state);
    mincgcreate(n, x, &state->cgstate, _state);
    ae_frame_leave(_state);
}


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
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;


    n = state->n;
    ae_assert(bndl->cnt>=n, "MinBLEICSetBC: Length(BndL)<N", _state);
    ae_assert(bndu->cnt>=n, "MinBLEICSetBC: Length(BndU)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_isfinite(bndl->ptr.p_double[i], _state)||ae_isneginf(bndl->ptr.p_double[i], _state), "MinBLEICSetBC: BndL contains NAN or +INF", _state);
        ae_assert(ae_isfinite(bndu->ptr.p_double[i], _state)||ae_isposinf(bndu->ptr.p_double[i], _state), "MinBLEICSetBC: BndL contains NAN or -INF", _state);
        state->bndl.ptr.p_double[i] = bndl->ptr.p_double[i];
        state->hasbndl.ptr.p_bool[i] = ae_isfinite(bndl->ptr.p_double[i], _state);
        state->bndu.ptr.p_double[i] = bndu->ptr.p_double[i];
        state->hasbndu.ptr.p_bool[i] = ae_isfinite(bndu->ptr.p_double[i], _state);
    }
}


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
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    ae_int_t idx;
    ae_bool b;
    double v;


    n = state->n;
    
    /*
     * First, check for errors in the inputs
     */
    ae_assert(k>=0, "MinBLEICSetLC: K<0", _state);
    ae_assert(c->cols>=n+1||k==0, "MinBLEICSetLC: Cols(C)<N+1", _state);
    ae_assert(c->rows>=k, "MinBLEICSetLC: Rows(C)<K", _state);
    ae_assert(ct->cnt>=k, "MinBLEICSetLC: Length(CT)<K", _state);
    ae_assert(apservisfinitematrix(c, k, n+1, _state), "MinBLEICSetLC: C contains infinite or NaN values!", _state);
    
    /*
     * Determine number of constraints,
     * allocate space and copy
     */
    state->cecnt = 0;
    state->cicnt = 0;
    for(i=0; i<=k-1; i++)
    {
        if( ct->ptr.p_int[i]!=0 )
        {
            state->cicnt = state->cicnt+1;
        }
        else
        {
            state->cecnt = state->cecnt+1;
        }
    }
    rmatrixsetlengthatleast(&state->ci, state->cicnt, n+1, _state);
    rmatrixsetlengthatleast(&state->ce, state->cecnt, n+1, _state);
    idx = 0;
    for(i=0; i<=k-1; i++)
    {
        if( ct->ptr.p_int[i]!=0 )
        {
            ae_v_move(&state->ci.ptr.pp_double[idx][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
            if( ct->ptr.p_int[i]<0 )
            {
                ae_v_muld(&state->ci.ptr.pp_double[idx][0], 1, ae_v_len(0,n), -1);
            }
            idx = idx+1;
        }
    }
    idx = 0;
    for(i=0; i<=k-1; i++)
    {
        if( ct->ptr.p_int[i]==0 )
        {
            ae_v_move(&state->ce.ptr.pp_double[idx][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
            idx = idx+1;
        }
    }
    
    /*
     * Calculate ortohognal basis of row space of CE.
     * Determine actual basis size, drop vectors corresponding to
     * small singular values.
     *
     * NOTE: it is important to use "W[I]>W[0]*Tol" form (strict
     * inequality instead of non-strict) because it allows us to
     * handle situations with zero CE in a natural and elegant way:
     * all singular values are zero, and even W[0] itself is not
     * greater than W[0]*tol.
     */
    if( state->cecnt>0 )
    {
        b = rmatrixsvd(&state->ce, state->cecnt, n, 1, 1, 2, &state->w, &state->cesvl, &state->cebasis, _state);
        ae_assert(b, "MinBLEIC: inconvergence of internal SVD", _state);
        state->cedim = 0;
        m = ae_minint(state->cecnt, n, _state);
        for(i=0; i<=m-1; i++)
        {
            if( ae_fp_greater(state->w.ptr.p_double[i],state->w.ptr.p_double[0]*minbleic_svdtol*ae_machineepsilon) )
            {
                state->cedim = state->cedim+1;
            }
        }
    }
    else
    {
        state->cedim = 0;
    }
    
    /*
     * Calculate XE: solution of CE*x = b.
     * Fill it with zeros if CEDim=0
     */
    if( state->cedim>0 )
    {
        rvectorsetlengthatleast(&state->tmp0, state->cedim, _state);
        for(i=0; i<=state->cedim-1; i++)
        {
            state->tmp0.ptr.p_double[i] = 0;
        }
        for(i=0; i<=state->cecnt-1; i++)
        {
            v = state->ce.ptr.pp_double[i][n];
            ae_v_addd(&state->tmp0.ptr.p_double[0], 1, &state->cesvl.ptr.pp_double[i][0], 1, ae_v_len(0,state->cedim-1), v);
        }
        for(i=0; i<=state->cedim-1; i++)
        {
            state->tmp0.ptr.p_double[i] = state->tmp0.ptr.p_double[i]/state->w.ptr.p_double[i];
        }
        for(i=0; i<=n-1; i++)
        {
            state->xe.ptr.p_double[i] = 0;
        }
        for(i=0; i<=state->cedim-1; i++)
        {
            v = state->tmp0.ptr.p_double[i];
            ae_v_addd(&state->xe.ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    else
    {
        
        /*
         * no constraints, fill with zeros
         */
        for(i=0; i<=n-1; i++)
        {
            state->xe.ptr.p_double[i] = 0;
        }
    }
}


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
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinBLEICSetInnerCond: EpsG is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinBLEICSetInnerCond: negative EpsG", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinBLEICSetInnerCond: EpsF is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinBLEICSetInnerCond: negative EpsF", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinBLEICSetInnerCond: EpsX is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinBLEICSetInnerCond: negative EpsX", _state);
    state->innerepsg = epsg;
    state->innerepsf = epsf;
    state->innerepsx = epsx;
}


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
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsx, _state), "MinBLEICSetOuterCond: EpsX is not finite number", _state);
    ae_assert(ae_fp_greater(epsx,0), "MinBLEICSetOuterCond: non-positive EpsX", _state);
    ae_assert(ae_isfinite(epsi, _state), "MinBLEICSetOuterCond: EpsI is not finite number", _state);
    ae_assert(ae_fp_greater(epsi,0), "MinBLEICSetOuterCond: non-positive EpsI", _state);
    state->outerepsx = epsx;
    state->outerepsi = epsi;
}


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
     ae_state *_state)
{


    ae_assert(ae_isfinite(mu, _state), "MinBLEICSetBarrierWidth: Mu is not finite number", _state);
    ae_assert(ae_fp_greater(mu,0), "MinBLEICSetBarrierWidth: non-positive Mu", _state);
    state->mustart = mu;
}


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
     ae_state *_state)
{


    ae_assert(ae_isfinite(mudecay, _state), "MinBLEICSetBarrierDecay: MuDecay is not finite number", _state);
    ae_assert(ae_fp_greater(mudecay,0), "MinBLEICSetBarrierDecay: non-positive MuDecay", _state);
    ae_assert(ae_fp_less_eq(mudecay,1), "MinBLEICSetBarrierDecay: MuDecay>1", _state);
    state->mudecay = mudecay;
}


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
     ae_state *_state)
{


    ae_assert(maxits>=0, "MinBLEICSetCond: negative MaxIts!", _state);
    state->maxits = maxits;
}


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
     ae_state *_state)
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
large  steps  which  lead   to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetstpmax(minbleicstate* state,
     double stpmax,
     ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinBLEICSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinBLEICSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool minbleiciteration(minbleicstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    double v;
    double vv;
    ae_bool b;
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
        i = state->rstate.ia.ptr.p_int[2];
        b = state->rstate.ba.ptr.p_bool[0];
        v = state->rstate.ra.ptr.p_double[0];
        vv = state->rstate.ra.ptr.p_double[1];
    }
    else
    {
        n = -983;
        m = -989;
        i = -834;
        b = ae_false;
        v = -287;
        vv = 364;
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
    
    /*
     * Routine body
     */
    
    /*
     * Prepare.
     */
    n = state->n;
    state->repterminationtype = 0;
    state->repinneriterationscount = 0;
    state->repouteriterationscount = 0;
    state->repnfev = 0;
    state->repdebugeqerr = 0.0;
    state->repdebugfs = _state->v_nan;
    state->repdebugff = _state->v_nan;
    state->repdebugdx = _state->v_nan;
    rvectorsetlengthatleast(&state->r, n, _state);
    rvectorsetlengthatleast(&state->tmp1, n, _state);
    
    /*
     * Check that equality constraints are consistent within EpsC.
     * If not - premature termination.
     */
    if( state->cedim>0 )
    {
        state->repdebugeqerr = 0.0;
        for(i=0; i<=state->cecnt-1; i++)
        {
            v = ae_v_dotproduct(&state->ce.ptr.pp_double[i][0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->repdebugeqerr = state->repdebugeqerr+ae_sqr(v-state->ce.ptr.pp_double[i][n], _state);
        }
        state->repdebugeqerr = ae_sqrt(state->repdebugeqerr, _state);
        if( ae_fp_greater(state->repdebugeqerr,state->outerepsi) )
        {
            state->repterminationtype = -3;
            result = ae_false;
            return result;
        }
    }
    
    /*
     * Find feasible point.
     *
     * We make up to 16*N iterations of nonlinear CG trying to
     * minimize penalty for violation of inequality constraints.
     *
     * if P is magnitude of violation, then penalty function has form:
     * * P*P+P for non-zero violation (P>=0)
     * * 0.0 for absence of violation (P<0)
     * Such function is non-smooth at P=0.0, but its nonsmoothness
     * allows us to rapidly converge to the feasible point.
     */
    minbleic_makeprojection(state, &state->xcur, &state->r, &vv, _state);
    mincgrestartfrom(&state->cgstate, &state->xcur, _state);
    mincgsetcond(&state->cgstate, 0.0, 0.0, 0.0, 16*n, _state);
    mincgsetxrep(&state->cgstate, ae_false, _state);
    while(mincgiteration(&state->cgstate, _state))
    {
        if( state->cgstate.needfg )
        {
            ae_v_move(&state->x.ptr.p_double[0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
            minbleic_penaltyfunction(state, &state->x, &state->cgstate.f, &state->cgstate.g, &state->r, &state->mba, &state->errfeas, _state);
            continue;
        }
    }
    mincgresults(&state->cgstate, &state->xcur, &state->cgrep, _state);
    ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_penaltyfunction(state, &state->tmp1, &state->f, &state->g, &state->r, &state->mba, &state->errfeas, _state);
    if( ae_fp_greater(state->errfeas,state->outerepsi) )
    {
        state->repterminationtype = -3;
        result = ae_false;
        return result;
    }
    
    /*
     * Initialize RepDebugFS with function value at initial point
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    state->repdebugfs = state->f;
    
    /*
     * Calculate number of inequality constraints and allocate
     * array for Lagrange multipliers.
     *
     * Initialize Lagrange multipliers and penalty term
     */
    state->lmcnt = state->cicnt;
    for(i=0; i<=n-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            state->lmcnt = state->lmcnt+1;
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            state->lmcnt = state->lmcnt+1;
        }
    }
    if( state->lmcnt>0 )
    {
        rvectorsetlengthatleast(&state->lm, state->lmcnt, _state);
    }
    for(i=0; i<=state->lmcnt-1; i++)
    {
        state->lm.ptr.p_double[i] = 1.0;
    }
    state->mu = ae_maxreal(state->mustart, state->outerepsi, _state);
    
    /*
     * BndMax is a maximum over right parts of inequality constraints.
     * It is used later to bound Mu from below.
     */
    state->bndmax = 0.0;
    for(i=0; i<=n-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            state->bndmax = ae_maxreal(state->bndmax, ae_fabs(state->bndl.ptr.p_double[i], _state), _state);
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            state->bndmax = ae_maxreal(state->bndmax, ae_fabs(state->bndu.ptr.p_double[i], _state), _state);
        }
    }
    for(i=0; i<=state->cicnt-1; i++)
    {
        state->bndmax = ae_maxreal(state->bndmax, ae_fabs(state->ci.ptr.pp_double[i][n], _state), _state);
    }
    
    /*
     * External cycle: optimization subject to current
     * estimate of Lagrange multipliers and penalty term
     */
    state->itsleft = state->maxits;
    state->mucounter = minbleic_mucountdownlen;
    ae_v_move(&state->xprev.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->xstart.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
lbl_5:
    if( ae_false )
    {
        goto lbl_6;
    }
    
    /*
     * Inner cycle: CG with projections and penalty functions
     */
    mincgrestartfrom(&state->cgstate, &state->xcur, _state);
    mincgsetcond(&state->cgstate, state->innerepsg, state->innerepsf, state->innerepsx, state->itsleft, _state);
    mincgsetxrep(&state->cgstate, state->xrep, _state);
    mincgsetdrep(&state->cgstate, ae_true, _state);
    mincgsetstpmax(&state->cgstate, state->stpmax, _state);
lbl_7:
    if( !mincgiteration(&state->cgstate, _state) )
    {
        goto lbl_8;
    }
    if( state->cgstate.lsstart )
    {
        
        /*
         * Beginning of the line search: set upper limit on step size
         * to prevent algo from leaving area where barrier function
         * is defined.
         *
         * We calculate State.CGState.CurStpMax in two steps:
         * * first, we calculate it as the distance from the current
         *   point to the boundary where modified barrier function
         *   overflows; distance is taken along State.CGState.D
         * * then we multiply it by 0.999 to make sure that we won't
         *   make step into the boundary due to the rounding noise
         */
        if( ae_fp_eq(state->cgstate.curstpmax,0) )
        {
            state->cgstate.curstpmax = 1.0E50;
        }
        state->boundary = -0.9*state->mu;
        state->closetobarrier = ae_false;
        for(i=0; i<=n-1; i++)
        {
            if( state->hasbndl.ptr.p_bool[i] )
            {
                v = state->cgstate.x.ptr.p_double[i]-state->bndl.ptr.p_double[i];
                if( ae_fp_less(state->cgstate.d.ptr.p_double[i],0) )
                {
                    state->cgstate.curstpmax = safeminposrv(v+state->mu, -state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                }
                if( ae_fp_greater(state->cgstate.d.ptr.p_double[i],0)&&ae_fp_less_eq(v,state->boundary) )
                {
                    state->cgstate.curstpmax = safeminposrv(-v, state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                    state->closetobarrier = ae_true;
                }
            }
            if( state->hasbndu.ptr.p_bool[i] )
            {
                v = state->bndu.ptr.p_double[i]-state->cgstate.x.ptr.p_double[i];
                if( ae_fp_greater(state->cgstate.d.ptr.p_double[i],0) )
                {
                    state->cgstate.curstpmax = safeminposrv(v+state->mu, state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                }
                if( ae_fp_less(state->cgstate.d.ptr.p_double[i],0)&&ae_fp_less_eq(v,state->boundary) )
                {
                    state->cgstate.curstpmax = safeminposrv(-v, -state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                    state->closetobarrier = ae_true;
                }
            }
        }
        for(i=0; i<=state->cicnt-1; i++)
        {
            v = ae_v_dotproduct(&state->ci.ptr.pp_double[i][0], 1, &state->cgstate.d.ptr.p_double[0], 1, ae_v_len(0,n-1));
            vv = ae_v_dotproduct(&state->ci.ptr.pp_double[i][0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
            vv = vv-state->ci.ptr.pp_double[i][n];
            if( ae_fp_less(v,0) )
            {
                state->cgstate.curstpmax = safeminposrv(vv+state->mu, -v, state->cgstate.curstpmax, _state);
            }
            if( ae_fp_greater(v,0)&&ae_fp_less_eq(vv,state->boundary) )
            {
                state->cgstate.curstpmax = safeminposrv(-vv, v, state->cgstate.curstpmax, _state);
                state->closetobarrier = ae_true;
            }
        }
        state->cgstate.curstpmax = 0.999*state->cgstate.curstpmax;
        if( state->closetobarrier )
        {
            state->cgstate.stp = 0.5*state->cgstate.curstpmax;
        }
        goto lbl_7;
    }
    if( !state->cgstate.needfg )
    {
        goto lbl_9;
    }
    
    /*
     * * get X from CG
     * * project X into equality constrained subspace
     * * RComm (note: X is stored in Tmp1 to prevent accidental corruption by user)
     * * modify target function with barriers and penalties
     * * pass F/G back to nonlinear CG
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_makeprojection(state, &state->x, &state->r, &vv, _state);
    ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->needfg = ae_false;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->tmp1.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_modifytargetfunction(state, &state->x, &state->r, vv, &state->f, &state->g, &state->gnorm, &state->mpgnorm, &state->mba, &state->errfeas, &state->errslack, _state);
    state->cgstate.f = state->f;
    ae_v_move(&state->cgstate.g.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    goto lbl_7;
lbl_9:
    if( !state->cgstate.xupdated )
    {
        goto lbl_11;
    }
    
    /*
     * Report
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->cgstate.f;
    minbleic_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->xupdated = ae_false;
    goto lbl_7;
lbl_11:
    goto lbl_7;
lbl_8:
    mincgresults(&state->cgstate, &state->xcur, &state->cgrep, _state);
    state->repinneriterationscount = state->repinneriterationscount+state->cgrep.iterationscount;
    state->repouteriterationscount = state->repouteriterationscount+1;
    state->repnfev = state->repnfev+state->cgrep.nfev;
    
    /*
     * Update RepDebugFF with function value at current point
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    state->repdebugff = state->f;
    
    /*
     * Update Lagrange multipliers and Mu
     *
     * We calculate three values during update:
     * * LMDif      -   absolute differense between new and old multipliers
     * * LMNorm     -   inf-norm of Lagrange vector
     * * LMGrowth   -   maximum componentwise relative growth of Lagrange vector
     *
     * We limit growth of Lagrange multipliers by MaxLMGrowth,
     * it allows us to stabilize algorithm when it moves deep into
     * the infeasible area.
     *
     * We calculate modified target function at XCur in order to get
     * information about overall problem properties. Some values
     * calculated here will be used later:
     * * ErrFeas  - feasibility error
     * * ErrSlack - complementary slackness error
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_makeprojection(state, &state->x, &state->r, &vv, _state);
    ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->needfg = ae_false;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->tmp1.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minbleic_modifytargetfunction(state, &state->x, &state->r, vv, &state->f, &state->g, &state->gnorm, &state->mpgnorm, &state->mba, &state->errfeas, &state->errslack, _state);
    m = 0;
    state->lmdif = 0;
    state->lmnorm = 0;
    state->lmgrowth = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            barrierfunc(state->xcur.ptr.p_double[i]-state->bndl.ptr.p_double[i], state->mu, &state->v0, &state->v1, &state->v2, _state);
            v = state->lm.ptr.p_double[m];
            vv = ae_minreal(minbleic_maxlmgrowth, -state->mu*state->v1, _state);
            state->lm.ptr.p_double[m] = vv*state->lm.ptr.p_double[m];
            state->lmnorm = ae_maxreal(state->lmnorm, ae_fabs(v, _state), _state);
            state->lmdif = ae_maxreal(state->lmdif, ae_fabs(state->lm.ptr.p_double[m]-v, _state), _state);
            state->lmgrowth = ae_maxreal(state->lmgrowth, vv, _state);
            m = m+1;
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            barrierfunc(state->bndu.ptr.p_double[i]-state->xcur.ptr.p_double[i], state->mu, &state->v0, &state->v1, &state->v2, _state);
            v = state->lm.ptr.p_double[m];
            vv = ae_minreal(minbleic_maxlmgrowth, -state->mu*state->v1, _state);
            state->lm.ptr.p_double[m] = vv*state->lm.ptr.p_double[m];
            state->lmnorm = ae_maxreal(state->lmnorm, ae_fabs(v, _state), _state);
            state->lmdif = ae_maxreal(state->lmdif, ae_fabs(state->lm.ptr.p_double[m]-v, _state), _state);
            state->lmgrowth = ae_maxreal(state->lmgrowth, vv, _state);
            m = m+1;
        }
    }
    for(i=0; i<=state->cicnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ci.ptr.pp_double[i][0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = v-state->ci.ptr.pp_double[i][n];
        barrierfunc(v, state->mu, &state->v0, &state->v1, &state->v2, _state);
        v = state->lm.ptr.p_double[m];
        vv = ae_minreal(minbleic_maxlmgrowth, -state->mu*state->v1, _state);
        state->lm.ptr.p_double[m] = vv*state->lm.ptr.p_double[m];
        state->lmnorm = ae_maxreal(state->lmnorm, ae_fabs(v, _state), _state);
        state->lmdif = ae_maxreal(state->lmdif, ae_fabs(state->lm.ptr.p_double[m]-v, _state), _state);
        state->lmgrowth = ae_maxreal(state->lmgrowth, vv, _state);
        m = m+1;
    }
    if( ae_fp_greater(state->mba,-0.2*state->mudecay*state->mu) )
    {
        state->mu = ae_maxreal(state->mudecay*state->mu, 1.0E6*ae_machineepsilon*state->bndmax, _state);
    }
    if( ae_fp_less(state->mu,state->outerepsi) )
    {
        state->mucounter = state->mucounter-1;
        state->mu = state->outerepsi;
    }
    
    /*
     * Check for stopping:
     * * "normal", outer step size is small enough, infeasibility is within bounds
     * * "inconsistent",  if Lagrange multipliers increased beyond threshold given by MaxLagrangeMul
     * * "too stringent", in other cases
     */
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->xcur.ptr.p_double[i]-state->xprev.ptr.p_double[i], _state);
    }
    v = ae_sqrt(v, _state);
    if( ae_fp_less_eq(state->errfeas,state->outerepsi)&&ae_fp_less_eq(v,state->outerepsx) )
    {
        state->repterminationtype = 4;
        goto lbl_6;
    }
    if( state->maxits>0 )
    {
        state->itsleft = state->itsleft-state->cgrep.iterationscount;
        if( state->itsleft<=0 )
        {
            state->repterminationtype = 5;
            goto lbl_6;
        }
    }
    if( ae_fp_greater_eq(state->repouteriterationscount,minbleic_maxouterits) )
    {
        state->repterminationtype = 5;
        goto lbl_6;
    }
    if( ae_fp_less(state->lmdif,minbleic_lmtol*ae_machineepsilon*state->lmnorm)||ae_fp_less(state->lmnorm,minbleic_minlagrangemul) )
    {
        state->repterminationtype = 7;
        goto lbl_6;
    }
    if( state->mucounter<=0 )
    {
        state->repterminationtype = 7;
        goto lbl_6;
    }
    if( ae_fp_greater(state->lmnorm,minbleic_maxlagrangemul) )
    {
        state->repterminationtype = -3;
        goto lbl_6;
    }
    
    /*
     * Next iteration
     */
    ae_v_move(&state->xprev.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,n-1));
    goto lbl_5;
lbl_6:
    
    /*
     * We've stopped, fill debug information
     */
    state->repdebugeqerr = 0.0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ce.ptr.pp_double[i][0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->repdebugeqerr = state->repdebugeqerr+ae_sqr(v-state->ce.ptr.pp_double[i][n], _state);
    }
    state->repdebugeqerr = ae_sqrt(state->repdebugeqerr, _state);
    state->repdebugdx = 0;
    for(i=0; i<=n-1; i++)
    {
        state->repdebugdx = state->repdebugdx+ae_sqr(state->xcur.ptr.p_double[i]-state->xstart.ptr.p_double[i], _state);
    }
    state->repdebugdx = ae_sqrt(state->repdebugdx, _state);
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = i;
    state->rstate.ba.ptr.p_bool[0] = b;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = vv;
    return result;
}


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
     ae_state *_state)
{

    ae_vector_clear(x);
    _minbleicreport_clear(rep);

    minbleicresultsbuf(state, x, rep, _state);
}


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
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->inneriterationscount = state->repinneriterationscount;
    rep->outeriterationscount = state->repouteriterationscount;
    rep->nfev = state->repnfev;
    rep->terminationtype = state->repterminationtype;
    rep->debugeqerr = state->repdebugeqerr;
    rep->debugfs = state->repdebugfs;
    rep->debugff = state->repdebugff;
    rep->debugdx = state->repdebugdx;
}


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
     ae_state *_state)
{
    ae_int_t n;


    n = state->n;
    
    /*
     * First, check for errors in the inputs
     */
    ae_assert(x->cnt>=n, "MinBLEICRestartFrom: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinBLEICRestartFrom: X contains infinite or NaN values!", _state);
    
    /*
     * Set XC
     */
    ae_v_move(&state->xcur.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * prepare RComm facilities
     */
    ae_vector_set_length(&state->rstate.ia, 2+1, _state);
    ae_vector_set_length(&state->rstate.ba, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
    minbleic_clearrequestfields(state, _state);
}


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
     ae_state *_state)
{
    double c0;
    double c1;
    double c2;
    double xpmu;

    *f = 0;
    *df = 0;
    *d2f = 0;

    xpmu = x+mu;
    if( ae_fp_greater(xpmu,0.5*mu) )
    {
        *f = -ae_log(x/mu+1, _state);
        *df = -1/(x+mu);
        *d2f = 1/(xpmu*xpmu);
        return;
    }
    if( ae_fp_greater(xpmu,0) )
    {
        c0 = -ae_log(0.5, _state)-0.5;
        c1 = -1/mu;
        c2 = mu/4;
        *f = c0+c1*(x+0.5*mu)+c2/xpmu;
        *df = c1-c2/((x+mu)*(x+mu));
        *d2f = 2*c2/((x+mu)*(x+mu)*(x+mu));
        return;
    }
    *f = ae_maxrealnumber;
    *df = 0;
    *d2f = 0;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forget to clear something)
*************************************************************************/
static void minbleic_clearrequestfields(minbleicstate* state,
     ae_state *_state)
{


    state->needfg = ae_false;
    state->xupdated = ae_false;
}


/*************************************************************************
This function makes projection of X into equality constrained subspace.
It calculates set of additional values which are used later for
modification of the target function F.

INPUT PARAMETERS:
    State   -   optimizer state (we use its fields to get information
                about constraints)
    X       -   vector being projected
    R       -   preallocated buffer, used to store residual from projection

OUTPUT PARAMETERS:
    X       -   projection of input X
    R       -   residual
    RNorm   -   residual norm squared, used later to modify target function
*************************************************************************/
static void minbleic_makeprojection(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* r,
     double* rnorm2,
     ae_state *_state)
{
    double v;
    ae_int_t i;
    ae_int_t n;

    *rnorm2 = 0;

    n = state->n;
    
    /*
     * * subtract XE from X
     * * project X
     * * calculate norm of deviation from null space, store it in VV
     * * calculate residual from projection, store it in R
     * * add XE to X
     */
    ae_v_sub(&x->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
    *rnorm2 = 0;
    for(i=0; i<=n-1; i++)
    {
        r->ptr.p_double[i] = 0;
    }
    for(i=0; i<=state->cedim-1; i++)
    {
        v = ae_v_dotproduct(&x->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        ae_v_subd(&x->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        ae_v_addd(&r->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        *rnorm2 = *rnorm2+ae_sqr(v, _state);
    }
    ae_v_add(&x->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
}


/*************************************************************************
This subroutine applies modifications to the target function given by
its value F and gradient G at the projected point X which lies in the
equality constrained subspace.

Following modifications are applied:
* modified barrier functions to handle inequality constraints
  (both F and G are modified)
* projection of gradient into equality constrained subspace
  (only G is modified)
* quadratic penalty for deviations from equality constrained subspace
  (both F and G are modified)

It also calculates gradient norm (three different norms for three
different types of gradient), feasibility and complementary slackness
errors.

INPUT PARAMETERS:
    State   -   optimizer state (we use its fields to get information
                about constraints)
    X       -   point (projected into equality constrained subspace)
    R       -   residual from projection
    RNorm2  -   residual norm squared
    F       -   function value at X
    G       -   function gradient at X

OUTPUT PARAMETERS:
    F       -   modified function value at X
    G       -   modified function gradient at X
    GNorm   -   2-norm of unmodified G
    MPGNorm -   2-norm of modified G
    MBA     -   minimum argument of barrier functions.
                If X is strictly feasible, it is greater than zero.
                If X lies on a boundary, it is zero.
                It is negative for infeasible X.
    FIErr   -   2-norm of feasibility error with respect to
                inequality/bound constraints
    CSErr   -   2-norm of complementarity slackness error
*************************************************************************/
static void minbleic_modifytargetfunction(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* r,
     double rnorm2,
     double* f,
     /* Real    */ ae_vector* g,
     double* gnorm,
     double* mpgnorm,
     double* mba,
     double* fierr,
     double* cserr,
     ae_state *_state)
{
    double v;
    double vv;
    double t;
    ae_int_t i;
    ae_int_t n;
    ae_int_t m;
    double v0;
    double v1;
    double v2;
    ae_bool hasconstraints;

    *gnorm = 0;
    *mpgnorm = 0;
    *mba = 0;
    *fierr = 0;
    *cserr = 0;

    n = state->n;
    *mba = ae_maxrealnumber;
    hasconstraints = ae_false;
    
    /*
     * GNorm
     */
    v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &g->ptr.p_double[0], 1, ae_v_len(0,n-1));
    *gnorm = ae_sqrt(v, _state);
    
    /*
     * Process bound and inequality constraints.
     * Bound constraints with +-INF are ignored.
     * Here M is used to store number of constraints processed.
     */
    m = 0;
    *fierr = 0;
    *cserr = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            v = x->ptr.p_double[i]-state->bndl.ptr.p_double[i];
            *mba = ae_minreal(v, *mba, _state);
            barrierfunc(v, state->mu, &v0, &v1, &v2, _state);
            *f = *f+state->mu*state->lm.ptr.p_double[m]*v0;
            g->ptr.p_double[i] = g->ptr.p_double[i]+state->mu*state->lm.ptr.p_double[m]*v1;
            if( ae_fp_less(v,0) )
            {
                *fierr = *fierr+v*v;
            }
            t = -state->lm.ptr.p_double[m]*v;
            *cserr = *cserr+t*t;
            m = m+1;
            hasconstraints = ae_true;
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            v = state->bndu.ptr.p_double[i]-x->ptr.p_double[i];
            *mba = ae_minreal(v, *mba, _state);
            barrierfunc(v, state->mu, &v0, &v1, &v2, _state);
            *f = *f+state->mu*state->lm.ptr.p_double[m]*v0;
            g->ptr.p_double[i] = g->ptr.p_double[i]-state->mu*state->lm.ptr.p_double[m]*v1;
            if( ae_fp_less(v,0) )
            {
                *fierr = *fierr+v*v;
            }
            t = -state->lm.ptr.p_double[m]*v;
            *cserr = *cserr+t*t;
            m = m+1;
            hasconstraints = ae_true;
        }
    }
    for(i=0; i<=state->cicnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ci.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = v-state->ci.ptr.pp_double[i][n];
        *mba = ae_minreal(v, *mba, _state);
        barrierfunc(v, state->mu, &v0, &v1, &v2, _state);
        *f = *f+state->mu*state->lm.ptr.p_double[m]*v0;
        vv = state->mu*state->lm.ptr.p_double[m]*v1;
        ae_v_addd(&g->ptr.p_double[0], 1, &state->ci.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), vv);
        if( ae_fp_less(v,0) )
        {
            *fierr = *fierr+v*v;
        }
        t = -state->lm.ptr.p_double[m]*v;
        *cserr = *cserr+t*t;
        m = m+1;
        hasconstraints = ae_true;
    }
    *fierr = ae_sqrt(*fierr, _state);
    *cserr = ae_sqrt(*cserr, _state);
    if( !hasconstraints )
    {
        *mba = 0.0;
    }
    
    /*
     * Process equality constraints:
     * * modify F to handle penalty term for equality constraints
     * * project gradient on null space of equality constraints
     * * add penalty term for equality constraints to gradient
     */
    *f = *f+rnorm2;
    for(i=0; i<=state->cedim-1; i++)
    {
        v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        ae_v_subd(&g->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
    }
    ae_v_addd(&g->ptr.p_double[0], 1, &r->ptr.p_double[0], 1, ae_v_len(0,n-1), 2);
    
    /*
     * MPGNorm
     */
    v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &g->ptr.p_double[0], 1, ae_v_len(0,n-1));
    *mpgnorm = ae_sqrt(v, _state);
}


/*************************************************************************
This subroutine calculates penalty for violation of equality/inequality
constraints. It is used to find feasible point.

Following modifications are applied:
* quadratic penalty for deviations from inequality constrained subspace
  (both F and G are modified)
* projection of gradient into equality constrained subspace
  (only G is modified)
* quadratic penalty for deviations from equality constrained subspace
  (both F and G are modified)

INPUT PARAMETERS:
    State   -   optimizer state (we use its fields to get information
                about constraints)
    X       -   point (modified by function)
    G       -   preallocated array[N]
    R       -   preallocated array[N]

OUTPUT PARAMETERS:
    X       -   projection of X into equality constrained subspace
    F       -   modified function value at X
    G       -   modified function gradient at X
    MBA     -   minimum argument of barrier functions.
                If X is strictly feasible, it is greater than zero.
                If X lies on a boundary, it is zero.
                It is negative for infeasible X.
    FIErr   -   2-norm of feasibility error with respect to
                inequality/bound constraints
*************************************************************************/
static void minbleic_penaltyfunction(minbleicstate* state,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* r,
     double* mba,
     double* fierr,
     ae_state *_state)
{
    double v;
    double vv;
    ae_int_t i;
    ae_int_t m;
    ae_int_t n;
    double rnorm2;
    ae_bool hasconstraints;

    *mba = 0;
    *fierr = 0;

    n = state->n;
    *mba = ae_maxrealnumber;
    hasconstraints = ae_false;
    
    /*
     * Initialize F/G
     */
    *f = 0.0;
    for(i=0; i<=n-1; i++)
    {
        g->ptr.p_double[i] = 0.0;
    }
    
    /*
     * Calculate projection of X:
     * * subtract XE from X
     * * project X
     * * calculate norm of deviation from null space, store it in VV
     * * calculate residual from projection, store it in R
     * * add XE to X
     */
    ae_v_sub(&x->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
    rnorm2 = 0;
    for(i=0; i<=n-1; i++)
    {
        r->ptr.p_double[i] = 0;
    }
    for(i=0; i<=state->cedim-1; i++)
    {
        v = ae_v_dotproduct(&x->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        ae_v_subd(&x->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        ae_v_addd(&r->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        rnorm2 = rnorm2+ae_sqr(v, _state);
    }
    ae_v_add(&x->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Process bound and inequality constraints.
     * Bound constraints with +-INF are ignored.
     * Here M is used to store number of constraints processed.
     */
    m = 0;
    *fierr = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            v = x->ptr.p_double[i]-state->bndl.ptr.p_double[i];
            *mba = ae_minreal(v, *mba, _state);
            if( ae_fp_less(v,0) )
            {
                *f = *f+v*v;
                g->ptr.p_double[i] = g->ptr.p_double[i]+2*v;
                *fierr = *fierr+v*v;
            }
            m = m+1;
            hasconstraints = ae_true;
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            v = state->bndu.ptr.p_double[i]-x->ptr.p_double[i];
            *mba = ae_minreal(v, *mba, _state);
            if( ae_fp_less(v,0) )
            {
                *f = *f+v*v;
                g->ptr.p_double[i] = g->ptr.p_double[i]-2*v;
                *fierr = *fierr+v*v;
            }
            m = m+1;
            hasconstraints = ae_true;
        }
    }
    for(i=0; i<=state->cicnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ci.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = v-state->ci.ptr.pp_double[i][n];
        *mba = ae_minreal(v, *mba, _state);
        if( ae_fp_less(v,0) )
        {
            *f = *f+v*v;
            vv = 2*v;
            ae_v_addd(&g->ptr.p_double[0], 1, &state->ci.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), vv);
            *fierr = *fierr+v*v;
        }
        m = m+1;
        hasconstraints = ae_true;
    }
    *fierr = ae_sqrt(*fierr, _state);
    if( !hasconstraints )
    {
        *mba = 0.0;
    }
    
    /*
     * Process equality constraints:
     * * modify F to handle penalty term for equality constraints
     * * project gradient on null space of equality constraints
     * * add penalty term for equality constraints to gradient
     */
    *f = *f+rnorm2;
    for(i=0; i<=state->cedim-1; i++)
    {
        v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        ae_v_subd(&g->ptr.p_double[0], 1, &state->cebasis.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
    }
    ae_v_addd(&g->ptr.p_double[0], 1, &r->ptr.p_double[0], 1, ae_v_len(0,n-1), 2);
}


ae_bool _minbleicstate_init(minbleicstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xcur, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xprev, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xstart, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ce, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ci, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->cebasis, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->cesvl, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xe, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hasbndl, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hasbndu, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lm, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->w, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->r, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_linminstate_init(&p->lstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgstate_init(&p->cgstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgreport_init(&p->cgrep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minbleicstate_init_copy(minbleicstate* dst, minbleicstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->innerepsg = src->innerepsg;
    dst->innerepsf = src->innerepsf;
    dst->innerepsx = src->innerepsx;
    dst->outerepsx = src->outerepsx;
    dst->outerepsi = src->outerepsi;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->cgtype = src->cgtype;
    dst->mustart = src->mustart;
    dst->mudecay = src->mudecay;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repinneriterationscount = src->repinneriterationscount;
    dst->repouteriterationscount = src->repouteriterationscount;
    dst->repnfev = src->repnfev;
    dst->repterminationtype = src->repterminationtype;
    dst->repdebugeqerr = src->repdebugeqerr;
    dst->repdebugfs = src->repdebugfs;
    dst->repdebugff = src->repdebugff;
    dst->repdebugdx = src->repdebugdx;
    if( !ae_vector_init_copy(&dst->xcur, &src->xcur, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xprev, &src->xprev, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xstart, &src->xstart, _state, make_automatic) )
        return ae_false;
    dst->itsleft = src->itsleft;
    dst->mucounter = src->mucounter;
    if( !ae_matrix_init_copy(&dst->ce, &src->ce, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->ci, &src->ci, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->cebasis, &src->cebasis, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->cesvl, &src->cesvl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xe, &src->xe, _state, make_automatic) )
        return ae_false;
    dst->cecnt = src->cecnt;
    dst->cicnt = src->cicnt;
    dst->cedim = src->cedim;
    if( !ae_vector_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hasbndl, &src->hasbndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hasbndu, &src->hasbndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lm, &src->lm, _state, make_automatic) )
        return ae_false;
    dst->lmcnt = src->lmcnt;
    dst->mu = src->mu;
    if( !ae_vector_init_copy(&dst->w, &src->w, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp1, &src->tmp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->r, &src->r, _state, make_automatic) )
        return ae_false;
    dst->v0 = src->v0;
    dst->v1 = src->v1;
    dst->v2 = src->v2;
    dst->t = src->t;
    dst->errfeas = src->errfeas;
    dst->errslack = src->errslack;
    dst->gnorm = src->gnorm;
    dst->mpgnorm = src->mpgnorm;
    dst->lmdif = src->lmdif;
    dst->lmnorm = src->lmnorm;
    dst->lmgrowth = src->lmgrowth;
    dst->mba = src->mba;
    dst->boundary = src->boundary;
    dst->closetobarrier = src->closetobarrier;
    dst->bndmax = src->bndmax;
    if( !_linminstate_init_copy(&dst->lstate, &src->lstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgstate_init_copy(&dst->cgstate, &src->cgstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgreport_init_copy(&dst->cgrep, &src->cgrep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _minbleicstate_clear(minbleicstate* p)
{
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    ae_vector_clear(&p->xcur);
    ae_vector_clear(&p->xprev);
    ae_vector_clear(&p->xstart);
    ae_matrix_clear(&p->ce);
    ae_matrix_clear(&p->ci);
    ae_matrix_clear(&p->cebasis);
    ae_matrix_clear(&p->cesvl);
    ae_vector_clear(&p->xe);
    ae_vector_clear(&p->bndl);
    ae_vector_clear(&p->hasbndl);
    ae_vector_clear(&p->bndu);
    ae_vector_clear(&p->hasbndu);
    ae_vector_clear(&p->lm);
    ae_vector_clear(&p->w);
    ae_vector_clear(&p->tmp0);
    ae_vector_clear(&p->tmp1);
    ae_vector_clear(&p->r);
    _linminstate_clear(&p->lstate);
    _mincgstate_clear(&p->cgstate);
    _mincgreport_clear(&p->cgrep);
}


ae_bool _minbleicreport_init(minbleicreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _minbleicreport_init_copy(minbleicreport* dst, minbleicreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->debugeqerr = src->debugeqerr;
    dst->debugfs = src->debugfs;
    dst->debugff = src->debugff;
    dst->debugdx = src->debugdx;
    return ae_true;
}


void _minbleicreport_clear(minbleicreport* p)
{
}


/*$ End $*/
