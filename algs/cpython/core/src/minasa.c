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
#include "minasa.h"


/*$ Declarations $*/
static ae_int_t minasa_n1 = 2;
static ae_int_t minasa_n2 = 2;
static double minasa_stpmin = 1.0E-300;
static double minasa_gpaftol = 0.0001;
static double minasa_gpadecay = 0.5;
static double minasa_asarho = 0.5;
static double minasa_asaboundedantigradnorm(minasastate* state,
     ae_state *_state);
static double minasa_asaginorm(minasastate* state, ae_state *_state);
static double minasa_asad1norm(minasastate* state, ae_state *_state);
static ae_bool minasa_asauisempty(minasastate* state, ae_state *_state);
static void minasa_clearrequestfields(minasastate* state,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
              NONLINEAR BOUND CONSTRAINED OPTIMIZATION USING
                      MODIFIED ACTIVE SET ALGORITHM
                   WILLIAM W. HAGER AND HONGCHAO ZHANG

DESCRIPTION:
The  subroutine  minimizes  function  F(x)  of  N  arguments  with   bound
constraints: BndL[i] <= x[i] <= BndU[i]

This method is  globally  convergent  as  long  as  grad(f)  is  Lipschitz
continuous on a level set: L = { x : f(x)<=f(x0) }.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinASACreate() call
2. User tunes solver parameters with MinASASetCond() MinASASetStpMax() and
   other functions
3. User calls MinASAOptimize() function which takes algorithm  state   and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinASAResults() to get solution
5. Optionally, user may call MinASARestartFrom() to solve another  problem
   with same N but another starting point and/or another function.
   MinASARestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from sizes of
                  X/BndL/BndU.
    X       -   starting point, array[0..N-1].
    BndL    -   lower bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very small number as bound: -1000,  -1.0E6
                or -1.0E300, or something like that.
    BndU    -   upper bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very large number as bound: +1000,  +1.0E6
                or +1.0E300, or something like that.

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

NOTES:

1. you may tune stopping conditions with MinASASetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinASASetStpMax() function to bound algorithm's steps.
3. this function does NOT support infinite/NaN values in X, BndL, BndU.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void minasacreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     minasastate* state,
     ae_state *_state)
{
    ae_int_t i;

    _minasastate_clear(state);

    ae_assert(n>=1, "MinASA: N too small!", _state);
    ae_assert(x->cnt>=n, "MinCGCreate: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinCGCreate: X contains infinite or NaN values!", _state);
    ae_assert(bndl->cnt>=n, "MinCGCreate: Length(BndL)<N!", _state);
    ae_assert(isfinitevector(bndl, n, _state), "MinCGCreate: BndL contains infinite or NaN values!", _state);
    ae_assert(bndu->cnt>=n, "MinCGCreate: Length(BndU)<N!", _state);
    ae_assert(isfinitevector(bndu, n, _state), "MinCGCreate: BndU contains infinite or NaN values!", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_fp_less_eq(bndl->ptr.p_double[i],bndu->ptr.p_double[i]), "MinASA: inconsistent bounds!", _state);
        ae_assert(ae_fp_less_eq(bndl->ptr.p_double[i],x->ptr.p_double[i]), "MinASA: infeasible X!", _state);
        ae_assert(ae_fp_less_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]), "MinASA: infeasible X!", _state);
    }
    
    /*
     * Initialize
     */
    state->n = n;
    minasasetcond(state, 0, 0, 0, 0, _state);
    minasasetxrep(state, ae_false, _state);
    minasasetstpmax(state, 0, _state);
    minasasetalgorithm(state, -1, _state);
    ae_vector_set_length(&state->bndl, n, _state);
    ae_vector_set_length(&state->bndu, n, _state);
    ae_vector_set_length(&state->ak, n, _state);
    ae_vector_set_length(&state->xk, n, _state);
    ae_vector_set_length(&state->dk, n, _state);
    ae_vector_set_length(&state->an, n, _state);
    ae_vector_set_length(&state->xn, n, _state);
    ae_vector_set_length(&state->dn, n, _state);
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->d, n, _state);
    ae_vector_set_length(&state->g, n, _state);
    ae_vector_set_length(&state->gc, n, _state);
    ae_vector_set_length(&state->work, n, _state);
    ae_vector_set_length(&state->yk, n, _state);
    minasarestartfrom(state, x, bndl, bndu, _state);
}


/*************************************************************************
This function sets stopping conditions for the ASA optimization algorithm.

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
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetcond(minasastate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinASASetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinASASetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinASASetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinASASetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinASASetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinASASetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinASASetCond: negative MaxIts!", _state);
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
provided to MinASAOptimize().

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetxrep(minasastate* state, ae_bool needxrep, ae_state *_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
This function sets optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm stat
    UAType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetalgorithm(minasastate* state,
     ae_int_t algotype,
     ae_state *_state)
{


    ae_assert(algotype>=-1&&algotype<=1, "MinASASetAlgorithm: incorrect AlgoType!", _state);
    if( algotype==-1 )
    {
        algotype = 1;
    }
    state->cgtype = algotype;
}


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length (zero by default).

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetstpmax(minasastate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinASASetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinASASetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool minasaiteration(minasastate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double betak;
    double v;
    double vv;
    ae_int_t mcinfo;
    ae_bool b;
    ae_bool stepfound;
    ae_int_t diffcnt;
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
        i = state->rstate.ia.ptr.p_int[1];
        mcinfo = state->rstate.ia.ptr.p_int[2];
        diffcnt = state->rstate.ia.ptr.p_int[3];
        b = state->rstate.ba.ptr.p_bool[0];
        stepfound = state->rstate.ba.ptr.p_bool[1];
        betak = state->rstate.ra.ptr.p_double[0];
        v = state->rstate.ra.ptr.p_double[1];
        vv = state->rstate.ra.ptr.p_double[2];
    }
    else
    {
        n = -983;
        i = -989;
        mcinfo = -834;
        diffcnt = 900;
        b = ae_true;
        stepfound = ae_false;
        betak = 214;
        v = -338;
        vv = -686;
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
    
    /*
     * Routine body
     */
    
    /*
     * Prepare
     */
    n = state->n;
    state->repterminationtype = 0;
    state->repiterationscount = 0;
    state->repnfev = 0;
    state->debugrestartscount = 0;
    state->cgtype = 1;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->xk.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xk.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->ak.ptr.p_double[i] = 0;
        }
        else
        {
            state->ak.ptr.p_double[i] = 1;
        }
    }
    state->mu = 0.1;
    state->curalgo = 0;
    
    /*
     * Calculate F/G, initialize algorithm
     */
    minasa_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needfg = ae_false;
    if( !state->xrep )
    {
        goto lbl_15;
    }
    
    /*
     * progress report
     */
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->xupdated = ae_false;
lbl_15:
    if( ae_fp_less_eq(minasa_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    state->repnfev = state->repnfev+1;
    
    /*
     * Main cycle
     *
     * At the beginning of new iteration:
     * * CurAlgo stores current algorithm selector
     * * State.XK, State.F and State.G store current X/F/G
     * * State.AK stores current set of active constraints
     */
lbl_17:
    if( ae_false )
    {
        goto lbl_18;
    }
    
    /*
     * GPA algorithm
     */
    if( state->curalgo!=0 )
    {
        goto lbl_19;
    }
    state->k = 0;
    state->acount = 0;
lbl_21:
    if( ae_false )
    {
        goto lbl_22;
    }
    
    /*
     * Determine Dk = proj(xk - gk)-xk
     */
    for(i=0; i<=n-1; i++)
    {
        state->d.ptr.p_double[i] = boundval(state->xk.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state)-state->xk.ptr.p_double[i];
    }
    
    /*
     * Armijo line search.
     * * exact search with alpha=1 is tried first,
     *   'exact' means that we evaluate f() EXACTLY at
     *   bound(x-g,bndl,bndu), without intermediate floating
     *   point operations.
     * * alpha<1 are tried if explicit search wasn't successful
     * Result is placed into XN.
     *
     * Two types of search are needed because we can't
     * just use second type with alpha=1 because in finite
     * precision arithmetics (x1-x0)+x0 may differ from x1.
     * So while x1 is correctly bounded (it lie EXACTLY on
     * boundary, if it is active), (x1-x0)+x0 may be
     * not bounded.
     */
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->dginit = v;
    state->finit = state->f;
    if( !(ae_fp_less_eq(minasa_asad1norm(state, _state),state->stpmax)||ae_fp_eq(state->stpmax,0)) )
    {
        goto lbl_23;
    }
    
    /*
     * Try alpha=1 step first
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = boundval(state->xk.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    minasa_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    stepfound = ae_fp_less_eq(state->f,state->finit+minasa_gpaftol*state->dginit);
    goto lbl_24;
lbl_23:
    stepfound = ae_false;
lbl_24:
    if( !stepfound )
    {
        goto lbl_25;
    }
    
    /*
     * we are at the boundary(ies)
     */
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->stp = 1;
    goto lbl_26;
lbl_25:
    
    /*
     * alpha=1 is too large, try smaller values
     */
    state->stp = 1;
    linminnormalized(&state->d, &state->stp, n, _state);
    state->dginit = state->dginit/state->stp;
    state->stp = minasa_gpadecay*state->stp;
    if( ae_fp_greater(state->stpmax,0) )
    {
        state->stp = ae_minreal(state->stp, state->stpmax, _state);
    }
lbl_27:
    if( ae_false )
    {
        goto lbl_28;
    }
    v = state->stp;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    minasa_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    if( ae_fp_less_eq(state->stp,minasa_stpmin) )
    {
        goto lbl_28;
    }
    if( ae_fp_less_eq(state->f,state->finit+state->stp*minasa_gpaftol*state->dginit) )
    {
        goto lbl_28;
    }
    state->stp = state->stp*minasa_gpadecay;
    goto lbl_27;
lbl_28:
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
lbl_26:
    state->repiterationscount = state->repiterationscount+1;
    if( !state->xrep )
    {
        goto lbl_29;
    }
    
    /*
     * progress report
     */
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->xupdated = ae_false;
lbl_29:
    
    /*
     * Calculate new set of active constraints.
     * Reset counter if active set was changed.
     * Prepare for the new iteration
     */
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->xn.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xn.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->an.ptr.p_double[i] = 0;
        }
        else
        {
            state->an.ptr.p_double[i] = 1;
        }
    }
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(state->ak.ptr.p_double[i],state->an.ptr.p_double[i]) )
        {
            state->acount = -1;
            break;
        }
    }
    state->acount = state->acount+1;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->ak.ptr.p_double[0], 1, &state->an.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Stopping conditions
     */
    if( !(state->repiterationscount>=state->maxits&&state->maxits>0) )
    {
        goto lbl_31;
    }
    
    /*
     * Too many iterations
     */
    state->repterminationtype = 5;
    if( !state->xrep )
    {
        goto lbl_33;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->xupdated = ae_false;
lbl_33:
    result = ae_false;
    return result;
lbl_31:
    if( ae_fp_greater(minasa_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        goto lbl_35;
    }
    
    /*
     * Gradient is small enough
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_37;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->xupdated = ae_false;
lbl_37:
    result = ae_false;
    return result;
lbl_35:
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( ae_fp_greater(ae_sqrt(v, _state)*state->stp,state->epsx) )
    {
        goto lbl_39;
    }
    
    /*
     * Step size is too small, no further improvement is
     * possible
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_41;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->xupdated = ae_false;
lbl_41:
    result = ae_false;
    return result;
lbl_39:
    if( ae_fp_greater(state->finit-state->f,state->epsf*ae_maxreal(ae_fabs(state->finit, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
    {
        goto lbl_43;
    }
    
    /*
     * F(k+1)-F(k) is small enough
     */
    state->repterminationtype = 1;
    if( !state->xrep )
    {
        goto lbl_45;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->xupdated = ae_false;
lbl_45:
    result = ae_false;
    return result;
lbl_43:
    
    /*
     * Decide - should we switch algorithm or not
     */
    if( minasa_asauisempty(state, _state) )
    {
        if( ae_fp_greater_eq(minasa_asaginorm(state, _state),state->mu*minasa_asad1norm(state, _state)) )
        {
            state->curalgo = 1;
            goto lbl_22;
        }
        else
        {
            state->mu = state->mu*minasa_asarho;
        }
    }
    else
    {
        if( state->acount==minasa_n1 )
        {
            if( ae_fp_greater_eq(minasa_asaginorm(state, _state),state->mu*minasa_asad1norm(state, _state)) )
            {
                state->curalgo = 1;
                goto lbl_22;
            }
        }
    }
    
    /*
     * Next iteration
     */
    state->k = state->k+1;
    goto lbl_21;
lbl_22:
lbl_19:
    
    /*
     * CG algorithm
     */
    if( state->curalgo!=1 )
    {
        goto lbl_47;
    }
    
    /*
     * first, check that there are non-active constraints.
     * move to GPA algorithm, if all constraints are active
     */
    b = ae_true;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(state->ak.ptr.p_double[i],0) )
        {
            b = ae_false;
            break;
        }
    }
    if( b )
    {
        state->curalgo = 0;
        goto lbl_17;
    }
    
    /*
     * CG iterations
     */
    state->fold = state->f;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        state->dk.ptr.p_double[i] = -state->g.ptr.p_double[i]*state->ak.ptr.p_double[i];
        state->gc.ptr.p_double[i] = state->g.ptr.p_double[i]*state->ak.ptr.p_double[i];
    }
lbl_49:
    if( ae_false )
    {
        goto lbl_50;
    }
    
    /*
     * Store G[k] for later calculation of Y[k]
     */
    for(i=0; i<=n-1; i++)
    {
        state->yk.ptr.p_double[i] = -state->gc.ptr.p_double[i];
    }
    
    /*
     * Make a CG step in direction given by DK[]:
     * * calculate step. Step projection into feasible set
     *   is used. It has several benefits: a) step may be
     *   found with usual line search, b) multiple constraints
     *   may be activated with one step, c) activated constraints
     *   are detected in a natural way - just compare x[i] with
     *   bounds
     * * update active set, set B to True, if there
     *   were changes in the set.
     */
    ae_v_move(&state->d.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->mcstage = 0;
    state->stp = 1;
    linminnormalized(&state->d, &state->stp, n, _state);
    if( ae_fp_neq(state->laststep,0) )
    {
        state->stp = state->laststep;
    }
    mcsrch(n, &state->xn, &state->f, &state->gc, &state->d, &state->stp, state->stpmax, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
lbl_51:
    if( state->mcstage==0 )
    {
        goto lbl_52;
    }
    
    /*
     * preprocess data: bound State.XN so it belongs to the
     * feasible set and store it in the State.X
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = boundval(state->xn.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    
    /*
     * RComm
     */
    minasa_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->needfg = ae_false;
    
    /*
     * postprocess data: zero components of G corresponding to
     * the active constraints
     */
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->gc.ptr.p_double[i] = 0;
        }
        else
        {
            state->gc.ptr.p_double[i] = state->g.ptr.p_double[i];
        }
    }
    mcsrch(n, &state->xn, &state->f, &state->gc, &state->d, &state->stp, state->stpmax, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
    goto lbl_51;
lbl_52:
    diffcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * XN contains unprojected result, project it,
         * save copy to X (will be used for progress reporting)
         */
        state->xn.ptr.p_double[i] = boundval(state->xn.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
        
        /*
         * update active set
         */
        if( ae_fp_eq(state->xn.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xn.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->an.ptr.p_double[i] = 0;
        }
        else
        {
            state->an.ptr.p_double[i] = 1;
        }
        if( ae_fp_neq(state->an.ptr.p_double[i],state->ak.ptr.p_double[i]) )
        {
            diffcnt = diffcnt+1;
        }
        state->ak.ptr.p_double[i] = state->an.ptr.p_double[i];
    }
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->repnfev = state->repnfev+state->nfev;
    state->repiterationscount = state->repiterationscount+1;
    if( !state->xrep )
    {
        goto lbl_53;
    }
    
    /*
     * progress report
     */
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->xupdated = ae_false;
lbl_53:
    
    /*
     * Update info about step length
     */
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->laststep = ae_sqrt(v, _state)*state->stp;
    
    /*
     * Check stopping conditions.
     */
    if( ae_fp_greater(minasa_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        goto lbl_55;
    }
    
    /*
     * Gradient is small enough
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_57;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->xupdated = ae_false;
lbl_57:
    result = ae_false;
    return result;
lbl_55:
    if( !(state->repiterationscount>=state->maxits&&state->maxits>0) )
    {
        goto lbl_59;
    }
    
    /*
     * Too many iterations
     */
    state->repterminationtype = 5;
    if( !state->xrep )
    {
        goto lbl_61;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->xupdated = ae_false;
lbl_61:
    result = ae_false;
    return result;
lbl_59:
    if( !(ae_fp_greater_eq(minasa_asaginorm(state, _state),state->mu*minasa_asad1norm(state, _state))&&diffcnt==0) )
    {
        goto lbl_63;
    }
    
    /*
     * These conditions (EpsF/EpsX) are explicitly or implicitly
     * related to the current step size and influenced
     * by changes in the active constraints.
     *
     * For these reasons they are checked only when we don't
     * want to 'unstick' at the end of the iteration and there
     * were no changes in the active set.
     *
     * NOTE: consition |G|>=Mu*|D1| must be exactly opposite
     * to the condition used to switch back to GPA. At least
     * one inequality must be strict, otherwise infinite cycle
     * may occur when |G|=Mu*|D1| (we DON'T test stopping
     * conditions and we DON'T switch to GPA, so we cycle
     * indefinitely).
     */
    if( ae_fp_greater(state->fold-state->f,state->epsf*ae_maxreal(ae_fabs(state->fold, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
    {
        goto lbl_65;
    }
    
    /*
     * F(k+1)-F(k) is small enough
     */
    state->repterminationtype = 1;
    if( !state->xrep )
    {
        goto lbl_67;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->xupdated = ae_false;
lbl_67:
    result = ae_false;
    return result;
lbl_65:
    if( ae_fp_greater(state->laststep,state->epsx) )
    {
        goto lbl_69;
    }
    
    /*
     * X(k+1)-X(k) is small enough
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_71;
    }
    minasa_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->xupdated = ae_false;
lbl_71:
    result = ae_false;
    return result;
lbl_69:
lbl_63:
    
    /*
     * Check conditions for switching
     */
    if( ae_fp_less(minasa_asaginorm(state, _state),state->mu*minasa_asad1norm(state, _state)) )
    {
        state->curalgo = 0;
        goto lbl_50;
    }
    if( diffcnt>0 )
    {
        if( minasa_asauisempty(state, _state)||diffcnt>=minasa_n2 )
        {
            state->curalgo = 1;
        }
        else
        {
            state->curalgo = 0;
        }
        goto lbl_50;
    }
    
    /*
     * Calculate D(k+1)
     *
     * Line search may result in:
     * * maximum feasible step being taken (already processed)
     * * point satisfying Wolfe conditions
     * * some kind of error (CG is restarted by assigning 0.0 to Beta)
     */
    if( mcinfo==1 )
    {
        
        /*
         * Standard Wolfe conditions are satisfied:
         * * calculate Y[K] and BetaK
         */
        ae_v_add(&state->yk.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
        vv = ae_v_dotproduct(&state->yk.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = ae_v_dotproduct(&state->gc.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->betady = v/vv;
        v = ae_v_dotproduct(&state->gc.ptr.p_double[0], 1, &state->yk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->betahs = v/vv;
        if( state->cgtype==0 )
        {
            betak = state->betady;
        }
        if( state->cgtype==1 )
        {
            betak = ae_maxreal(0, ae_minreal(state->betady, state->betahs, _state), _state);
        }
    }
    else
    {
        
        /*
         * Something is wrong (may be function is too wild or too flat).
         *
         * We'll set BetaK=0, which will restart CG algorithm.
         * We can stop later (during normal checks) if stopping conditions are met.
         */
        betak = 0;
        state->debugrestartscount = state->debugrestartscount+1;
    }
    ae_v_moveneg(&state->dn.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->dn.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1), betak);
    ae_v_move(&state->dk.ptr.p_double[0], 1, &state->dn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * update other information
     */
    state->fold = state->f;
    state->k = state->k+1;
    goto lbl_49;
lbl_50:
lbl_47:
    goto lbl_17;
lbl_18:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = i;
    state->rstate.ia.ptr.p_int[2] = mcinfo;
    state->rstate.ia.ptr.p_int[3] = diffcnt;
    state->rstate.ba.ptr.p_bool[0] = b;
    state->rstate.ba.ptr.p_bool[1] = stepfound;
    state->rstate.ra.ptr.p_double[0] = betak;
    state->rstate.ra.ptr.p_double[1] = v;
    state->rstate.ra.ptr.p_double[2] = vv;
    return result;
}


/*************************************************************************
ASA results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations
                * ActiveConstraints contains number of active constraints

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresults(minasastate* state,
     /* Real    */ ae_vector* x,
     minasareport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minasareport_clear(rep);

    minasaresultsbuf(state, x, rep, _state);
}


/*************************************************************************
ASA results

Buffered implementation of MinASAResults() which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresultsbuf(minasastate* state,
     /* Real    */ ae_vector* x,
     minasareport* rep,
     ae_state *_state)
{
    ae_int_t i;


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->nfev = state->repnfev;
    rep->terminationtype = state->repterminationtype;
    rep->activeconstraints = 0;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_eq(state->ak.ptr.p_double[i],0) )
        {
            rep->activeconstraints = rep->activeconstraints+1;
        }
    }
}


/*************************************************************************
This  subroutine  restarts  CG  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure previously allocated with MinCGCreate call.
    X       -   new starting point.
    BndL    -   new lower bounds
    BndU    -   new upper bounds

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minasarestartfrom(minasastate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinASARestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinASARestartFrom: X contains infinite or NaN values!", _state);
    ae_assert(bndl->cnt>=state->n, "MinASARestartFrom: Length(BndL)<N!", _state);
    ae_assert(isfinitevector(bndl, state->n, _state), "MinASARestartFrom: BndL contains infinite or NaN values!", _state);
    ae_assert(bndu->cnt>=state->n, "MinASARestartFrom: Length(BndU)<N!", _state);
    ae_assert(isfinitevector(bndu, state->n, _state), "MinASARestartFrom: BndU contains infinite or NaN values!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_v_move(&state->bndl.ptr.p_double[0], 1, &bndl->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_v_move(&state->bndu.ptr.p_double[0], 1, &bndu->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    state->laststep = 0;
    ae_vector_set_length(&state->rstate.ia, 3+1, _state);
    ae_vector_set_length(&state->rstate.ba, 1+1, _state);
    ae_vector_set_length(&state->rstate.ra, 2+1, _state);
    state->rstate.stage = -1;
    minasa_clearrequestfields(state, _state);
}


/*************************************************************************
Returns norm of bounded anti-gradient.

Bounded antigradient is a vector obtained from  anti-gradient  by  zeroing
components which point outwards:
    result = norm(v)
    v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or
                  ((-g[i]>0)and(x[i]=bndu[i]))
    v[i]=-g[i] otherwise

This function may be used to check a stopping criterion.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double minasa_asaboundedantigradnorm(minasastate* state,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        v = -state->g.ptr.p_double[i];
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])&&ae_fp_less(-state->g.ptr.p_double[i],0) )
        {
            v = 0;
        }
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i])&&ae_fp_greater(-state->g.ptr.p_double[i],0) )
        {
            v = 0;
        }
        result = result+ae_sqr(v, _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns norm of GI(x).

GI(x) is  a  gradient  vector  whose  components  associated  with  active
constraints are zeroed. It  differs  from  bounded  anti-gradient  because
components  of   GI(x)   are   zeroed  independently  of  sign(g[i]),  and
anti-gradient's components are zeroed with respect to both constraint  and
sign.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double minasa_asaginorm(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_neq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])&&ae_fp_neq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            result = result+ae_sqr(state->g.ptr.p_double[i], _state);
        }
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns norm(D1(State.X))

For a meaning of D1 see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double minasa_asad1norm(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        result = result+ae_sqr(boundval(state->x.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state)-state->x.ptr.p_double[i], _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns True, if U set is empty.

* State.X is used as point,
* State.G - as gradient,
* D is calculated within function (because State.D may have different
  meaning depending on current optimization algorithm)

For a meaning of U see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static ae_bool minasa_asauisempty(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double d;
    double d2;
    double d32;
    ae_bool result;


    d = minasa_asad1norm(state, _state);
    d2 = ae_sqrt(d, _state);
    d32 = d*d2;
    result = ae_true;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_greater_eq(ae_fabs(state->g.ptr.p_double[i], _state),d2)&&ae_fp_greater_eq(ae_minreal(state->x.ptr.p_double[i]-state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i]-state->x.ptr.p_double[i], _state),d32) )
        {
            result = ae_false;
            return result;
        }
    }
    return result;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void minasa_clearrequestfields(minasastate* state,
     ae_state *_state)
{


    state->needfg = ae_false;
    state->xupdated = ae_false;
}


ae_bool _minasastate_init(minasastate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->bndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ak, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->an, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->d, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->work, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->yk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !_linminstate_init(&p->lstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minasastate_init_copy(minasastate* dst, minasastate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->cgtype = src->cgtype;
    dst->k = src->k;
    dst->nfev = src->nfev;
    dst->mcstage = src->mcstage;
    if( !ae_vector_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    dst->curalgo = src->curalgo;
    dst->acount = src->acount;
    dst->mu = src->mu;
    dst->finit = src->finit;
    dst->dginit = src->dginit;
    if( !ae_vector_init_copy(&dst->ak, &src->ak, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xk, &src->xk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dk, &src->dk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->an, &src->an, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xn, &src->xn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dn, &src->dn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->d, &src->d, _state, make_automatic) )
        return ae_false;
    dst->fold = src->fold;
    dst->stp = src->stp;
    if( !ae_vector_init_copy(&dst->work, &src->work, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->yk, &src->yk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gc, &src->gc, _state, make_automatic) )
        return ae_false;
    dst->laststep = src->laststep;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repiterationscount = src->repiterationscount;
    dst->repnfev = src->repnfev;
    dst->repterminationtype = src->repterminationtype;
    dst->debugrestartscount = src->debugrestartscount;
    if( !_linminstate_init_copy(&dst->lstate, &src->lstate, _state, make_automatic) )
        return ae_false;
    dst->betahs = src->betahs;
    dst->betady = src->betady;
    return ae_true;
}


void _minasastate_clear(minasastate* p)
{
    ae_vector_clear(&p->bndl);
    ae_vector_clear(&p->bndu);
    ae_vector_clear(&p->ak);
    ae_vector_clear(&p->xk);
    ae_vector_clear(&p->dk);
    ae_vector_clear(&p->an);
    ae_vector_clear(&p->xn);
    ae_vector_clear(&p->dn);
    ae_vector_clear(&p->d);
    ae_vector_clear(&p->work);
    ae_vector_clear(&p->yk);
    ae_vector_clear(&p->gc);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    _linminstate_clear(&p->lstate);
}


ae_bool _minasareport_init(minasareport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _minasareport_init_copy(minasareport* dst, minasareport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->activeconstraints = src->activeconstraints;
    return ae_true;
}


void _minasareport_clear(minasareport* p)
{
}


/*$ End $*/
