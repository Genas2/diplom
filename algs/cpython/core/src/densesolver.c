/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
#include "densesolver.h"


/*$ Declarations $*/
static void densesolver_rmatrixlusolveinternal(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_bool havea,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
static void densesolver_spdmatrixcholeskysolveinternal(/* Real    */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* a,
     ae_bool havea,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
static void densesolver_cmatrixlusolveinternal(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_bool havea,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
static void densesolver_hpdmatrixcholeskysolveinternal(/* Complex */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* a,
     ae_bool havea,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
static ae_int_t densesolver_densesolverrfsmax(ae_int_t n,
     double r1,
     double rinf,
     ae_state *_state);
static ae_int_t densesolver_densesolverrfsmaxv2(ae_int_t n,
     double r2,
     ae_state *_state);
static void densesolver_rbasiclusolve(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Real    */ ae_vector* xb,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static void densesolver_spdbasiccholeskysolve(/* Real    */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* xb,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static void densesolver_cbasiclusolve(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Complex */ ae_vector* xb,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static void densesolver_hpdbasiccholeskysolve(/* Complex */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* xb,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix, x and b are vectors.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1], it contains:
                * solution of A*x=b if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R1        reciprocal of condition number: 1/cond(A), 1-norm.
* RInf      reciprocal of condition number: 1/cond(A), inf-norm.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_REAL, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_move(&bm.ptr.pp_double[0][0], bm.stride, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    rmatrixsolvem(a, n, &bm, 1, ae_true, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_move(&x->ptr.p_double[0], 1, &xm.ptr.pp_double[0][0], xm.stride, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* optional iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix da;
    ae_matrix emptya;
    ae_vector p;
    double scalea;
    ae_int_t i;
    ae_int_t j;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&da, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&emptya, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&p, 0, DT_INT, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&da, n, n, _state);
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&da.ptr.pp_double[i][0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
    }
    rmatrixlu(&da, n, n, &p, _state);
    if( rfs )
    {
        densesolver_rmatrixlusolveinternal(&da, &p, scalea, n, a, ae_true, b, m, info, rep, x, _state);
    }
    else
    {
        densesolver_rmatrixlusolveinternal(&da, &p, scalea, n, &emptya, ae_false, b, m, info, rep, x, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*X=B,  where A is NxN non-denegerate
real matrix given by its LU decomposition, X and B are NxM real matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve
    
  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolve(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_REAL, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_move(&bm.ptr.pp_double[0][0], bm.stride, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    rmatrixlusolvem(lua, p, n, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_move(&x->ptr.p_double[0], 1, &xm.ptr.pp_double[0][0], xm.stride, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves task with multiple right parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvem(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix emptya;
    ae_int_t i;
    ae_int_t j;
    double scalea;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&emptya, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * 1. scale matrix, max(|U[i,j]|)
     *    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
     * 2. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=i; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_fabs(lua->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    densesolver_rmatrixlusolveinternal(lua, p, scalea, n, &emptya, ae_false, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
both A and its LU decomposition.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolve(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_REAL, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_move(&bm.ptr.pp_double[0][0], bm.stride, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    rmatrixmixedsolvem(a, lua, p, n, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_move(&x->ptr.p_double[0], 1, &xm.ptr.pp_double[0][0], xm.stride, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolvem(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    double scalea;
    ae_int_t i;
    ae_int_t j;

    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        return;
    }
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    densesolver_rmatrixlusolveinternal(lua, p, scalea, n, a, ae_true, b, m, info, rep, x, _state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix da;
    ae_matrix emptya;
    ae_vector p;
    double scalea;
    ae_int_t i;
    ae_int_t j;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&da, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&emptya, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&p, 0, DT_INT, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&da, n, n, _state);
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_c_abs(a->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    for(i=0; i<=n-1; i++)
    {
        ae_v_cmove(&da.ptr.pp_complex[i][0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,n-1));
    }
    cmatrixlu(&da, n, n, &p, _state);
    if( rfs )
    {
        densesolver_cmatrixlusolveinternal(&da, &p, scalea, n, a, ae_true, b, m, info, rep, x, _state);
    }
    else
    {
        densesolver_cmatrixlusolveinternal(&da, &p, scalea, n, &emptya, ae_false, b, m, info, rep, x, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_COMPLEX, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_cmove(&bm.ptr.pp_complex[0][0], bm.stride, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
    cmatrixsolvem(a, n, &bm, 1, ae_true, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_cmove(&x->ptr.p_complex[0], 1, &xm.ptr.pp_complex[0][0], xm.stride, "N", ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvem(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix emptya;
    ae_int_t i;
    ae_int_t j;
    double scalea;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&emptya, 0, 0, DT_COMPLEX, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * 1. scale matrix, max(|U[i,j]|)
     *    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
     * 2. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=i; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_c_abs(lua->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    densesolver_cmatrixlusolveinternal(lua, p, scalea, n, &emptya, ae_false, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolve(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_COMPLEX, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_cmove(&bm.ptr.pp_complex[0][0], bm.stride, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
    cmatrixlusolvem(lua, p, n, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_cmove(&x->ptr.p_complex[0], 1, &xm.ptr.pp_complex[0][0], xm.stride, "N", ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolvem(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    double scalea;
    ae_int_t i;
    ae_int_t j;

    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        return;
    }
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    scalea = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            scalea = ae_maxreal(scalea, ae_c_abs(a->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    densesolver_cmatrixlusolveinternal(lua, p, scalea, n, a, ae_true, b, m, info, rep, x, _state);
}


/*************************************************************************
Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolve(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_COMPLEX, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_cmove(&bm.ptr.pp_complex[0][0], bm.stride, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
    cmatrixmixedsolvem(a, lua, p, n, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_cmove(&x->ptr.p_complex[0], 1, &xm.ptr.pp_complex[0][0], xm.stride, "N", ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for symmetric positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix da;
    double sqrtscalea;
    ae_int_t i;
    ae_int_t j;
    ae_int_t j1;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&da, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&da, n, n, _state);
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    sqrtscalea = 0;
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
            sqrtscalea = ae_maxreal(sqrtscalea, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    sqrtscalea = ae_sqrt(sqrtscalea, _state);
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
        ae_v_move(&da.ptr.pp_double[i][j1], 1, &a->ptr.pp_double[i][j1], 1, ae_v_len(j1,j2));
    }
    if( !spdmatrixcholesky(&da, n, isupper, _state) )
    {
        ae_matrix_set_length(x, n, m, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_double[i][j] = 0;
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    densesolver_spdmatrixcholeskysolveinternal(&da, sqrtscalea, n, isupper, a, ae_true, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(), but for SPD matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_REAL, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_move(&bm.ptr.pp_double[0][0], bm.stride, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    spdmatrixsolvem(a, n, isupper, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_move(&x->ptr.p_double[0], 1, &xm.ptr.pp_double[0][0], xm.stride, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix emptya;
    double sqrtscalea;
    ae_int_t i;
    ae_int_t j;
    ae_int_t j1;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&emptya, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * 1. scale matrix, max(|U[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    sqrtscalea = 0;
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
            sqrtscalea = ae_maxreal(sqrtscalea, ae_fabs(cha->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    densesolver_spdmatrixcholeskysolveinternal(cha, sqrtscalea, n, isupper, &emptya, ae_false, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolve(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_REAL, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_move(&bm.ptr.pp_double[0][0], bm.stride, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    spdmatrixcholeskysolvem(cha, n, isupper, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_move(&x->ptr.p_double[0], 1, &xm.ptr.pp_double[0][0], xm.stride, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix da;
    double sqrtscalea;
    ae_int_t i;
    ae_int_t j;
    ae_int_t j1;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&da, 0, 0, DT_COMPLEX, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&da, n, n, _state);
    
    /*
     * 1. scale matrix, max(|A[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    sqrtscalea = 0;
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
            sqrtscalea = ae_maxreal(sqrtscalea, ae_c_abs(a->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    sqrtscalea = ae_sqrt(sqrtscalea, _state);
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
        ae_v_cmove(&da.ptr.pp_complex[i][j1], 1, &a->ptr.pp_complex[i][j1], 1, "N", ae_v_len(j1,j2));
    }
    if( !hpdmatrixcholesky(&da, n, isupper, _state) )
    {
        ae_matrix_set_length(x, n, m, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_complex[i][j] = ae_complex_from_d(0);
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    densesolver_hpdmatrixcholeskysolveinternal(&da, sqrtscalea, n, isupper, a, ae_true, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(),  but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_COMPLEX, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_cmove(&bm.ptr.pp_complex[0][0], bm.stride, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
    hpdmatrixsolvem(a, n, isupper, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_cmove(&x->ptr.p_complex[0], 1, &xm.ptr.pp_complex[0][0], xm.stride, "N", ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix emptya;
    double sqrtscalea;
    ae_int_t i;
    ae_int_t j;
    ae_int_t j1;
    ae_int_t j2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_matrix_init(&emptya, 0, 0, DT_COMPLEX, _state, ae_true);

    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * 1. scale matrix, max(|U[i,j]|)
     * 2. factorize scaled matrix
     * 3. solve
     */
    sqrtscalea = 0;
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
            sqrtscalea = ae_maxreal(sqrtscalea, ae_c_abs(cha->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    densesolver_hpdmatrixcholeskysolveinternal(cha, sqrtscalea, n, isupper, &emptya, ae_false, b, m, info, rep, x, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolve(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix bm;
    ae_matrix xm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_vector_clear(x);
    ae_matrix_init(&bm, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&xm, 0, 0, DT_COMPLEX, _state, ae_true);

    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(&bm, n, 1, _state);
    ae_v_cmove(&bm.ptr.pp_complex[0][0], bm.stride, &b->ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
    hpdmatrixcholeskysolvem(cha, n, isupper, &bm, 1, info, rep, &xm, _state);
    ae_vector_set_length(x, n, _state);
    ae_v_cmove(&x->ptr.p_complex[0], 1, &xm.ptr.pp_complex[0][0], xm.stride, "N", ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the  usual sense is returned

Algorithm features:
* automatic detection of degenerate cases
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.
                
OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvels(/* Real    */ ae_matrix* a,
     ae_int_t nrows,
     ae_int_t ncols,
     /* Real    */ ae_vector* b,
     double threshold,
     ae_int_t* info,
     densesolverlsreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector sv;
    ae_matrix u;
    ae_matrix vt;
    ae_vector rp;
    ae_vector utb;
    ae_vector sutb;
    ae_vector tmp;
    ae_vector ta;
    ae_vector tx;
    ae_vector buf;
    ae_vector w;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nsv;
    ae_int_t kernelidx;
    double v;
    double verr;
    ae_bool svdfailed;
    ae_bool zeroa;
    ae_int_t rfs;
    ae_int_t nrfs;
    ae_bool terminatenexttime;
    ae_bool smallerr;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverlsreport_clear(rep);
    ae_vector_clear(x);
    ae_vector_init(&sv, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&rp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&utb, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sutb, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ta, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);

    if( (nrows<=0||ncols<=0)||ae_fp_less(threshold,0) )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( ae_fp_eq(threshold,0) )
    {
        threshold = 1000*ae_machineepsilon;
    }
    
    /*
     * Factorize A first
     */
    svdfailed = !rmatrixsvd(a, nrows, ncols, 1, 2, 2, &sv, &u, &vt, _state);
    zeroa = ae_fp_eq(sv.ptr.p_double[0],0);
    if( svdfailed||zeroa )
    {
        if( svdfailed )
        {
            *info = -4;
        }
        else
        {
            *info = 1;
        }
        ae_vector_set_length(x, ncols, _state);
        for(i=0; i<=ncols-1; i++)
        {
            x->ptr.p_double[i] = 0;
        }
        rep->n = ncols;
        rep->k = ncols;
        ae_matrix_set_length(&rep->cx, ncols, ncols, _state);
        for(i=0; i<=ncols-1; i++)
        {
            for(j=0; j<=ncols-1; j++)
            {
                if( i==j )
                {
                    rep->cx.ptr.pp_double[i][j] = 1;
                }
                else
                {
                    rep->cx.ptr.pp_double[i][j] = 0;
                }
            }
        }
        rep->r2 = 0;
        ae_frame_leave(_state);
        return;
    }
    nsv = ae_minint(ncols, nrows, _state);
    if( nsv==ncols )
    {
        rep->r2 = sv.ptr.p_double[nsv-1]/sv.ptr.p_double[0];
    }
    else
    {
        rep->r2 = 0;
    }
    rep->n = ncols;
    *info = 1;
    
    /*
     * Iterative refinement of xc combined with solution:
     * 1. xc = 0
     * 2. calculate r = bc-A*xc using extra-precise dot product
     * 3. solve A*y = r
     * 4. update x:=x+r
     * 5. goto 2
     *
     * This cycle is executed until one of two things happens:
     * 1. maximum number of iterations reached
     * 2. last iteration decreased error to the lower limit
     */
    ae_vector_set_length(&utb, nsv, _state);
    ae_vector_set_length(&sutb, nsv, _state);
    ae_vector_set_length(x, ncols, _state);
    ae_vector_set_length(&tmp, ncols, _state);
    ae_vector_set_length(&ta, ncols+1, _state);
    ae_vector_set_length(&tx, ncols+1, _state);
    ae_vector_set_length(&buf, ncols+1, _state);
    for(i=0; i<=ncols-1; i++)
    {
        x->ptr.p_double[i] = 0;
    }
    kernelidx = nsv;
    for(i=0; i<=nsv-1; i++)
    {
        if( ae_fp_less_eq(sv.ptr.p_double[i],threshold*sv.ptr.p_double[0]) )
        {
            kernelidx = i;
            break;
        }
    }
    rep->k = ncols-kernelidx;
    nrfs = densesolver_densesolverrfsmaxv2(ncols, rep->r2, _state);
    terminatenexttime = ae_false;
    ae_vector_set_length(&rp, nrows, _state);
    for(rfs=0; rfs<=nrfs; rfs++)
    {
        if( terminatenexttime )
        {
            break;
        }
        
        /*
         * calculate right part
         */
        if( rfs==0 )
        {
            ae_v_move(&rp.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,nrows-1));
        }
        else
        {
            smallerr = ae_true;
            for(i=0; i<=nrows-1; i++)
            {
                ae_v_move(&ta.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,ncols-1));
                ta.ptr.p_double[ncols] = -1;
                ae_v_move(&tx.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,ncols-1));
                tx.ptr.p_double[ncols] = b->ptr.p_double[i];
                xdot(&ta, &tx, ncols+1, &buf, &v, &verr, _state);
                rp.ptr.p_double[i] = -v;
                smallerr = smallerr&&ae_fp_less(ae_fabs(v, _state),4*verr);
            }
            if( smallerr )
            {
                terminatenexttime = ae_true;
            }
        }
        
        /*
         * solve A*dx = rp
         */
        for(i=0; i<=ncols-1; i++)
        {
            tmp.ptr.p_double[i] = 0;
        }
        for(i=0; i<=nsv-1; i++)
        {
            utb.ptr.p_double[i] = 0;
        }
        for(i=0; i<=nrows-1; i++)
        {
            v = rp.ptr.p_double[i];
            ae_v_addd(&utb.ptr.p_double[0], 1, &u.ptr.pp_double[i][0], 1, ae_v_len(0,nsv-1), v);
        }
        for(i=0; i<=nsv-1; i++)
        {
            if( i<kernelidx )
            {
                sutb.ptr.p_double[i] = utb.ptr.p_double[i]/sv.ptr.p_double[i];
            }
            else
            {
                sutb.ptr.p_double[i] = 0;
            }
        }
        for(i=0; i<=nsv-1; i++)
        {
            v = sutb.ptr.p_double[i];
            ae_v_addd(&tmp.ptr.p_double[0], 1, &vt.ptr.pp_double[i][0], 1, ae_v_len(0,ncols-1), v);
        }
        
        /*
         * update x:  x:=x+dx
         */
        ae_v_add(&x->ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,ncols-1));
    }
    
    /*
     * fill CX
     */
    if( rep->k>0 )
    {
        ae_matrix_set_length(&rep->cx, ncols, rep->k, _state);
        for(i=0; i<=rep->k-1; i++)
        {
            ae_v_move(&rep->cx.ptr.pp_double[0][i], rep->cx.stride, &vt.ptr.pp_double[kernelidx+i][0], 1, ae_v_len(0,ncols-1));
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_rmatrixlusolveinternal(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_bool havea,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t rfs;
    ae_int_t nrfs;
    ae_vector xc;
    ae_vector y;
    ae_vector bc;
    ae_vector xa;
    ae_vector xb;
    ae_vector tx;
    double v;
    double verr;
    double mxb;
    double scaleright;
    ae_bool smallerr;
    ae_bool terminatenexttime;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xa, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xb, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tx, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_fp_greater(scalea,0), "Assertion failed", _state);
    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( p->ptr.p_int[i]>n-1||p->ptr.p_int[i]<i )
        {
            *info = -1;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(x, n, m, _state);
    ae_vector_set_length(&y, n, _state);
    ae_vector_set_length(&xc, n, _state);
    ae_vector_set_length(&bc, n, _state);
    ae_vector_set_length(&tx, n+1, _state);
    ae_vector_set_length(&xa, n+1, _state);
    ae_vector_set_length(&xb, n+1, _state);
    
    /*
     * estimate condition number, test for near singularity
     */
    rep->r1 = rmatrixlurcond1(lua, n, _state);
    rep->rinf = rmatrixlurcondinf(lua, n, _state);
    if( ae_fp_less(rep->r1,rcondthreshold(_state))||ae_fp_less(rep->rinf,rcondthreshold(_state)) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_double[i][j] = 0;
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * solve
     */
    for(k=0; k<=m-1; k++)
    {
        
        /*
         * copy B to contiguous storage
         */
        ae_v_move(&bc.ptr.p_double[0], 1, &b->ptr.pp_double[0][k], b->stride, ae_v_len(0,n-1));
        
        /*
         * Scale right part:
         * * MX stores max(|Bi|)
         * * ScaleRight stores actual scaling applied to B when solving systems
         *   it is chosen to make |scaleRight*b| close to 1.
         */
        mxb = 0;
        for(i=0; i<=n-1; i++)
        {
            mxb = ae_maxreal(mxb, ae_fabs(bc.ptr.p_double[i], _state), _state);
        }
        if( ae_fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        /*
         * First, non-iterative part of solution process.
         * We use separate code for this task because
         * XDot is quite slow and we want to save time.
         */
        ae_v_moved(&xc.ptr.p_double[0], 1, &bc.ptr.p_double[0], 1, ae_v_len(0,n-1), scaleright);
        densesolver_rbasiclusolve(lua, p, scalea, n, &xc, &tx, _state);
        
        /*
         * Iterative refinement of xc:
         * * calculate r = bc-A*xc using extra-precise dot product
         * * solve A*y = r
         * * update x:=x+r
         *
         * This cycle is executed until one of two things happens:
         * 1. maximum number of iterations reached
         * 2. last iteration decreased error to the lower limit
         */
        if( havea )
        {
            nrfs = densesolver_densesolverrfsmax(n, rep->r1, rep->rinf, _state);
            terminatenexttime = ae_false;
            for(rfs=0; rfs<=nrfs-1; rfs++)
            {
                if( terminatenexttime )
                {
                    break;
                }
                
                /*
                 * generate right part
                 */
                smallerr = ae_true;
                ae_v_move(&xb.ptr.p_double[0], 1, &xc.ptr.p_double[0], 1, ae_v_len(0,n-1));
                for(i=0; i<=n-1; i++)
                {
                    ae_v_moved(&xa.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), scalea);
                    xa.ptr.p_double[n] = -1;
                    xb.ptr.p_double[n] = scaleright*bc.ptr.p_double[i];
                    xdot(&xa, &xb, n+1, &tx, &v, &verr, _state);
                    y.ptr.p_double[i] = -v;
                    smallerr = smallerr&&ae_fp_less(ae_fabs(v, _state),4*verr);
                }
                if( smallerr )
                {
                    terminatenexttime = ae_true;
                }
                
                /*
                 * solve and update
                 */
                densesolver_rbasiclusolve(lua, p, scalea, n, &y, &tx, _state);
                ae_v_add(&xc.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,n-1));
            }
        }
        
        /*
         * Store xc.
         * Post-scale result.
         */
        v = scalea*mxb;
        ae_v_moved(&x->ptr.pp_double[0][k], x->stride, &xc.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_spdmatrixcholeskysolveinternal(/* Real    */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* a,
     ae_bool havea,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_vector xc;
    ae_vector y;
    ae_vector bc;
    ae_vector xa;
    ae_vector xb;
    ae_vector tx;
    double v;
    double mxb;
    double scaleright;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xa, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xb, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tx, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_fp_greater(sqrtscalea,0), "Assertion failed", _state);
    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(x, n, m, _state);
    ae_vector_set_length(&y, n, _state);
    ae_vector_set_length(&xc, n, _state);
    ae_vector_set_length(&bc, n, _state);
    ae_vector_set_length(&tx, n+1, _state);
    ae_vector_set_length(&xa, n+1, _state);
    ae_vector_set_length(&xb, n+1, _state);
    
    /*
     * estimate condition number, test for near singularity
     */
    rep->r1 = spdmatrixcholeskyrcond(cha, n, isupper, _state);
    rep->rinf = rep->r1;
    if( ae_fp_less(rep->r1,rcondthreshold(_state)) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_double[i][j] = 0;
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * solve
     */
    for(k=0; k<=m-1; k++)
    {
        
        /*
         * copy B to contiguous storage
         */
        ae_v_move(&bc.ptr.p_double[0], 1, &b->ptr.pp_double[0][k], b->stride, ae_v_len(0,n-1));
        
        /*
         * Scale right part:
         * * MX stores max(|Bi|)
         * * ScaleRight stores actual scaling applied to B when solving systems
         *   it is chosen to make |scaleRight*b| close to 1.
         */
        mxb = 0;
        for(i=0; i<=n-1; i++)
        {
            mxb = ae_maxreal(mxb, ae_fabs(bc.ptr.p_double[i], _state), _state);
        }
        if( ae_fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        /*
         * First, non-iterative part of solution process.
         * We use separate code for this task because
         * XDot is quite slow and we want to save time.
         */
        ae_v_moved(&xc.ptr.p_double[0], 1, &bc.ptr.p_double[0], 1, ae_v_len(0,n-1), scaleright);
        densesolver_spdbasiccholeskysolve(cha, sqrtscalea, n, isupper, &xc, &tx, _state);
        
        /*
         * Store xc.
         * Post-scale result.
         */
        v = ae_sqr(sqrtscalea, _state)*mxb;
        ae_v_moved(&x->ptr.pp_double[0][k], x->stride, &xc.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_cmatrixlusolveinternal(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_bool havea,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t rfs;
    ae_int_t nrfs;
    ae_vector xc;
    ae_vector y;
    ae_vector bc;
    ae_vector xa;
    ae_vector xb;
    ae_vector tx;
    ae_vector tmpbuf;
    ae_complex v;
    double verr;
    double mxb;
    double scaleright;
    ae_bool smallerr;
    ae_bool terminatenexttime;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_vector_init(&xc, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&y, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&bc, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&xa, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&xb, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&tx, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&tmpbuf, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_fp_greater(scalea,0), "Assertion failed", _state);
    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( p->ptr.p_int[i]>n-1||p->ptr.p_int[i]<i )
        {
            *info = -1;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(x, n, m, _state);
    ae_vector_set_length(&y, n, _state);
    ae_vector_set_length(&xc, n, _state);
    ae_vector_set_length(&bc, n, _state);
    ae_vector_set_length(&tx, n, _state);
    ae_vector_set_length(&xa, n+1, _state);
    ae_vector_set_length(&xb, n+1, _state);
    ae_vector_set_length(&tmpbuf, 2*n+2, _state);
    
    /*
     * estimate condition number, test for near singularity
     */
    rep->r1 = cmatrixlurcond1(lua, n, _state);
    rep->rinf = cmatrixlurcondinf(lua, n, _state);
    if( ae_fp_less(rep->r1,rcondthreshold(_state))||ae_fp_less(rep->rinf,rcondthreshold(_state)) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_complex[i][j] = ae_complex_from_d(0);
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * solve
     */
    for(k=0; k<=m-1; k++)
    {
        
        /*
         * copy B to contiguous storage
         */
        ae_v_cmove(&bc.ptr.p_complex[0], 1, &b->ptr.pp_complex[0][k], b->stride, "N", ae_v_len(0,n-1));
        
        /*
         * Scale right part:
         * * MX stores max(|Bi|)
         * * ScaleRight stores actual scaling applied to B when solving systems
         *   it is chosen to make |scaleRight*b| close to 1.
         */
        mxb = 0;
        for(i=0; i<=n-1; i++)
        {
            mxb = ae_maxreal(mxb, ae_c_abs(bc.ptr.p_complex[i], _state), _state);
        }
        if( ae_fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        /*
         * First, non-iterative part of solution process.
         * We use separate code for this task because
         * XDot is quite slow and we want to save time.
         */
        ae_v_cmoved(&xc.ptr.p_complex[0], 1, &bc.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1), scaleright);
        densesolver_cbasiclusolve(lua, p, scalea, n, &xc, &tx, _state);
        
        /*
         * Iterative refinement of xc:
         * * calculate r = bc-A*xc using extra-precise dot product
         * * solve A*y = r
         * * update x:=x+r
         *
         * This cycle is executed until one of two things happens:
         * 1. maximum number of iterations reached
         * 2. last iteration decreased error to the lower limit
         */
        if( havea )
        {
            nrfs = densesolver_densesolverrfsmax(n, rep->r1, rep->rinf, _state);
            terminatenexttime = ae_false;
            for(rfs=0; rfs<=nrfs-1; rfs++)
            {
                if( terminatenexttime )
                {
                    break;
                }
                
                /*
                 * generate right part
                 */
                smallerr = ae_true;
                ae_v_cmove(&xb.ptr.p_complex[0], 1, &xc.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
                for(i=0; i<=n-1; i++)
                {
                    ae_v_cmoved(&xa.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,n-1), scalea);
                    xa.ptr.p_complex[n] = ae_complex_from_d(-1);
                    xb.ptr.p_complex[n] = ae_c_mul_d(bc.ptr.p_complex[i],scaleright);
                    xcdot(&xa, &xb, n+1, &tmpbuf, &v, &verr, _state);
                    y.ptr.p_complex[i] = ae_c_neg(v);
                    smallerr = smallerr&&ae_fp_less(ae_c_abs(v, _state),4*verr);
                }
                if( smallerr )
                {
                    terminatenexttime = ae_true;
                }
                
                /*
                 * solve and update
                 */
                densesolver_cbasiclusolve(lua, p, scalea, n, &y, &tx, _state);
                ae_v_cadd(&xc.ptr.p_complex[0], 1, &y.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
            }
        }
        
        /*
         * Store xc.
         * Post-scale result.
         */
        v = ae_complex_from_d(scalea*mxb);
        ae_v_cmovec(&x->ptr.pp_complex[0][k], x->stride, &xc.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_hpdmatrixcholeskysolveinternal(/* Complex */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* a,
     ae_bool havea,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_vector xc;
    ae_vector y;
    ae_vector bc;
    ae_vector xa;
    ae_vector xb;
    ae_vector tx;
    double v;
    double mxb;
    double scaleright;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _densesolverreport_clear(rep);
    ae_matrix_clear(x);
    ae_vector_init(&xc, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&y, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&bc, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&xa, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&xb, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&tx, 0, DT_COMPLEX, _state, ae_true);

    ae_assert(ae_fp_greater(sqrtscalea,0), "Assertion failed", _state);
    
    /*
     * prepare: check inputs, allocate space...
     */
    if( n<=0||m<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(x, n, m, _state);
    ae_vector_set_length(&y, n, _state);
    ae_vector_set_length(&xc, n, _state);
    ae_vector_set_length(&bc, n, _state);
    ae_vector_set_length(&tx, n+1, _state);
    ae_vector_set_length(&xa, n+1, _state);
    ae_vector_set_length(&xb, n+1, _state);
    
    /*
     * estimate condition number, test for near singularity
     */
    rep->r1 = hpdmatrixcholeskyrcond(cha, n, isupper, _state);
    rep->rinf = rep->r1;
    if( ae_fp_less(rep->r1,rcondthreshold(_state)) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                x->ptr.pp_complex[i][j] = ae_complex_from_d(0);
            }
        }
        rep->r1 = 0;
        rep->rinf = 0;
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * solve
     */
    for(k=0; k<=m-1; k++)
    {
        
        /*
         * copy B to contiguous storage
         */
        ae_v_cmove(&bc.ptr.p_complex[0], 1, &b->ptr.pp_complex[0][k], b->stride, "N", ae_v_len(0,n-1));
        
        /*
         * Scale right part:
         * * MX stores max(|Bi|)
         * * ScaleRight stores actual scaling applied to B when solving systems
         *   it is chosen to make |scaleRight*b| close to 1.
         */
        mxb = 0;
        for(i=0; i<=n-1; i++)
        {
            mxb = ae_maxreal(mxb, ae_c_abs(bc.ptr.p_complex[i], _state), _state);
        }
        if( ae_fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        /*
         * First, non-iterative part of solution process.
         * We use separate code for this task because
         * XDot is quite slow and we want to save time.
         */
        ae_v_cmoved(&xc.ptr.p_complex[0], 1, &bc.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1), scaleright);
        densesolver_hpdbasiccholeskysolve(cha, sqrtscalea, n, isupper, &xc, &tx, _state);
        
        /*
         * Store xc.
         * Post-scale result.
         */
        v = ae_sqr(sqrtscalea, _state)*mxb;
        ae_v_cmoved(&x->ptr.pp_complex[0][k], x->stride, &xc.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static ae_int_t densesolver_densesolverrfsmax(ae_int_t n,
     double r1,
     double rinf,
     ae_state *_state)
{
    ae_int_t result;


    result = 5;
    return result;
}


/*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. norm-2 condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static ae_int_t densesolver_densesolverrfsmaxv2(ae_int_t n,
     double r2,
     ae_state *_state)
{
    ae_int_t result;


    result = densesolver_densesolverrfsmax(n, 0, 0, _state);
    return result;
}


/*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_rbasiclusolve(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Real    */ ae_vector* xb,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    double v;


    for(i=0; i<=n-1; i++)
    {
        if( p->ptr.p_int[i]!=i )
        {
            v = xb->ptr.p_double[i];
            xb->ptr.p_double[i] = xb->ptr.p_double[p->ptr.p_int[i]];
            xb->ptr.p_double[p->ptr.p_int[i]] = v;
        }
    }
    for(i=1; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&lua->ptr.pp_double[i][0], 1, &xb->ptr.p_double[0], 1, ae_v_len(0,i-1));
        xb->ptr.p_double[i] = xb->ptr.p_double[i]-v;
    }
    xb->ptr.p_double[n-1] = xb->ptr.p_double[n-1]/(scalea*lua->ptr.pp_double[n-1][n-1]);
    for(i=n-2; i>=0; i--)
    {
        ae_v_moved(&tmp->ptr.p_double[i+1], 1, &lua->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), scalea);
        v = ae_v_dotproduct(&tmp->ptr.p_double[i+1], 1, &xb->ptr.p_double[i+1], 1, ae_v_len(i+1,n-1));
        xb->ptr.p_double[i] = (xb->ptr.p_double[i]-v)/(scalea*lua->ptr.pp_double[i][i]);
    }
}


/*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_spdbasiccholeskysolve(/* Real    */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* xb,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    double v;


    
    /*
     * A = L*L' or A=U'*U
     */
    if( isupper )
    {
        
        /*
         * Solve U'*y=b first.
         */
        for(i=0; i<=n-1; i++)
        {
            xb->ptr.p_double[i] = xb->ptr.p_double[i]/(sqrtscalea*cha->ptr.pp_double[i][i]);
            if( i<n-1 )
            {
                v = xb->ptr.p_double[i];
                ae_v_moved(&tmp->ptr.p_double[i+1], 1, &cha->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sqrtscalea);
                ae_v_subd(&xb->ptr.p_double[i+1], 1, &tmp->ptr.p_double[i+1], 1, ae_v_len(i+1,n-1), v);
            }
        }
        
        /*
         * Solve U*x=y then.
         */
        for(i=n-1; i>=0; i--)
        {
            if( i<n-1 )
            {
                ae_v_moved(&tmp->ptr.p_double[i+1], 1, &cha->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sqrtscalea);
                v = ae_v_dotproduct(&tmp->ptr.p_double[i+1], 1, &xb->ptr.p_double[i+1], 1, ae_v_len(i+1,n-1));
                xb->ptr.p_double[i] = xb->ptr.p_double[i]-v;
            }
            xb->ptr.p_double[i] = xb->ptr.p_double[i]/(sqrtscalea*cha->ptr.pp_double[i][i]);
        }
    }
    else
    {
        
        /*
         * Solve L*y=b first
         */
        for(i=0; i<=n-1; i++)
        {
            if( i>0 )
            {
                ae_v_moved(&tmp->ptr.p_double[0], 1, &cha->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sqrtscalea);
                v = ae_v_dotproduct(&tmp->ptr.p_double[0], 1, &xb->ptr.p_double[0], 1, ae_v_len(0,i-1));
                xb->ptr.p_double[i] = xb->ptr.p_double[i]-v;
            }
            xb->ptr.p_double[i] = xb->ptr.p_double[i]/(sqrtscalea*cha->ptr.pp_double[i][i]);
        }
        
        /*
         * Solve L'*x=y then.
         */
        for(i=n-1; i>=0; i--)
        {
            xb->ptr.p_double[i] = xb->ptr.p_double[i]/(sqrtscalea*cha->ptr.pp_double[i][i]);
            if( i>0 )
            {
                v = xb->ptr.p_double[i];
                ae_v_moved(&tmp->ptr.p_double[0], 1, &cha->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sqrtscalea);
                ae_v_subd(&xb->ptr.p_double[0], 1, &tmp->ptr.p_double[0], 1, ae_v_len(0,i-1), v);
            }
        }
    }
}


/*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_cbasiclusolve(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     double scalea,
     ae_int_t n,
     /* Complex */ ae_vector* xb,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_complex v;


    for(i=0; i<=n-1; i++)
    {
        if( p->ptr.p_int[i]!=i )
        {
            v = xb->ptr.p_complex[i];
            xb->ptr.p_complex[i] = xb->ptr.p_complex[p->ptr.p_int[i]];
            xb->ptr.p_complex[p->ptr.p_int[i]] = v;
        }
    }
    for(i=1; i<=n-1; i++)
    {
        v = ae_v_cdotproduct(&lua->ptr.pp_complex[i][0], 1, "N", &xb->ptr.p_complex[0], 1, "N", ae_v_len(0,i-1));
        xb->ptr.p_complex[i] = ae_c_sub(xb->ptr.p_complex[i],v);
    }
    xb->ptr.p_complex[n-1] = ae_c_div(xb->ptr.p_complex[n-1],ae_c_mul_d(lua->ptr.pp_complex[n-1][n-1],scalea));
    for(i=n-2; i>=0; i--)
    {
        ae_v_cmoved(&tmp->ptr.p_complex[i+1], 1, &lua->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), scalea);
        v = ae_v_cdotproduct(&tmp->ptr.p_complex[i+1], 1, "N", &xb->ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1));
        xb->ptr.p_complex[i] = ae_c_div(ae_c_sub(xb->ptr.p_complex[i],v),ae_c_mul_d(lua->ptr.pp_complex[i][i],scalea));
    }
}


/*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void densesolver_hpdbasiccholeskysolve(/* Complex */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* xb,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_complex v;


    
    /*
     * A = L*L' or A=U'*U
     */
    if( isupper )
    {
        
        /*
         * Solve U'*y=b first.
         */
        for(i=0; i<=n-1; i++)
        {
            xb->ptr.p_complex[i] = ae_c_div(xb->ptr.p_complex[i],ae_c_mul_d(ae_c_conj(cha->ptr.pp_complex[i][i], _state),sqrtscalea));
            if( i<n-1 )
            {
                v = xb->ptr.p_complex[i];
                ae_v_cmoved(&tmp->ptr.p_complex[i+1], 1, &cha->ptr.pp_complex[i][i+1], 1, "Conj", ae_v_len(i+1,n-1), sqrtscalea);
                ae_v_csubc(&xb->ptr.p_complex[i+1], 1, &tmp->ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1), v);
            }
        }
        
        /*
         * Solve U*x=y then.
         */
        for(i=n-1; i>=0; i--)
        {
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp->ptr.p_complex[i+1], 1, &cha->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), sqrtscalea);
                v = ae_v_cdotproduct(&tmp->ptr.p_complex[i+1], 1, "N", &xb->ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1));
                xb->ptr.p_complex[i] = ae_c_sub(xb->ptr.p_complex[i],v);
            }
            xb->ptr.p_complex[i] = ae_c_div(xb->ptr.p_complex[i],ae_c_mul_d(cha->ptr.pp_complex[i][i],sqrtscalea));
        }
    }
    else
    {
        
        /*
         * Solve L*y=b first
         */
        for(i=0; i<=n-1; i++)
        {
            if( i>0 )
            {
                ae_v_cmoved(&tmp->ptr.p_complex[0], 1, &cha->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i-1), sqrtscalea);
                v = ae_v_cdotproduct(&tmp->ptr.p_complex[0], 1, "N", &xb->ptr.p_complex[0], 1, "N", ae_v_len(0,i-1));
                xb->ptr.p_complex[i] = ae_c_sub(xb->ptr.p_complex[i],v);
            }
            xb->ptr.p_complex[i] = ae_c_div(xb->ptr.p_complex[i],ae_c_mul_d(cha->ptr.pp_complex[i][i],sqrtscalea));
        }
        
        /*
         * Solve L'*x=y then.
         */
        for(i=n-1; i>=0; i--)
        {
            xb->ptr.p_complex[i] = ae_c_div(xb->ptr.p_complex[i],ae_c_mul_d(ae_c_conj(cha->ptr.pp_complex[i][i], _state),sqrtscalea));
            if( i>0 )
            {
                v = xb->ptr.p_complex[i];
                ae_v_cmoved(&tmp->ptr.p_complex[0], 1, &cha->ptr.pp_complex[i][0], 1, "Conj", ae_v_len(0,i-1), sqrtscalea);
                ae_v_csubc(&xb->ptr.p_complex[0], 1, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(0,i-1), v);
            }
        }
    }
}


ae_bool _densesolverreport_init(densesolverreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _densesolverreport_init_copy(densesolverreport* dst, densesolverreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->r1 = src->r1;
    dst->rinf = src->rinf;
    return ae_true;
}


void _densesolverreport_clear(densesolverreport* p)
{
}


ae_bool _densesolverlsreport_init(densesolverlsreport* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_matrix_init(&p->cx, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _densesolverlsreport_init_copy(densesolverlsreport* dst, densesolverlsreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->r2 = src->r2;
    if( !ae_matrix_init_copy(&dst->cx, &src->cx, _state, make_automatic) )
        return ae_false;
    dst->n = src->n;
    dst->k = src->k;
    return ae_true;
}


void _densesolverlsreport_clear(densesolverlsreport* p)
{
    ae_matrix_clear(&p->cx);
}


/*$ End $*/
