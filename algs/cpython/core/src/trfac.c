/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
#include "trfac.h"


/*$ Declarations $*/
static void trfac_cmatrixluprec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static void trfac_rmatrixluprec(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static void trfac_cmatrixplurec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static void trfac_rmatrixplurec(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static void trfac_cmatrixlup2(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static void trfac_rmatrixlup2(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static void trfac_cmatrixplu2(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static void trfac_rmatrixplu2(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static ae_bool trfac_hpdmatrixcholeskyrec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static ae_bool trfac_hpdmatrixcholesky2(/* Complex */ ae_matrix* aaa,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tmp,
     ae_state *_state);
static ae_bool trfac_spdmatrixcholesky2(/* Real    */ ae_matrix* aaa,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
LU decomposition of a general real matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition.
It is optimized for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixlu(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{

    ae_vector_clear(pivots);

    ae_assert(m>0, "RMatrixLU: incorrect M!", _state);
    ae_assert(n>0, "RMatrixLU: incorrect N!", _state);
    rmatrixplu(a, m, n, pivots, _state);
}


/*************************************************************************
LU decomposition of a general complex matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition. It is optimized
for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixlu(/* Complex */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{

    ae_vector_clear(pivots);

    ae_assert(m>0, "CMatrixLU: incorrect M!", _state);
    ae_assert(n>0, "CMatrixLU: incorrect N!", _state);
    cmatrixplu(a, m, n, pivots, _state);
}


/*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U'*U  or A=L*L' (here X' detones conj(X^T)).

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U'*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
ae_bool hpdmatrixcholesky(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_COMPLEX, _state, ae_true);

    if( n<1 )
    {
        result = ae_false;
        ae_frame_leave(_state);
        return result;
    }
    result = trfac_hpdmatrixcholeskyrec(a, 0, n, isupper, &tmp, _state);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U^T*U  or A=L*L^T

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U^T*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
ae_bool spdmatrixcholesky(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    if( n<1 )
    {
        result = ae_false;
        ae_frame_leave(_state);
        return result;
    }
    result = spdmatrixcholeskyrec(a, 0, n, isupper, &tmp, _state);
    ae_frame_leave(_state);
    return result;
}


void rmatrixlup(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    double mx;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(pivots);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    
    /*
     * Internal LU decomposition subroutine.
     * Never call it directly.
     */
    ae_assert(m>0, "RMatrixLUP: incorrect M!", _state);
    ae_assert(n>0, "RMatrixLUP: incorrect N!", _state);
    
    /*
     * Scale matrix to avoid overflows,
     * decompose it, then scale back.
     */
    mx = 0;
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            mx = ae_maxreal(mx, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i=0; i<=m-1; i++)
        {
            ae_v_muld(&a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    ae_vector_set_length(pivots, ae_minint(m, n, _state), _state);
    ae_vector_set_length(&tmp, 2*ae_maxint(m, n, _state), _state);
    trfac_rmatrixluprec(a, 0, m, n, pivots, &tmp, _state);
    if( ae_fp_neq(mx,0) )
    {
        v = mx;
        for(i=0; i<=m-1; i++)
        {
            ae_v_muld(&a->ptr.pp_double[i][0], 1, ae_v_len(0,ae_minint(i, n-1, _state)), v);
        }
    }
    ae_frame_leave(_state);
}


void cmatrixlup(/* Complex */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    double mx;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(pivots);
    ae_vector_init(&tmp, 0, DT_COMPLEX, _state, ae_true);

    
    /*
     * Internal LU decomposition subroutine.
     * Never call it directly.
     */
    ae_assert(m>0, "CMatrixLUP: incorrect M!", _state);
    ae_assert(n>0, "CMatrixLUP: incorrect N!", _state);
    
    /*
     * Scale matrix to avoid overflows,
     * decompose it, then scale back.
     */
    mx = 0;
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            mx = ae_maxreal(mx, ae_c_abs(a->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i=0; i<=m-1; i++)
        {
            ae_v_cmuld(&a->ptr.pp_complex[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    ae_vector_set_length(pivots, ae_minint(m, n, _state), _state);
    ae_vector_set_length(&tmp, 2*ae_maxint(m, n, _state), _state);
    trfac_cmatrixluprec(a, 0, m, n, pivots, &tmp, _state);
    if( ae_fp_neq(mx,0) )
    {
        v = mx;
        for(i=0; i<=m-1; i++)
        {
            ae_v_cmuld(&a->ptr.pp_complex[i][0], 1, ae_v_len(0,ae_minint(i, n-1, _state)), v);
        }
    }
    ae_frame_leave(_state);
}


void rmatrixplu(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    double mx;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(pivots);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    
    /*
     * Internal LU decomposition subroutine.
     * Never call it directly.
     */
    ae_assert(m>0, "RMatrixPLU: incorrect M!", _state);
    ae_assert(n>0, "RMatrixPLU: incorrect N!", _state);
    ae_vector_set_length(&tmp, 2*ae_maxint(m, n, _state), _state);
    ae_vector_set_length(pivots, ae_minint(m, n, _state), _state);
    
    /*
     * Scale matrix to avoid overflows,
     * decompose it, then scale back.
     */
    mx = 0;
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            mx = ae_maxreal(mx, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
        }
    }
    if( ae_fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i=0; i<=m-1; i++)
        {
            ae_v_muld(&a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    trfac_rmatrixplurec(a, 0, m, n, pivots, &tmp, _state);
    if( ae_fp_neq(mx,0) )
    {
        v = mx;
        for(i=0; i<=ae_minint(m, n, _state)-1; i++)
        {
            ae_v_muld(&a->ptr.pp_double[i][i], 1, ae_v_len(i,n-1), v);
        }
    }
    ae_frame_leave(_state);
}


void cmatrixplu(/* Complex */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    double mx;
    ae_complex v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(pivots);
    ae_vector_init(&tmp, 0, DT_COMPLEX, _state, ae_true);

    
    /*
     * Internal LU decomposition subroutine.
     * Never call it directly.
     */
    ae_assert(m>0, "CMatrixPLU: incorrect M!", _state);
    ae_assert(n>0, "CMatrixPLU: incorrect N!", _state);
    ae_vector_set_length(&tmp, 2*ae_maxint(m, n, _state), _state);
    ae_vector_set_length(pivots, ae_minint(m, n, _state), _state);
    
    /*
     * Scale matrix to avoid overflows,
     * decompose it, then scale back.
     */
    mx = 0;
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            mx = ae_maxreal(mx, ae_c_abs(a->ptr.pp_complex[i][j], _state), _state);
        }
    }
    if( ae_fp_neq(mx,0) )
    {
        v = ae_complex_from_d(1/mx);
        for(i=0; i<=m-1; i++)
        {
            ae_v_cmulc(&a->ptr.pp_complex[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    trfac_cmatrixplurec(a, 0, m, n, pivots, &tmp, _state);
    if( ae_fp_neq(mx,0) )
    {
        v = ae_complex_from_d(mx);
        for(i=0; i<=ae_minint(m, n, _state)-1; i++)
        {
            ae_v_cmulc(&a->ptr.pp_complex[i][i], 1, ae_v_len(i,n-1), v);
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Recursive computational subroutine for SPDMatrixCholesky.

INPUT PARAMETERS:
    A       -   matrix given by upper or lower triangle
    Offs    -   offset of diagonal block to decompose
    N       -   diagonal block size
    IsUpper -   what half is given
    Tmp     -   temporary array; allocated by function, if its size is too
                small; can be reused on subsequent calls.
                
OUTPUT PARAMETERS:
    A       -   upper (or lower) triangle contains Cholesky decomposition

RESULT:
    True, on success
    False, on failure

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
ae_bool spdmatrixcholeskyrec(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t n1;
    ae_int_t n2;
    ae_bool result;


    
    /*
     * check N
     */
    if( n<1 )
    {
        result = ae_false;
        return result;
    }
    
    /*
     * Prepare buffer
     */
    if( tmp->cnt<2*n )
    {
        ae_vector_set_length(tmp, 2*n, _state);
    }
    
    /*
     * special cases
     */
    if( n==1 )
    {
        if( ae_fp_greater(a->ptr.pp_double[offs][offs],0) )
        {
            a->ptr.pp_double[offs][offs] = ae_sqrt(a->ptr.pp_double[offs][offs], _state);
            result = ae_true;
        }
        else
        {
            result = ae_false;
        }
        return result;
    }
    if( n<=ablasblocksize(a, _state) )
    {
        result = trfac_spdmatrixcholesky2(a, offs, n, isupper, tmp, _state);
        return result;
    }
    
    /*
     * general case: split task in cache-oblivious manner
     */
    result = ae_true;
    ablassplitlength(a, n, &n1, &n2, _state);
    result = spdmatrixcholeskyrec(a, offs, n1, isupper, tmp, _state);
    if( !result )
    {
        return result;
    }
    if( n2>0 )
    {
        if( isupper )
        {
            rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, ae_false, 1, a, offs, offs+n1, _state);
            rmatrixsyrk(n2, n1, -1.0, a, offs, offs+n1, 1, 1.0, a, offs+n1, offs+n1, isupper, _state);
        }
        else
        {
            rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, ae_false, 1, a, offs+n1, offs, _state);
            rmatrixsyrk(n2, n1, -1.0, a, offs+n1, offs, 0, 1.0, a, offs+n1, offs+n1, isupper, _state);
        }
        result = spdmatrixcholeskyrec(a, offs+n1, n2, isupper, tmp, _state);
        if( !result )
        {
            return result;
        }
    }
    return result;
}


/*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_cmatrixluprec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t m1;
    ae_int_t m2;


    
    /*
     * Kernel case
     */
    if( ae_minint(m, n, _state)<=ablascomplexblocksize(a, _state) )
    {
        trfac_cmatrixlup2(a, offs, m, n, pivots, tmp, _state);
        return;
    }
    
    /*
     * Preliminary step, make N>=M
     *
     *     ( A1 )
     * A = (    ), where A1 is square
     *     ( A2 )
     *
     * Factorize A1, update A2
     */
    if( m>n )
    {
        trfac_cmatrixluprec(a, offs, n, n, pivots, tmp, _state);
        for(i=0; i<=n-1; i++)
        {
            ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+n][offs+i], a->stride, "N", ae_v_len(0,m-n-1));
            ae_v_cmove(&a->ptr.pp_complex[offs+n][offs+i], a->stride, &a->ptr.pp_complex[offs+n][pivots->ptr.p_int[offs+i]], a->stride, "N", ae_v_len(offs+n,offs+m-1));
            ae_v_cmove(&a->ptr.pp_complex[offs+n][pivots->ptr.p_int[offs+i]], a->stride, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs+n,offs+m-1));
        }
        cmatrixrighttrsm(m-n, n, a, offs, offs, ae_true, ae_true, 0, a, offs+n, offs, _state);
        return;
    }
    
    /*
     * Non-kernel case
     */
    ablascomplexsplitlength(a, m, &m1, &m2, _state);
    trfac_cmatrixluprec(a, offs, m1, n, pivots, tmp, _state);
    if( m2>0 )
    {
        for(i=0; i<=m1-1; i++)
        {
            if( offs+i!=pivots->ptr.p_int[offs+i] )
            {
                ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+m1][offs+i], a->stride, "N", ae_v_len(0,m2-1));
                ae_v_cmove(&a->ptr.pp_complex[offs+m1][offs+i], a->stride, &a->ptr.pp_complex[offs+m1][pivots->ptr.p_int[offs+i]], a->stride, "N", ae_v_len(offs+m1,offs+m-1));
                ae_v_cmove(&a->ptr.pp_complex[offs+m1][pivots->ptr.p_int[offs+i]], a->stride, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs+m1,offs+m-1));
            }
        }
        cmatrixrighttrsm(m2, m1, a, offs, offs, ae_true, ae_true, 0, a, offs+m1, offs, _state);
        cmatrixgemm(m-m1, n-m1, m1, ae_complex_from_d(-1.0), a, offs+m1, offs, 0, a, offs, offs+m1, 0, ae_complex_from_d(1.0), a, offs+m1, offs+m1, _state);
        trfac_cmatrixluprec(a, offs+m1, m-m1, n-m1, pivots, tmp, _state);
        for(i=0; i<=m2-1; i++)
        {
            if( offs+m1+i!=pivots->ptr.p_int[offs+m1+i] )
            {
                ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs][offs+m1+i], a->stride, "N", ae_v_len(0,m1-1));
                ae_v_cmove(&a->ptr.pp_complex[offs][offs+m1+i], a->stride, &a->ptr.pp_complex[offs][pivots->ptr.p_int[offs+m1+i]], a->stride, "N", ae_v_len(offs,offs+m1-1));
                ae_v_cmove(&a->ptr.pp_complex[offs][pivots->ptr.p_int[offs+m1+i]], a->stride, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs,offs+m1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_rmatrixluprec(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t m1;
    ae_int_t m2;


    
    /*
     * Kernel case
     */
    if( ae_minint(m, n, _state)<=ablasblocksize(a, _state) )
    {
        trfac_rmatrixlup2(a, offs, m, n, pivots, tmp, _state);
        return;
    }
    
    /*
     * Preliminary step, make N>=M
     *
     *     ( A1 )
     * A = (    ), where A1 is square
     *     ( A2 )
     *
     * Factorize A1, update A2
     */
    if( m>n )
    {
        trfac_rmatrixluprec(a, offs, n, n, pivots, tmp, _state);
        for(i=0; i<=n-1; i++)
        {
            if( offs+i!=pivots->ptr.p_int[offs+i] )
            {
                ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+n][offs+i], a->stride, ae_v_len(0,m-n-1));
                ae_v_move(&a->ptr.pp_double[offs+n][offs+i], a->stride, &a->ptr.pp_double[offs+n][pivots->ptr.p_int[offs+i]], a->stride, ae_v_len(offs+n,offs+m-1));
                ae_v_move(&a->ptr.pp_double[offs+n][pivots->ptr.p_int[offs+i]], a->stride, &tmp->ptr.p_double[0], 1, ae_v_len(offs+n,offs+m-1));
            }
        }
        rmatrixrighttrsm(m-n, n, a, offs, offs, ae_true, ae_true, 0, a, offs+n, offs, _state);
        return;
    }
    
    /*
     * Non-kernel case
     */
    ablassplitlength(a, m, &m1, &m2, _state);
    trfac_rmatrixluprec(a, offs, m1, n, pivots, tmp, _state);
    if( m2>0 )
    {
        for(i=0; i<=m1-1; i++)
        {
            if( offs+i!=pivots->ptr.p_int[offs+i] )
            {
                ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+m1][offs+i], a->stride, ae_v_len(0,m2-1));
                ae_v_move(&a->ptr.pp_double[offs+m1][offs+i], a->stride, &a->ptr.pp_double[offs+m1][pivots->ptr.p_int[offs+i]], a->stride, ae_v_len(offs+m1,offs+m-1));
                ae_v_move(&a->ptr.pp_double[offs+m1][pivots->ptr.p_int[offs+i]], a->stride, &tmp->ptr.p_double[0], 1, ae_v_len(offs+m1,offs+m-1));
            }
        }
        rmatrixrighttrsm(m2, m1, a, offs, offs, ae_true, ae_true, 0, a, offs+m1, offs, _state);
        rmatrixgemm(m-m1, n-m1, m1, -1.0, a, offs+m1, offs, 0, a, offs, offs+m1, 0, 1.0, a, offs+m1, offs+m1, _state);
        trfac_rmatrixluprec(a, offs+m1, m-m1, n-m1, pivots, tmp, _state);
        for(i=0; i<=m2-1; i++)
        {
            if( offs+m1+i!=pivots->ptr.p_int[offs+m1+i] )
            {
                ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs][offs+m1+i], a->stride, ae_v_len(0,m1-1));
                ae_v_move(&a->ptr.pp_double[offs][offs+m1+i], a->stride, &a->ptr.pp_double[offs][pivots->ptr.p_int[offs+m1+i]], a->stride, ae_v_len(offs,offs+m1-1));
                ae_v_move(&a->ptr.pp_double[offs][pivots->ptr.p_int[offs+m1+i]], a->stride, &tmp->ptr.p_double[0], 1, ae_v_len(offs,offs+m1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_cmatrixplurec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;


    
    /*
     * Kernel case
     */
    if( ae_minint(m, n, _state)<=ablascomplexblocksize(a, _state) )
    {
        trfac_cmatrixplu2(a, offs, m, n, pivots, tmp, _state);
        return;
    }
    
    /*
     * Preliminary step, make M>=N.
     *
     * A = (A1 A2), where A1 is square
     * Factorize A1, update A2
     */
    if( n>m )
    {
        trfac_cmatrixplurec(a, offs, m, m, pivots, tmp, _state);
        for(i=0; i<=m-1; i++)
        {
            ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+i][offs+m], 1, "N", ae_v_len(0,n-m-1));
            ae_v_cmove(&a->ptr.pp_complex[offs+i][offs+m], 1, &a->ptr.pp_complex[pivots->ptr.p_int[offs+i]][offs+m], 1, "N", ae_v_len(offs+m,offs+n-1));
            ae_v_cmove(&a->ptr.pp_complex[pivots->ptr.p_int[offs+i]][offs+m], 1, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs+m,offs+n-1));
        }
        cmatrixlefttrsm(m, n-m, a, offs, offs, ae_false, ae_true, 0, a, offs, offs+m, _state);
        return;
    }
    
    /*
     * Non-kernel case
     */
    ablascomplexsplitlength(a, n, &n1, &n2, _state);
    trfac_cmatrixplurec(a, offs, m, n1, pivots, tmp, _state);
    if( n2>0 )
    {
        for(i=0; i<=n1-1; i++)
        {
            if( offs+i!=pivots->ptr.p_int[offs+i] )
            {
                ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+i][offs+n1], 1, "N", ae_v_len(0,n2-1));
                ae_v_cmove(&a->ptr.pp_complex[offs+i][offs+n1], 1, &a->ptr.pp_complex[pivots->ptr.p_int[offs+i]][offs+n1], 1, "N", ae_v_len(offs+n1,offs+n-1));
                ae_v_cmove(&a->ptr.pp_complex[pivots->ptr.p_int[offs+i]][offs+n1], 1, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs+n1,offs+n-1));
            }
        }
        cmatrixlefttrsm(n1, n2, a, offs, offs, ae_false, ae_true, 0, a, offs, offs+n1, _state);
        cmatrixgemm(m-n1, n-n1, n1, ae_complex_from_d(-1.0), a, offs+n1, offs, 0, a, offs, offs+n1, 0, ae_complex_from_d(1.0), a, offs+n1, offs+n1, _state);
        trfac_cmatrixplurec(a, offs+n1, m-n1, n-n1, pivots, tmp, _state);
        for(i=0; i<=n2-1; i++)
        {
            if( offs+n1+i!=pivots->ptr.p_int[offs+n1+i] )
            {
                ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+n1+i][offs], 1, "N", ae_v_len(0,n1-1));
                ae_v_cmove(&a->ptr.pp_complex[offs+n1+i][offs], 1, &a->ptr.pp_complex[pivots->ptr.p_int[offs+n1+i]][offs], 1, "N", ae_v_len(offs,offs+n1-1));
                ae_v_cmove(&a->ptr.pp_complex[pivots->ptr.p_int[offs+n1+i]][offs], 1, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs,offs+n1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_rmatrixplurec(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;


    
    /*
     * Kernel case
     */
    if( ae_minint(m, n, _state)<=ablasblocksize(a, _state) )
    {
        trfac_rmatrixplu2(a, offs, m, n, pivots, tmp, _state);
        return;
    }
    
    /*
     * Preliminary step, make M>=N.
     *
     * A = (A1 A2), where A1 is square
     * Factorize A1, update A2
     */
    if( n>m )
    {
        trfac_rmatrixplurec(a, offs, m, m, pivots, tmp, _state);
        for(i=0; i<=m-1; i++)
        {
            ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+i][offs+m], 1, ae_v_len(0,n-m-1));
            ae_v_move(&a->ptr.pp_double[offs+i][offs+m], 1, &a->ptr.pp_double[pivots->ptr.p_int[offs+i]][offs+m], 1, ae_v_len(offs+m,offs+n-1));
            ae_v_move(&a->ptr.pp_double[pivots->ptr.p_int[offs+i]][offs+m], 1, &tmp->ptr.p_double[0], 1, ae_v_len(offs+m,offs+n-1));
        }
        rmatrixlefttrsm(m, n-m, a, offs, offs, ae_false, ae_true, 0, a, offs, offs+m, _state);
        return;
    }
    
    /*
     * Non-kernel case
     */
    ablassplitlength(a, n, &n1, &n2, _state);
    trfac_rmatrixplurec(a, offs, m, n1, pivots, tmp, _state);
    if( n2>0 )
    {
        for(i=0; i<=n1-1; i++)
        {
            if( offs+i!=pivots->ptr.p_int[offs+i] )
            {
                ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+i][offs+n1], 1, ae_v_len(0,n2-1));
                ae_v_move(&a->ptr.pp_double[offs+i][offs+n1], 1, &a->ptr.pp_double[pivots->ptr.p_int[offs+i]][offs+n1], 1, ae_v_len(offs+n1,offs+n-1));
                ae_v_move(&a->ptr.pp_double[pivots->ptr.p_int[offs+i]][offs+n1], 1, &tmp->ptr.p_double[0], 1, ae_v_len(offs+n1,offs+n-1));
            }
        }
        rmatrixlefttrsm(n1, n2, a, offs, offs, ae_false, ae_true, 0, a, offs, offs+n1, _state);
        rmatrixgemm(m-n1, n-n1, n1, -1.0, a, offs+n1, offs, 0, a, offs, offs+n1, 0, 1.0, a, offs+n1, offs+n1, _state);
        trfac_rmatrixplurec(a, offs+n1, m-n1, n-n1, pivots, tmp, _state);
        for(i=0; i<=n2-1; i++)
        {
            if( offs+n1+i!=pivots->ptr.p_int[offs+n1+i] )
            {
                ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+n1+i][offs], 1, ae_v_len(0,n1-1));
                ae_v_move(&a->ptr.pp_double[offs+n1+i][offs], 1, &a->ptr.pp_double[pivots->ptr.p_int[offs+n1+i]][offs], 1, ae_v_len(offs,offs+n1-1));
                ae_v_move(&a->ptr.pp_double[pivots->ptr.p_int[offs+n1+i]][offs], 1, &tmp->ptr.p_double[0], 1, ae_v_len(offs,offs+n1-1));
            }
        }
    }
}


/*************************************************************************
Complex LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_cmatrixlup2(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t jp;
    ae_complex s;


    
    /*
     * Quick return if possible
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * main cycle
     */
    for(j=0; j<=ae_minint(m-1, n-1, _state); j++)
    {
        
        /*
         * Find pivot, swap columns
         */
        jp = j;
        for(i=j+1; i<=n-1; i++)
        {
            if( ae_fp_greater(ae_c_abs(a->ptr.pp_complex[offs+j][offs+i], _state),ae_c_abs(a->ptr.pp_complex[offs+j][offs+jp], _state)) )
            {
                jp = i;
            }
        }
        pivots->ptr.p_int[offs+j] = offs+jp;
        if( jp!=j )
        {
            ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs][offs+j], a->stride, "N", ae_v_len(0,m-1));
            ae_v_cmove(&a->ptr.pp_complex[offs][offs+j], a->stride, &a->ptr.pp_complex[offs][offs+jp], a->stride, "N", ae_v_len(offs,offs+m-1));
            ae_v_cmove(&a->ptr.pp_complex[offs][offs+jp], a->stride, &tmp->ptr.p_complex[0], 1, "N", ae_v_len(offs,offs+m-1));
        }
        
        /*
         * LU decomposition of 1x(N-J) matrix
         */
        if( ae_c_neq_d(a->ptr.pp_complex[offs+j][offs+j],0)&&j+1<=n-1 )
        {
            s = ae_c_d_div(1,a->ptr.pp_complex[offs+j][offs+j]);
            ae_v_cmulc(&a->ptr.pp_complex[offs+j][offs+j+1], 1, ae_v_len(offs+j+1,offs+n-1), s);
        }
        
        /*
         * Update trailing (M-J-1)x(N-J-1) matrix
         */
        if( j<ae_minint(m-1, n-1, _state) )
        {
            ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+j+1][offs+j], a->stride, "N", ae_v_len(0,m-j-2));
            ae_v_cmoveneg(&tmp->ptr.p_complex[m], 1, &a->ptr.pp_complex[offs+j][offs+j+1], 1, "N", ae_v_len(m,m+n-j-2));
            cmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m, _state);
        }
    }
}


/*************************************************************************
Real LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
static void trfac_rmatrixlup2(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t jp;
    double s;


    
    /*
     * Quick return if possible
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * main cycle
     */
    for(j=0; j<=ae_minint(m-1, n-1, _state); j++)
    {
        
        /*
         * Find pivot, swap columns
         */
        jp = j;
        for(i=j+1; i<=n-1; i++)
        {
            if( ae_fp_greater(ae_fabs(a->ptr.pp_double[offs+j][offs+i], _state),ae_fabs(a->ptr.pp_double[offs+j][offs+jp], _state)) )
            {
                jp = i;
            }
        }
        pivots->ptr.p_int[offs+j] = offs+jp;
        if( jp!=j )
        {
            ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs][offs+j], a->stride, ae_v_len(0,m-1));
            ae_v_move(&a->ptr.pp_double[offs][offs+j], a->stride, &a->ptr.pp_double[offs][offs+jp], a->stride, ae_v_len(offs,offs+m-1));
            ae_v_move(&a->ptr.pp_double[offs][offs+jp], a->stride, &tmp->ptr.p_double[0], 1, ae_v_len(offs,offs+m-1));
        }
        
        /*
         * LU decomposition of 1x(N-J) matrix
         */
        if( ae_fp_neq(a->ptr.pp_double[offs+j][offs+j],0)&&j+1<=n-1 )
        {
            s = 1/a->ptr.pp_double[offs+j][offs+j];
            ae_v_muld(&a->ptr.pp_double[offs+j][offs+j+1], 1, ae_v_len(offs+j+1,offs+n-1), s);
        }
        
        /*
         * Update trailing (M-J-1)x(N-J-1) matrix
         */
        if( j<ae_minint(m-1, n-1, _state) )
        {
            ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+j+1][offs+j], a->stride, ae_v_len(0,m-j-2));
            ae_v_moveneg(&tmp->ptr.p_double[m], 1, &a->ptr.pp_double[offs+j][offs+j+1], 1, ae_v_len(m,m+n-j-2));
            rmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m, _state);
        }
    }
}


/*************************************************************************
Complex PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
static void trfac_cmatrixplu2(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t jp;
    ae_complex s;


    
    /*
     * Quick return if possible
     */
    if( m==0||n==0 )
    {
        return;
    }
    for(j=0; j<=ae_minint(m-1, n-1, _state); j++)
    {
        
        /*
         * Find pivot and test for singularity.
         */
        jp = j;
        for(i=j+1; i<=m-1; i++)
        {
            if( ae_fp_greater(ae_c_abs(a->ptr.pp_complex[offs+i][offs+j], _state),ae_c_abs(a->ptr.pp_complex[offs+jp][offs+j], _state)) )
            {
                jp = i;
            }
        }
        pivots->ptr.p_int[offs+j] = offs+jp;
        if( ae_c_neq_d(a->ptr.pp_complex[offs+jp][offs+j],0) )
        {
            
            /*
             *Apply the interchange to rows
             */
            if( jp!=j )
            {
                for(i=0; i<=n-1; i++)
                {
                    s = a->ptr.pp_complex[offs+j][offs+i];
                    a->ptr.pp_complex[offs+j][offs+i] = a->ptr.pp_complex[offs+jp][offs+i];
                    a->ptr.pp_complex[offs+jp][offs+i] = s;
                }
            }
            
            /*
             *Compute elements J+1:M of J-th column.
             */
            if( j+1<=m-1 )
            {
                s = ae_c_d_div(1,a->ptr.pp_complex[offs+j][offs+j]);
                ae_v_cmulc(&a->ptr.pp_complex[offs+j+1][offs+j], a->stride, ae_v_len(offs+j+1,offs+m-1), s);
            }
        }
        if( j<ae_minint(m, n, _state)-1 )
        {
            
            /*
             *Update trailing submatrix.
             */
            ae_v_cmove(&tmp->ptr.p_complex[0], 1, &a->ptr.pp_complex[offs+j+1][offs+j], a->stride, "N", ae_v_len(0,m-j-2));
            ae_v_cmoveneg(&tmp->ptr.p_complex[m], 1, &a->ptr.pp_complex[offs+j][offs+j+1], 1, "N", ae_v_len(m,m+n-j-2));
            cmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m, _state);
        }
    }
}


/*************************************************************************
Real PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
static void trfac_rmatrixplu2(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t jp;
    double s;


    
    /*
     * Quick return if possible
     */
    if( m==0||n==0 )
    {
        return;
    }
    for(j=0; j<=ae_minint(m-1, n-1, _state); j++)
    {
        
        /*
         * Find pivot and test for singularity.
         */
        jp = j;
        for(i=j+1; i<=m-1; i++)
        {
            if( ae_fp_greater(ae_fabs(a->ptr.pp_double[offs+i][offs+j], _state),ae_fabs(a->ptr.pp_double[offs+jp][offs+j], _state)) )
            {
                jp = i;
            }
        }
        pivots->ptr.p_int[offs+j] = offs+jp;
        if( ae_fp_neq(a->ptr.pp_double[offs+jp][offs+j],0) )
        {
            
            /*
             *Apply the interchange to rows
             */
            if( jp!=j )
            {
                for(i=0; i<=n-1; i++)
                {
                    s = a->ptr.pp_double[offs+j][offs+i];
                    a->ptr.pp_double[offs+j][offs+i] = a->ptr.pp_double[offs+jp][offs+i];
                    a->ptr.pp_double[offs+jp][offs+i] = s;
                }
            }
            
            /*
             *Compute elements J+1:M of J-th column.
             */
            if( j+1<=m-1 )
            {
                s = 1/a->ptr.pp_double[offs+j][offs+j];
                ae_v_muld(&a->ptr.pp_double[offs+j+1][offs+j], a->stride, ae_v_len(offs+j+1,offs+m-1), s);
            }
        }
        if( j<ae_minint(m, n, _state)-1 )
        {
            
            /*
             *Update trailing submatrix.
             */
            ae_v_move(&tmp->ptr.p_double[0], 1, &a->ptr.pp_double[offs+j+1][offs+j], a->stride, ae_v_len(0,m-j-2));
            ae_v_moveneg(&tmp->ptr.p_double[m], 1, &a->ptr.pp_double[offs+j][offs+j+1], 1, ae_v_len(m,m+n-j-2));
            rmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m, _state);
        }
    }
}


/*************************************************************************
Recursive computational subroutine for HPDMatrixCholesky

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
static ae_bool trfac_hpdmatrixcholeskyrec(/* Complex */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t n1;
    ae_int_t n2;
    ae_bool result;


    
    /*
     * check N
     */
    if( n<1 )
    {
        result = ae_false;
        return result;
    }
    
    /*
     * Prepare buffer
     */
    if( tmp->cnt<2*n )
    {
        ae_vector_set_length(tmp, 2*n, _state);
    }
    
    /*
     * special cases
     */
    if( n==1 )
    {
        if( ae_fp_greater(a->ptr.pp_complex[offs][offs].x,0) )
        {
            a->ptr.pp_complex[offs][offs] = ae_complex_from_d(ae_sqrt(a->ptr.pp_complex[offs][offs].x, _state));
            result = ae_true;
        }
        else
        {
            result = ae_false;
        }
        return result;
    }
    if( n<=ablascomplexblocksize(a, _state) )
    {
        result = trfac_hpdmatrixcholesky2(a, offs, n, isupper, tmp, _state);
        return result;
    }
    
    /*
     * general case: split task in cache-oblivious manner
     */
    result = ae_true;
    ablascomplexsplitlength(a, n, &n1, &n2, _state);
    result = trfac_hpdmatrixcholeskyrec(a, offs, n1, isupper, tmp, _state);
    if( !result )
    {
        return result;
    }
    if( n2>0 )
    {
        if( isupper )
        {
            cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, ae_false, 2, a, offs, offs+n1, _state);
            cmatrixsyrk(n2, n1, -1.0, a, offs, offs+n1, 2, 1.0, a, offs+n1, offs+n1, isupper, _state);
        }
        else
        {
            cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, ae_false, 2, a, offs+n1, offs, _state);
            cmatrixsyrk(n2, n1, -1.0, a, offs+n1, offs, 0, 1.0, a, offs+n1, offs+n1, isupper, _state);
        }
        result = trfac_hpdmatrixcholeskyrec(a, offs+n1, n2, isupper, tmp, _state);
        if( !result )
        {
            return result;
        }
    }
    return result;
}


/*************************************************************************
Level-2 Hermitian Cholesky subroutine.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static ae_bool trfac_hpdmatrixcholesky2(/* Complex */ ae_matrix* aaa,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double ajj;
    ae_complex v;
    double r;
    ae_bool result;


    result = ae_true;
    if( n<0 )
    {
        result = ae_false;
        return result;
    }
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return result;
    }
    if( isupper )
    {
        
        /*
         * Compute the Cholesky factorization A = U'*U.
         */
        for(j=0; j<=n-1; j++)
        {
            
            /*
             * Compute U(J,J) and test for non-positive-definiteness.
             */
            v = ae_v_cdotproduct(&aaa->ptr.pp_complex[offs][offs+j], aaa->stride, "Conj", &aaa->ptr.pp_complex[offs][offs+j], aaa->stride, "N", ae_v_len(offs,offs+j-1));
            ajj = ae_c_sub(aaa->ptr.pp_complex[offs+j][offs+j],v).x;
            if( ae_fp_less_eq(ajj,0) )
            {
                aaa->ptr.pp_complex[offs+j][offs+j] = ae_complex_from_d(ajj);
                result = ae_false;
                return result;
            }
            ajj = ae_sqrt(ajj, _state);
            aaa->ptr.pp_complex[offs+j][offs+j] = ae_complex_from_d(ajj);
            
            /*
             * Compute elements J+1:N-1 of row J.
             */
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ae_v_cmoveneg(&tmp->ptr.p_complex[0], 1, &aaa->ptr.pp_complex[offs][offs+j], aaa->stride, "Conj", ae_v_len(0,j-1));
                    cmatrixmv(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n, _state);
                    ae_v_cadd(&aaa->ptr.pp_complex[offs+j][offs+j+1], 1, &tmp->ptr.p_complex[n], 1, "N", ae_v_len(offs+j+1,offs+n-1));
                }
                r = 1/ajj;
                ae_v_cmuld(&aaa->ptr.pp_complex[offs+j][offs+j+1], 1, ae_v_len(offs+j+1,offs+n-1), r);
            }
        }
    }
    else
    {
        
        /*
         * Compute the Cholesky factorization A = L*L'.
         */
        for(j=0; j<=n-1; j++)
        {
            
            /*
             * Compute L(J+1,J+1) and test for non-positive-definiteness.
             */
            v = ae_v_cdotproduct(&aaa->ptr.pp_complex[offs+j][offs], 1, "Conj", &aaa->ptr.pp_complex[offs+j][offs], 1, "N", ae_v_len(offs,offs+j-1));
            ajj = ae_c_sub(aaa->ptr.pp_complex[offs+j][offs+j],v).x;
            if( ae_fp_less_eq(ajj,0) )
            {
                aaa->ptr.pp_complex[offs+j][offs+j] = ae_complex_from_d(ajj);
                result = ae_false;
                return result;
            }
            ajj = ae_sqrt(ajj, _state);
            aaa->ptr.pp_complex[offs+j][offs+j] = ae_complex_from_d(ajj);
            
            /*
             * Compute elements J+1:N of column J.
             */
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ae_v_cmove(&tmp->ptr.p_complex[0], 1, &aaa->ptr.pp_complex[offs+j][offs], 1, "Conj", ae_v_len(0,j-1));
                    cmatrixmv(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n, _state);
                    for(i=0; i<=n-j-2; i++)
                    {
                        aaa->ptr.pp_complex[offs+j+1+i][offs+j] = ae_c_div_d(ae_c_sub(aaa->ptr.pp_complex[offs+j+1+i][offs+j],tmp->ptr.p_complex[n+i]),ajj);
                    }
                }
                else
                {
                    for(i=0; i<=n-j-2; i++)
                    {
                        aaa->ptr.pp_complex[offs+j+1+i][offs+j] = ae_c_div_d(aaa->ptr.pp_complex[offs+j+1+i][offs+j],ajj);
                    }
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Level-2 Cholesky subroutine

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static ae_bool trfac_spdmatrixcholesky2(/* Real    */ ae_matrix* aaa,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double ajj;
    double v;
    double r;
    ae_bool result;


    result = ae_true;
    if( n<0 )
    {
        result = ae_false;
        return result;
    }
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return result;
    }
    if( isupper )
    {
        
        /*
         * Compute the Cholesky factorization A = U'*U.
         */
        for(j=0; j<=n-1; j++)
        {
            
            /*
             * Compute U(J,J) and test for non-positive-definiteness.
             */
            v = ae_v_dotproduct(&aaa->ptr.pp_double[offs][offs+j], aaa->stride, &aaa->ptr.pp_double[offs][offs+j], aaa->stride, ae_v_len(offs,offs+j-1));
            ajj = aaa->ptr.pp_double[offs+j][offs+j]-v;
            if( ae_fp_less_eq(ajj,0) )
            {
                aaa->ptr.pp_double[offs+j][offs+j] = ajj;
                result = ae_false;
                return result;
            }
            ajj = ae_sqrt(ajj, _state);
            aaa->ptr.pp_double[offs+j][offs+j] = ajj;
            
            /*
             * Compute elements J+1:N-1 of row J.
             */
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ae_v_moveneg(&tmp->ptr.p_double[0], 1, &aaa->ptr.pp_double[offs][offs+j], aaa->stride, ae_v_len(0,j-1));
                    rmatrixmv(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n, _state);
                    ae_v_add(&aaa->ptr.pp_double[offs+j][offs+j+1], 1, &tmp->ptr.p_double[n], 1, ae_v_len(offs+j+1,offs+n-1));
                }
                r = 1/ajj;
                ae_v_muld(&aaa->ptr.pp_double[offs+j][offs+j+1], 1, ae_v_len(offs+j+1,offs+n-1), r);
            }
        }
    }
    else
    {
        
        /*
         * Compute the Cholesky factorization A = L*L'.
         */
        for(j=0; j<=n-1; j++)
        {
            
            /*
             * Compute L(J+1,J+1) and test for non-positive-definiteness.
             */
            v = ae_v_dotproduct(&aaa->ptr.pp_double[offs+j][offs], 1, &aaa->ptr.pp_double[offs+j][offs], 1, ae_v_len(offs,offs+j-1));
            ajj = aaa->ptr.pp_double[offs+j][offs+j]-v;
            if( ae_fp_less_eq(ajj,0) )
            {
                aaa->ptr.pp_double[offs+j][offs+j] = ajj;
                result = ae_false;
                return result;
            }
            ajj = ae_sqrt(ajj, _state);
            aaa->ptr.pp_double[offs+j][offs+j] = ajj;
            
            /*
             * Compute elements J+1:N of column J.
             */
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ae_v_move(&tmp->ptr.p_double[0], 1, &aaa->ptr.pp_double[offs+j][offs], 1, ae_v_len(0,j-1));
                    rmatrixmv(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n, _state);
                    for(i=0; i<=n-j-2; i++)
                    {
                        aaa->ptr.pp_double[offs+j+1+i][offs+j] = (aaa->ptr.pp_double[offs+j+1+i][offs+j]-tmp->ptr.p_double[n+i])/ajj;
                    }
                }
                else
                {
                    for(i=0; i<=n-j-2; i++)
                    {
                        aaa->ptr.pp_double[offs+j+1+i][offs+j] = aaa->ptr.pp_double[offs+j+1+i][offs+j]/ajj;
                    }
                }
            }
        }
    }
    return result;
}


/*$ End $*/
