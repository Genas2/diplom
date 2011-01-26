/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

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
#include "bdsvd.h"


/*$ Declarations $*/
static ae_bool bdsvd_bidiagonalsvddecompositioninternal(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     ae_bool isfractionalaccuracyrequired,
     /* Real    */ ae_matrix* u,
     ae_int_t ustart,
     ae_int_t nru,
     /* Real    */ ae_matrix* c,
     ae_int_t cstart,
     ae_int_t ncc,
     /* Real    */ ae_matrix* vt,
     ae_int_t vstart,
     ae_int_t ncvt,
     ae_state *_state);
static double bdsvd_extsignbdsqr(double a, double b, ae_state *_state);
static void bdsvd_svd2x2(double f,
     double g,
     double h,
     double* ssmin,
     double* ssmax,
     ae_state *_state);
static void bdsvd_svdv2x2(double f,
     double g,
     double h,
     double* ssmin,
     double* ssmax,
     double* snr,
     double* csr,
     double* snl,
     double* csl,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Singular value decomposition of a bidiagonal matrix (extended algorithm)

The algorithm performs the singular value decomposition  of  a  bidiagonal
matrix B (upper or lower) representing it as B = Q*S*P^T, where Q and  P -
orthogonal matrices, S - diagonal matrix with non-negative elements on the
main diagonal, in descending order.

The  algorithm  finds  singular  values.  In  addition,  the algorithm can
calculate  matrices  Q  and P (more precisely, not the matrices, but their
product  with  given  matrices U and VT - U*Q and (P^T)*VT)).  Of  course,
matrices U and VT can be of any type, including identity. Furthermore, the
algorithm can calculate Q'*C (this product is calculated more  effectively
than U*Q,  because  this calculation operates with rows instead  of matrix
columns).

The feature of the algorithm is its ability to find  all  singular  values
including those which are arbitrarily close to 0  with  relative  accuracy
close to  machine precision. If the parameter IsFractionalAccuracyRequired
is set to True, all singular values will have high relative accuracy close
to machine precision. If the parameter is set to False, only  the  biggest
singular value will have relative accuracy  close  to  machine  precision.
The absolute error of other singular values is equal to the absolute error
of the biggest singular value.

Input parameters:
    D       -   main diagonal of matrix B.
                Array whose index ranges within [0..N-1].
    E       -   superdiagonal (or subdiagonal) of matrix B.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix B.
    IsUpper -   True, if the matrix is upper bidiagonal.
    IsFractionalAccuracyRequired -
                accuracy to search singular values with.
    U       -   matrix to be multiplied by Q.
                Array whose indexes range within [0..NRU-1, 0..N-1].
                The matrix can be bigger, in that case only the  submatrix
                [0..NRU-1, 0..N-1] will be multiplied by Q.
    NRU     -   number of rows in matrix U.
    C       -   matrix to be multiplied by Q'.
                Array whose indexes range within [0..N-1, 0..NCC-1].
                The matrix can be bigger, in that case only the  submatrix
                [0..N-1, 0..NCC-1] will be multiplied by Q'.
    NCC     -   number of columns in matrix C.
    VT      -   matrix to be multiplied by P^T.
                Array whose indexes range within [0..N-1, 0..NCVT-1].
                The matrix can be bigger, in that case only the  submatrix
                [0..N-1, 0..NCVT-1] will be multiplied by P^T.
    NCVT    -   number of columns in matrix VT.

Output parameters:
    D       -   singular values of matrix B in descending order.
    U       -   if NRU>0, contains matrix U*Q.
    VT      -   if NCVT>0, contains matrix (P^T)*VT.
    C       -   if NCC>0, contains matrix Q'*C.

Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged (rare case).

Additional information:
    The type of convergence is controlled by the internal  parameter  TOL.
    If the parameter is greater than 0, the singular values will have
    relative accuracy TOL. If TOL<0, the singular values will have
    absolute accuracy ABS(TOL)*norm(B).
    By default, |TOL| falls within the range of 10*Epsilon and 100*Epsilon,
    where Epsilon is the machine precision. It is not  recommended  to  use
    TOL less than 10*Epsilon since this will  considerably  slow  down  the
    algorithm and may not lead to error decreasing.
History:
    * 31 March, 2007.
        changed MAXITR from 6 to 12.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1999.
*************************************************************************/
ae_bool rmatrixbdsvd(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     ae_bool isfractionalaccuracyrequired,
     /* Real    */ ae_matrix* u,
     ae_int_t nru,
     /* Real    */ ae_matrix* c,
     ae_int_t ncc,
     /* Real    */ ae_matrix* vt,
     ae_int_t ncvt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _e;
    ae_vector d1;
    ae_vector e1;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_e, e, _state, ae_true);
    e = &_e;
    ae_vector_init(&d1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&e1, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&d1, n+1, _state);
    ae_v_move(&d1.ptr.p_double[1], 1, &d->ptr.p_double[0], 1, ae_v_len(1,n));
    if( n>1 )
    {
        ae_vector_set_length(&e1, n-1+1, _state);
        ae_v_move(&e1.ptr.p_double[1], 1, &e->ptr.p_double[0], 1, ae_v_len(1,n-1));
    }
    result = bdsvd_bidiagonalsvddecompositioninternal(&d1, &e1, n, isupper, isfractionalaccuracyrequired, u, 0, nru, c, 0, ncc, vt, 0, ncvt, _state);
    ae_v_move(&d->ptr.p_double[0], 1, &d1.ptr.p_double[1], 1, ae_v_len(0,n-1));
    ae_frame_leave(_state);
    return result;
}


ae_bool bidiagonalsvddecomposition(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     ae_bool isfractionalaccuracyrequired,
     /* Real    */ ae_matrix* u,
     ae_int_t nru,
     /* Real    */ ae_matrix* c,
     ae_int_t ncc,
     /* Real    */ ae_matrix* vt,
     ae_int_t ncvt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _e;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_e, e, _state, ae_true);
    e = &_e;

    result = bdsvd_bidiagonalsvddecompositioninternal(d, e, n, isupper, isfractionalaccuracyrequired, u, 1, nru, c, 1, ncc, vt, 1, ncvt, _state);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Internal working subroutine for bidiagonal decomposition
*************************************************************************/
static ae_bool bdsvd_bidiagonalsvddecompositioninternal(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     ae_bool isfractionalaccuracyrequired,
     /* Real    */ ae_matrix* u,
     ae_int_t ustart,
     ae_int_t nru,
     /* Real    */ ae_matrix* c,
     ae_int_t cstart,
     ae_int_t ncc,
     /* Real    */ ae_matrix* vt,
     ae_int_t vstart,
     ae_int_t ncvt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _e;
    ae_int_t i;
    ae_int_t idir;
    ae_int_t isub;
    ae_int_t iter;
    ae_int_t j;
    ae_int_t ll;
    ae_int_t lll;
    ae_int_t m;
    ae_int_t maxit;
    ae_int_t oldll;
    ae_int_t oldm;
    double abse;
    double abss;
    double cosl;
    double cosr;
    double cs;
    double eps;
    double f;
    double g;
    double h;
    double mu;
    double oldcs;
    double oldsn;
    double r;
    double shift;
    double sigmn;
    double sigmx;
    double sinl;
    double sinr;
    double sll;
    double smax;
    double smin;
    double sminl;
    double sminlo;
    double sminoa;
    double sn;
    double thresh;
    double tol;
    double tolmul;
    double unfl;
    ae_vector work0;
    ae_vector work1;
    ae_vector work2;
    ae_vector work3;
    ae_int_t maxitr;
    ae_bool matrixsplitflag;
    ae_bool iterflag;
    ae_vector utemp;
    ae_vector vttemp;
    ae_vector ctemp;
    ae_vector etemp;
    ae_bool rightside;
    ae_bool fwddir;
    double tmp;
    ae_int_t mm1;
    ae_int_t mm0;
    ae_bool bchangedir;
    ae_int_t uend;
    ae_int_t cend;
    ae_int_t vend;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_e, e, _state, ae_true);
    e = &_e;
    ae_vector_init(&work0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work3, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&utemp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&vttemp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ctemp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&etemp, 0, DT_REAL, _state, ae_true);

    result = ae_true;
    if( n==0 )
    {
        ae_frame_leave(_state);
        return result;
    }
    if( n==1 )
    {
        if( ae_fp_less(d->ptr.p_double[1],0) )
        {
            d->ptr.p_double[1] = -d->ptr.p_double[1];
            if( ncvt>0 )
            {
                ae_v_muld(&vt->ptr.pp_double[vstart][vstart], 1, ae_v_len(vstart,vstart+ncvt-1), -1);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * these initializers are not really necessary,
     * but without them compiler complains about uninitialized locals
     */
    ll = 0;
    oldsn = 0;
    
    /*
     * init
     */
    ae_vector_set_length(&work0, n-1+1, _state);
    ae_vector_set_length(&work1, n-1+1, _state);
    ae_vector_set_length(&work2, n-1+1, _state);
    ae_vector_set_length(&work3, n-1+1, _state);
    uend = ustart+ae_maxint(nru-1, 0, _state);
    vend = vstart+ae_maxint(ncvt-1, 0, _state);
    cend = cstart+ae_maxint(ncc-1, 0, _state);
    ae_vector_set_length(&utemp, uend+1, _state);
    ae_vector_set_length(&vttemp, vend+1, _state);
    ae_vector_set_length(&ctemp, cend+1, _state);
    maxitr = 12;
    rightside = ae_true;
    fwddir = ae_true;
    
    /*
     * resize E from N-1 to N
     */
    ae_vector_set_length(&etemp, n+1, _state);
    for(i=1; i<=n-1; i++)
    {
        etemp.ptr.p_double[i] = e->ptr.p_double[i];
    }
    ae_vector_set_length(e, n+1, _state);
    for(i=1; i<=n-1; i++)
    {
        e->ptr.p_double[i] = etemp.ptr.p_double[i];
    }
    e->ptr.p_double[n] = 0;
    idir = 0;
    
    /*
     * Get machine constants
     */
    eps = ae_machineepsilon;
    unfl = ae_minrealnumber;
    
    /*
     * If matrix lower bidiagonal, rotate to be upper bidiagonal
     * by applying Givens rotations on the left
     */
    if( !isupper )
    {
        for(i=1; i<=n-1; i++)
        {
            generaterotation(d->ptr.p_double[i], e->ptr.p_double[i], &cs, &sn, &r, _state);
            d->ptr.p_double[i] = r;
            e->ptr.p_double[i] = sn*d->ptr.p_double[i+1];
            d->ptr.p_double[i+1] = cs*d->ptr.p_double[i+1];
            work0.ptr.p_double[i] = cs;
            work1.ptr.p_double[i] = sn;
        }
        
        /*
         * Update singular vectors if desired
         */
        if( nru>0 )
        {
            applyrotationsfromtheright(fwddir, ustart, uend, 1+ustart-1, n+ustart-1, &work0, &work1, u, &utemp, _state);
        }
        if( ncc>0 )
        {
            applyrotationsfromtheleft(fwddir, 1+cstart-1, n+cstart-1, cstart, cend, &work0, &work1, c, &ctemp, _state);
        }
    }
    
    /*
     * Compute singular values to relative accuracy TOL
     * (By setting TOL to be negative, algorithm will compute
     * singular values to absolute accuracy ABS(TOL)*norm(input matrix))
     */
    tolmul = ae_maxreal(10, ae_minreal(100, ae_pow(eps, -0.125, _state), _state), _state);
    tol = tolmul*eps;
    if( !isfractionalaccuracyrequired )
    {
        tol = -tol;
    }
    
    /*
     * Compute approximate maximum, minimum singular values
     */
    smax = 0;
    for(i=1; i<=n; i++)
    {
        smax = ae_maxreal(smax, ae_fabs(d->ptr.p_double[i], _state), _state);
    }
    for(i=1; i<=n-1; i++)
    {
        smax = ae_maxreal(smax, ae_fabs(e->ptr.p_double[i], _state), _state);
    }
    sminl = 0;
    if( ae_fp_greater_eq(tol,0) )
    {
        
        /*
         * Relative accuracy desired
         */
        sminoa = ae_fabs(d->ptr.p_double[1], _state);
        if( ae_fp_neq(sminoa,0) )
        {
            mu = sminoa;
            for(i=2; i<=n; i++)
            {
                mu = ae_fabs(d->ptr.p_double[i], _state)*(mu/(mu+ae_fabs(e->ptr.p_double[i-1], _state)));
                sminoa = ae_minreal(sminoa, mu, _state);
                if( ae_fp_eq(sminoa,0) )
                {
                    break;
                }
            }
        }
        sminoa = sminoa/ae_sqrt(n, _state);
        thresh = ae_maxreal(tol*sminoa, maxitr*n*n*unfl, _state);
    }
    else
    {
        
        /*
         * Absolute accuracy desired
         */
        thresh = ae_maxreal(ae_fabs(tol, _state)*smax, maxitr*n*n*unfl, _state);
    }
    
    /*
     * Prepare for main iteration loop for the singular values
     * (MAXIT is the maximum number of passes through the inner
     * loop permitted before nonconvergence signalled.)
     */
    maxit = maxitr*n*n;
    iter = 0;
    oldll = -1;
    oldm = -1;
    
    /*
     * M points to last element of unconverged part of matrix
     */
    m = n;
    
    /*
     * Begin main iteration loop
     */
    for(;;)
    {
        
        /*
         * Check for convergence or exceeding iteration count
         */
        if( m<=1 )
        {
            break;
        }
        if( iter>maxit )
        {
            result = ae_false;
            ae_frame_leave(_state);
            return result;
        }
        
        /*
         * Find diagonal block of matrix to work on
         */
        if( ae_fp_less(tol,0)&&ae_fp_less_eq(ae_fabs(d->ptr.p_double[m], _state),thresh) )
        {
            d->ptr.p_double[m] = 0;
        }
        smax = ae_fabs(d->ptr.p_double[m], _state);
        smin = smax;
        matrixsplitflag = ae_false;
        for(lll=1; lll<=m-1; lll++)
        {
            ll = m-lll;
            abss = ae_fabs(d->ptr.p_double[ll], _state);
            abse = ae_fabs(e->ptr.p_double[ll], _state);
            if( ae_fp_less(tol,0)&&ae_fp_less_eq(abss,thresh) )
            {
                d->ptr.p_double[ll] = 0;
            }
            if( ae_fp_less_eq(abse,thresh) )
            {
                matrixsplitflag = ae_true;
                break;
            }
            smin = ae_minreal(smin, abss, _state);
            smax = ae_maxreal(smax, ae_maxreal(abss, abse, _state), _state);
        }
        if( !matrixsplitflag )
        {
            ll = 0;
        }
        else
        {
            
            /*
             * Matrix splits since E(LL) = 0
             */
            e->ptr.p_double[ll] = 0;
            if( ll==m-1 )
            {
                
                /*
                 * Convergence of bottom singular value, return to top of loop
                 */
                m = m-1;
                continue;
            }
        }
        ll = ll+1;
        
        /*
         * E(LL) through E(M-1) are nonzero, E(LL-1) is zero
         */
        if( ll==m-1 )
        {
            
            /*
             * 2 by 2 block, handle separately
             */
            bdsvd_svdv2x2(d->ptr.p_double[m-1], e->ptr.p_double[m-1], d->ptr.p_double[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl, _state);
            d->ptr.p_double[m-1] = sigmx;
            e->ptr.p_double[m-1] = 0;
            d->ptr.p_double[m] = sigmn;
            
            /*
             * Compute singular vectors, if desired
             */
            if( ncvt>0 )
            {
                mm0 = m+(vstart-1);
                mm1 = m-1+(vstart-1);
                ae_v_moved(&vttemp.ptr.p_double[vstart], 1, &vt->ptr.pp_double[mm1][vstart], 1, ae_v_len(vstart,vend), cosr);
                ae_v_addd(&vttemp.ptr.p_double[vstart], 1, &vt->ptr.pp_double[mm0][vstart], 1, ae_v_len(vstart,vend), sinr);
                ae_v_muld(&vt->ptr.pp_double[mm0][vstart], 1, ae_v_len(vstart,vend), cosr);
                ae_v_subd(&vt->ptr.pp_double[mm0][vstart], 1, &vt->ptr.pp_double[mm1][vstart], 1, ae_v_len(vstart,vend), sinr);
                ae_v_move(&vt->ptr.pp_double[mm1][vstart], 1, &vttemp.ptr.p_double[vstart], 1, ae_v_len(vstart,vend));
            }
            if( nru>0 )
            {
                mm0 = m+ustart-1;
                mm1 = m-1+ustart-1;
                ae_v_moved(&utemp.ptr.p_double[ustart], 1, &u->ptr.pp_double[ustart][mm1], u->stride, ae_v_len(ustart,uend), cosl);
                ae_v_addd(&utemp.ptr.p_double[ustart], 1, &u->ptr.pp_double[ustart][mm0], u->stride, ae_v_len(ustart,uend), sinl);
                ae_v_muld(&u->ptr.pp_double[ustart][mm0], u->stride, ae_v_len(ustart,uend), cosl);
                ae_v_subd(&u->ptr.pp_double[ustart][mm0], u->stride, &u->ptr.pp_double[ustart][mm1], u->stride, ae_v_len(ustart,uend), sinl);
                ae_v_move(&u->ptr.pp_double[ustart][mm1], u->stride, &utemp.ptr.p_double[ustart], 1, ae_v_len(ustart,uend));
            }
            if( ncc>0 )
            {
                mm0 = m+cstart-1;
                mm1 = m-1+cstart-1;
                ae_v_moved(&ctemp.ptr.p_double[cstart], 1, &c->ptr.pp_double[mm1][cstart], 1, ae_v_len(cstart,cend), cosl);
                ae_v_addd(&ctemp.ptr.p_double[cstart], 1, &c->ptr.pp_double[mm0][cstart], 1, ae_v_len(cstart,cend), sinl);
                ae_v_muld(&c->ptr.pp_double[mm0][cstart], 1, ae_v_len(cstart,cend), cosl);
                ae_v_subd(&c->ptr.pp_double[mm0][cstart], 1, &c->ptr.pp_double[mm1][cstart], 1, ae_v_len(cstart,cend), sinl);
                ae_v_move(&c->ptr.pp_double[mm1][cstart], 1, &ctemp.ptr.p_double[cstart], 1, ae_v_len(cstart,cend));
            }
            m = m-2;
            continue;
        }
        
        /*
         * If working on new submatrix, choose shift direction
         * (from larger end diagonal element towards smaller)
         *
         * Previously was
         *     "if (LL>OLDM) or (M<OLDLL) then"
         * fixed thanks to Michael Rolle < m@rolle.name >
         * Very strange that LAPACK still contains it.
         */
        bchangedir = ae_false;
        if( idir==1&&ae_fp_less(ae_fabs(d->ptr.p_double[ll], _state),1.0E-3*ae_fabs(d->ptr.p_double[m], _state)) )
        {
            bchangedir = ae_true;
        }
        if( idir==2&&ae_fp_less(ae_fabs(d->ptr.p_double[m], _state),1.0E-3*ae_fabs(d->ptr.p_double[ll], _state)) )
        {
            bchangedir = ae_true;
        }
        if( (ll!=oldll||m!=oldm)||bchangedir )
        {
            if( ae_fp_greater_eq(ae_fabs(d->ptr.p_double[ll], _state),ae_fabs(d->ptr.p_double[m], _state)) )
            {
                
                /*
                 * Chase bulge from top (big end) to bottom (small end)
                 */
                idir = 1;
            }
            else
            {
                
                /*
                 * Chase bulge from bottom (big end) to top (small end)
                 */
                idir = 2;
            }
        }
        
        /*
         * Apply convergence tests
         */
        if( idir==1 )
        {
            
            /*
             * Run convergence test in forward direction
             * First apply standard test to bottom of matrix
             */
            if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[m-1], _state),ae_fabs(tol, _state)*ae_fabs(d->ptr.p_double[m], _state))||(ae_fp_less(tol,0)&&ae_fp_less_eq(ae_fabs(e->ptr.p_double[m-1], _state),thresh)) )
            {
                e->ptr.p_double[m-1] = 0;
                continue;
            }
            if( ae_fp_greater_eq(tol,0) )
            {
                
                /*
                 * If relative accuracy desired,
                 * apply convergence criterion forward
                 */
                mu = ae_fabs(d->ptr.p_double[ll], _state);
                sminl = mu;
                iterflag = ae_false;
                for(lll=ll; lll<=m-1; lll++)
                {
                    if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[lll], _state),tol*mu) )
                    {
                        e->ptr.p_double[lll] = 0;
                        iterflag = ae_true;
                        break;
                    }
                    sminlo = sminl;
                    mu = ae_fabs(d->ptr.p_double[lll+1], _state)*(mu/(mu+ae_fabs(e->ptr.p_double[lll], _state)));
                    sminl = ae_minreal(sminl, mu, _state);
                }
                if( iterflag )
                {
                    continue;
                }
            }
        }
        else
        {
            
            /*
             * Run convergence test in backward direction
             * First apply standard test to top of matrix
             */
            if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[ll], _state),ae_fabs(tol, _state)*ae_fabs(d->ptr.p_double[ll], _state))||(ae_fp_less(tol,0)&&ae_fp_less_eq(ae_fabs(e->ptr.p_double[ll], _state),thresh)) )
            {
                e->ptr.p_double[ll] = 0;
                continue;
            }
            if( ae_fp_greater_eq(tol,0) )
            {
                
                /*
                 * If relative accuracy desired,
                 * apply convergence criterion backward
                 */
                mu = ae_fabs(d->ptr.p_double[m], _state);
                sminl = mu;
                iterflag = ae_false;
                for(lll=m-1; lll>=ll; lll--)
                {
                    if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[lll], _state),tol*mu) )
                    {
                        e->ptr.p_double[lll] = 0;
                        iterflag = ae_true;
                        break;
                    }
                    sminlo = sminl;
                    mu = ae_fabs(d->ptr.p_double[lll], _state)*(mu/(mu+ae_fabs(e->ptr.p_double[lll], _state)));
                    sminl = ae_minreal(sminl, mu, _state);
                }
                if( iterflag )
                {
                    continue;
                }
            }
        }
        oldll = ll;
        oldm = m;
        
        /*
         * Compute shift.  First, test if shifting would ruin relative
         * accuracy, and if so set the shift to zero.
         */
        if( ae_fp_greater_eq(tol,0)&&ae_fp_less_eq(n*tol*(sminl/smax),ae_maxreal(eps, 0.01*tol, _state)) )
        {
            
            /*
             * Use a zero shift to avoid loss of relative accuracy
             */
            shift = 0;
        }
        else
        {
            
            /*
             * Compute the shift from 2-by-2 block at end of matrix
             */
            if( idir==1 )
            {
                sll = ae_fabs(d->ptr.p_double[ll], _state);
                bdsvd_svd2x2(d->ptr.p_double[m-1], e->ptr.p_double[m-1], d->ptr.p_double[m], &shift, &r, _state);
            }
            else
            {
                sll = ae_fabs(d->ptr.p_double[m], _state);
                bdsvd_svd2x2(d->ptr.p_double[ll], e->ptr.p_double[ll], d->ptr.p_double[ll+1], &shift, &r, _state);
            }
            
            /*
             * Test if shift negligible, and if so set to zero
             */
            if( ae_fp_greater(sll,0) )
            {
                if( ae_fp_less(ae_sqr(shift/sll, _state),eps) )
                {
                    shift = 0;
                }
            }
        }
        
        /*
         * Increment iteration count
         */
        iter = iter+m-ll;
        
        /*
         * If SHIFT = 0, do simplified QR iteration
         */
        if( ae_fp_eq(shift,0) )
        {
            if( idir==1 )
            {
                
                /*
                 * Chase bulge from top to bottom
                 * Save cosines and sines for later singular vector updates
                 */
                cs = 1;
                oldcs = 1;
                for(i=ll; i<=m-1; i++)
                {
                    generaterotation(d->ptr.p_double[i]*cs, e->ptr.p_double[i], &cs, &sn, &r, _state);
                    if( i>ll )
                    {
                        e->ptr.p_double[i-1] = oldsn*r;
                    }
                    generaterotation(oldcs*r, d->ptr.p_double[i+1]*sn, &oldcs, &oldsn, &tmp, _state);
                    d->ptr.p_double[i] = tmp;
                    work0.ptr.p_double[i-ll+1] = cs;
                    work1.ptr.p_double[i-ll+1] = sn;
                    work2.ptr.p_double[i-ll+1] = oldcs;
                    work3.ptr.p_double[i-ll+1] = oldsn;
                }
                h = d->ptr.p_double[m]*cs;
                d->ptr.p_double[m] = h*oldcs;
                e->ptr.p_double[m-1] = h*oldsn;
                
                /*
                 * Update singular vectors
                 */
                if( ncvt>0 )
                {
                    applyrotationsfromtheleft(fwddir, ll+vstart-1, m+vstart-1, vstart, vend, &work0, &work1, vt, &vttemp, _state);
                }
                if( nru>0 )
                {
                    applyrotationsfromtheright(fwddir, ustart, uend, ll+ustart-1, m+ustart-1, &work2, &work3, u, &utemp, _state);
                }
                if( ncc>0 )
                {
                    applyrotationsfromtheleft(fwddir, ll+cstart-1, m+cstart-1, cstart, cend, &work2, &work3, c, &ctemp, _state);
                }
                
                /*
                 * Test convergence
                 */
                if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[m-1], _state),thresh) )
                {
                    e->ptr.p_double[m-1] = 0;
                }
            }
            else
            {
                
                /*
                 * Chase bulge from bottom to top
                 * Save cosines and sines for later singular vector updates
                 */
                cs = 1;
                oldcs = 1;
                for(i=m; i>=ll+1; i--)
                {
                    generaterotation(d->ptr.p_double[i]*cs, e->ptr.p_double[i-1], &cs, &sn, &r, _state);
                    if( i<m )
                    {
                        e->ptr.p_double[i] = oldsn*r;
                    }
                    generaterotation(oldcs*r, d->ptr.p_double[i-1]*sn, &oldcs, &oldsn, &tmp, _state);
                    d->ptr.p_double[i] = tmp;
                    work0.ptr.p_double[i-ll] = cs;
                    work1.ptr.p_double[i-ll] = -sn;
                    work2.ptr.p_double[i-ll] = oldcs;
                    work3.ptr.p_double[i-ll] = -oldsn;
                }
                h = d->ptr.p_double[ll]*cs;
                d->ptr.p_double[ll] = h*oldcs;
                e->ptr.p_double[ll] = h*oldsn;
                
                /*
                 * Update singular vectors
                 */
                if( ncvt>0 )
                {
                    applyrotationsfromtheleft(!fwddir, ll+vstart-1, m+vstart-1, vstart, vend, &work2, &work3, vt, &vttemp, _state);
                }
                if( nru>0 )
                {
                    applyrotationsfromtheright(!fwddir, ustart, uend, ll+ustart-1, m+ustart-1, &work0, &work1, u, &utemp, _state);
                }
                if( ncc>0 )
                {
                    applyrotationsfromtheleft(!fwddir, ll+cstart-1, m+cstart-1, cstart, cend, &work0, &work1, c, &ctemp, _state);
                }
                
                /*
                 * Test convergence
                 */
                if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[ll], _state),thresh) )
                {
                    e->ptr.p_double[ll] = 0;
                }
            }
        }
        else
        {
            
            /*
             * Use nonzero shift
             */
            if( idir==1 )
            {
                
                /*
                 * Chase bulge from top to bottom
                 * Save cosines and sines for later singular vector updates
                 */
                f = (ae_fabs(d->ptr.p_double[ll], _state)-shift)*(bdsvd_extsignbdsqr(1, d->ptr.p_double[ll], _state)+shift/d->ptr.p_double[ll]);
                g = e->ptr.p_double[ll];
                for(i=ll; i<=m-1; i++)
                {
                    generaterotation(f, g, &cosr, &sinr, &r, _state);
                    if( i>ll )
                    {
                        e->ptr.p_double[i-1] = r;
                    }
                    f = cosr*d->ptr.p_double[i]+sinr*e->ptr.p_double[i];
                    e->ptr.p_double[i] = cosr*e->ptr.p_double[i]-sinr*d->ptr.p_double[i];
                    g = sinr*d->ptr.p_double[i+1];
                    d->ptr.p_double[i+1] = cosr*d->ptr.p_double[i+1];
                    generaterotation(f, g, &cosl, &sinl, &r, _state);
                    d->ptr.p_double[i] = r;
                    f = cosl*e->ptr.p_double[i]+sinl*d->ptr.p_double[i+1];
                    d->ptr.p_double[i+1] = cosl*d->ptr.p_double[i+1]-sinl*e->ptr.p_double[i];
                    if( i<m-1 )
                    {
                        g = sinl*e->ptr.p_double[i+1];
                        e->ptr.p_double[i+1] = cosl*e->ptr.p_double[i+1];
                    }
                    work0.ptr.p_double[i-ll+1] = cosr;
                    work1.ptr.p_double[i-ll+1] = sinr;
                    work2.ptr.p_double[i-ll+1] = cosl;
                    work3.ptr.p_double[i-ll+1] = sinl;
                }
                e->ptr.p_double[m-1] = f;
                
                /*
                 * Update singular vectors
                 */
                if( ncvt>0 )
                {
                    applyrotationsfromtheleft(fwddir, ll+vstart-1, m+vstart-1, vstart, vend, &work0, &work1, vt, &vttemp, _state);
                }
                if( nru>0 )
                {
                    applyrotationsfromtheright(fwddir, ustart, uend, ll+ustart-1, m+ustart-1, &work2, &work3, u, &utemp, _state);
                }
                if( ncc>0 )
                {
                    applyrotationsfromtheleft(fwddir, ll+cstart-1, m+cstart-1, cstart, cend, &work2, &work3, c, &ctemp, _state);
                }
                
                /*
                 * Test convergence
                 */
                if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[m-1], _state),thresh) )
                {
                    e->ptr.p_double[m-1] = 0;
                }
            }
            else
            {
                
                /*
                 * Chase bulge from bottom to top
                 * Save cosines and sines for later singular vector updates
                 */
                f = (ae_fabs(d->ptr.p_double[m], _state)-shift)*(bdsvd_extsignbdsqr(1, d->ptr.p_double[m], _state)+shift/d->ptr.p_double[m]);
                g = e->ptr.p_double[m-1];
                for(i=m; i>=ll+1; i--)
                {
                    generaterotation(f, g, &cosr, &sinr, &r, _state);
                    if( i<m )
                    {
                        e->ptr.p_double[i] = r;
                    }
                    f = cosr*d->ptr.p_double[i]+sinr*e->ptr.p_double[i-1];
                    e->ptr.p_double[i-1] = cosr*e->ptr.p_double[i-1]-sinr*d->ptr.p_double[i];
                    g = sinr*d->ptr.p_double[i-1];
                    d->ptr.p_double[i-1] = cosr*d->ptr.p_double[i-1];
                    generaterotation(f, g, &cosl, &sinl, &r, _state);
                    d->ptr.p_double[i] = r;
                    f = cosl*e->ptr.p_double[i-1]+sinl*d->ptr.p_double[i-1];
                    d->ptr.p_double[i-1] = cosl*d->ptr.p_double[i-1]-sinl*e->ptr.p_double[i-1];
                    if( i>ll+1 )
                    {
                        g = sinl*e->ptr.p_double[i-2];
                        e->ptr.p_double[i-2] = cosl*e->ptr.p_double[i-2];
                    }
                    work0.ptr.p_double[i-ll] = cosr;
                    work1.ptr.p_double[i-ll] = -sinr;
                    work2.ptr.p_double[i-ll] = cosl;
                    work3.ptr.p_double[i-ll] = -sinl;
                }
                e->ptr.p_double[ll] = f;
                
                /*
                 * Test convergence
                 */
                if( ae_fp_less_eq(ae_fabs(e->ptr.p_double[ll], _state),thresh) )
                {
                    e->ptr.p_double[ll] = 0;
                }
                
                /*
                 * Update singular vectors if desired
                 */
                if( ncvt>0 )
                {
                    applyrotationsfromtheleft(!fwddir, ll+vstart-1, m+vstart-1, vstart, vend, &work2, &work3, vt, &vttemp, _state);
                }
                if( nru>0 )
                {
                    applyrotationsfromtheright(!fwddir, ustart, uend, ll+ustart-1, m+ustart-1, &work0, &work1, u, &utemp, _state);
                }
                if( ncc>0 )
                {
                    applyrotationsfromtheleft(!fwddir, ll+cstart-1, m+cstart-1, cstart, cend, &work0, &work1, c, &ctemp, _state);
                }
            }
        }
        
        /*
         * QR iteration finished, go back and check convergence
         */
        continue;
    }
    
    /*
     * All singular values converged, so make them positive
     */
    for(i=1; i<=n; i++)
    {
        if( ae_fp_less(d->ptr.p_double[i],0) )
        {
            d->ptr.p_double[i] = -d->ptr.p_double[i];
            
            /*
             * Change sign of singular vectors, if desired
             */
            if( ncvt>0 )
            {
                ae_v_muld(&vt->ptr.pp_double[i+vstart-1][vstart], 1, ae_v_len(vstart,vend), -1);
            }
        }
    }
    
    /*
     * Sort the singular values into decreasing order (insertion sort on
     * singular values, but only one transposition per singular vector)
     */
    for(i=1; i<=n-1; i++)
    {
        
        /*
         * Scan for smallest D(I)
         */
        isub = 1;
        smin = d->ptr.p_double[1];
        for(j=2; j<=n+1-i; j++)
        {
            if( ae_fp_less_eq(d->ptr.p_double[j],smin) )
            {
                isub = j;
                smin = d->ptr.p_double[j];
            }
        }
        if( isub!=n+1-i )
        {
            
            /*
             * Swap singular values and vectors
             */
            d->ptr.p_double[isub] = d->ptr.p_double[n+1-i];
            d->ptr.p_double[n+1-i] = smin;
            if( ncvt>0 )
            {
                j = n+1-i;
                ae_v_move(&vttemp.ptr.p_double[vstart], 1, &vt->ptr.pp_double[isub+vstart-1][vstart], 1, ae_v_len(vstart,vend));
                ae_v_move(&vt->ptr.pp_double[isub+vstart-1][vstart], 1, &vt->ptr.pp_double[j+vstart-1][vstart], 1, ae_v_len(vstart,vend));
                ae_v_move(&vt->ptr.pp_double[j+vstart-1][vstart], 1, &vttemp.ptr.p_double[vstart], 1, ae_v_len(vstart,vend));
            }
            if( nru>0 )
            {
                j = n+1-i;
                ae_v_move(&utemp.ptr.p_double[ustart], 1, &u->ptr.pp_double[ustart][isub+ustart-1], u->stride, ae_v_len(ustart,uend));
                ae_v_move(&u->ptr.pp_double[ustart][isub+ustart-1], u->stride, &u->ptr.pp_double[ustart][j+ustart-1], u->stride, ae_v_len(ustart,uend));
                ae_v_move(&u->ptr.pp_double[ustart][j+ustart-1], u->stride, &utemp.ptr.p_double[ustart], 1, ae_v_len(ustart,uend));
            }
            if( ncc>0 )
            {
                j = n+1-i;
                ae_v_move(&ctemp.ptr.p_double[cstart], 1, &c->ptr.pp_double[isub+cstart-1][cstart], 1, ae_v_len(cstart,cend));
                ae_v_move(&c->ptr.pp_double[isub+cstart-1][cstart], 1, &c->ptr.pp_double[j+cstart-1][cstart], 1, ae_v_len(cstart,cend));
                ae_v_move(&c->ptr.pp_double[j+cstart-1][cstart], 1, &ctemp.ptr.p_double[cstart], 1, ae_v_len(cstart,cend));
            }
        }
    }
    ae_frame_leave(_state);
    return result;
}


static double bdsvd_extsignbdsqr(double a, double b, ae_state *_state)
{
    double result;


    if( ae_fp_greater_eq(b,0) )
    {
        result = ae_fabs(a, _state);
    }
    else
    {
        result = -ae_fabs(a, _state);
    }
    return result;
}


static void bdsvd_svd2x2(double f,
     double g,
     double h,
     double* ssmin,
     double* ssmax,
     ae_state *_state)
{
    double aas;
    double at;
    double au;
    double c;
    double fa;
    double fhmn;
    double fhmx;
    double ga;
    double ha;

    *ssmin = 0;
    *ssmax = 0;

    fa = ae_fabs(f, _state);
    ga = ae_fabs(g, _state);
    ha = ae_fabs(h, _state);
    fhmn = ae_minreal(fa, ha, _state);
    fhmx = ae_maxreal(fa, ha, _state);
    if( ae_fp_eq(fhmn,0) )
    {
        *ssmin = 0;
        if( ae_fp_eq(fhmx,0) )
        {
            *ssmax = ga;
        }
        else
        {
            *ssmax = ae_maxreal(fhmx, ga, _state)*ae_sqrt(1+ae_sqr(ae_minreal(fhmx, ga, _state)/ae_maxreal(fhmx, ga, _state), _state), _state);
        }
    }
    else
    {
        if( ae_fp_less(ga,fhmx) )
        {
            aas = 1+fhmn/fhmx;
            at = (fhmx-fhmn)/fhmx;
            au = ae_sqr(ga/fhmx, _state);
            c = 2/(ae_sqrt(aas*aas+au, _state)+ae_sqrt(at*at+au, _state));
            *ssmin = fhmn*c;
            *ssmax = fhmx/c;
        }
        else
        {
            au = fhmx/ga;
            if( ae_fp_eq(au,0) )
            {
                
                /*
                 * Avoid possible harmful underflow if exponent range
                 * asymmetric (true SSMIN may not underflow even if
                 * AU underflows)
                 */
                *ssmin = fhmn*fhmx/ga;
                *ssmax = ga;
            }
            else
            {
                aas = 1+fhmn/fhmx;
                at = (fhmx-fhmn)/fhmx;
                c = 1/(ae_sqrt(1+ae_sqr(aas*au, _state), _state)+ae_sqrt(1+ae_sqr(at*au, _state), _state));
                *ssmin = fhmn*c*au;
                *ssmin = *ssmin+(*ssmin);
                *ssmax = ga/(c+c);
            }
        }
    }
}


static void bdsvd_svdv2x2(double f,
     double g,
     double h,
     double* ssmin,
     double* ssmax,
     double* snr,
     double* csr,
     double* snl,
     double* csl,
     ae_state *_state)
{
    ae_bool gasmal;
    ae_bool swp;
    ae_int_t pmax;
    double a;
    double clt;
    double crt;
    double d;
    double fa;
    double ft;
    double ga;
    double gt;
    double ha;
    double ht;
    double l;
    double m;
    double mm;
    double r;
    double s;
    double slt;
    double srt;
    double t;
    double temp;
    double tsign;
    double tt;
    double v;

    *ssmin = 0;
    *ssmax = 0;
    *snr = 0;
    *csr = 0;
    *snl = 0;
    *csl = 0;

    ft = f;
    fa = ae_fabs(ft, _state);
    ht = h;
    ha = ae_fabs(h, _state);
    
    /*
     * these initializers are not really necessary,
     * but without them compiler complains about uninitialized locals
     */
    clt = 0;
    crt = 0;
    slt = 0;
    srt = 0;
    tsign = 0;
    
    /*
     * PMAX points to the maximum absolute element of matrix
     *  PMAX = 1 if F largest in absolute values
     *  PMAX = 2 if G largest in absolute values
     *  PMAX = 3 if H largest in absolute values
     */
    pmax = 1;
    swp = ae_fp_greater(ha,fa);
    if( swp )
    {
        
        /*
         * Now FA .ge. HA
         */
        pmax = 3;
        temp = ft;
        ft = ht;
        ht = temp;
        temp = fa;
        fa = ha;
        ha = temp;
    }
    gt = g;
    ga = ae_fabs(gt, _state);
    if( ae_fp_eq(ga,0) )
    {
        
        /*
         * Diagonal matrix
         */
        *ssmin = ha;
        *ssmax = fa;
        clt = 1;
        crt = 1;
        slt = 0;
        srt = 0;
    }
    else
    {
        gasmal = ae_true;
        if( ae_fp_greater(ga,fa) )
        {
            pmax = 2;
            if( ae_fp_less(fa/ga,ae_machineepsilon) )
            {
                
                /*
                 * Case of very large GA
                 */
                gasmal = ae_false;
                *ssmax = ga;
                if( ae_fp_greater(ha,1) )
                {
                    v = ga/ha;
                    *ssmin = fa/v;
                }
                else
                {
                    v = fa/ga;
                    *ssmin = v*ha;
                }
                clt = 1;
                slt = ht/gt;
                srt = 1;
                crt = ft/gt;
            }
        }
        if( gasmal )
        {
            
            /*
             * Normal case
             */
            d = fa-ha;
            if( ae_fp_eq(d,fa) )
            {
                l = 1;
            }
            else
            {
                l = d/fa;
            }
            m = gt/ft;
            t = 2-l;
            mm = m*m;
            tt = t*t;
            s = ae_sqrt(tt+mm, _state);
            if( ae_fp_eq(l,0) )
            {
                r = ae_fabs(m, _state);
            }
            else
            {
                r = ae_sqrt(l*l+mm, _state);
            }
            a = 0.5*(s+r);
            *ssmin = ha/a;
            *ssmax = fa*a;
            if( ae_fp_eq(mm,0) )
            {
                
                /*
                 * Note that M is very tiny
                 */
                if( ae_fp_eq(l,0) )
                {
                    t = bdsvd_extsignbdsqr(2, ft, _state)*bdsvd_extsignbdsqr(1, gt, _state);
                }
                else
                {
                    t = gt/bdsvd_extsignbdsqr(d, ft, _state)+m/t;
                }
            }
            else
            {
                t = (m/(s+t)+m/(r+l))*(1+a);
            }
            l = ae_sqrt(t*t+4, _state);
            crt = 2/l;
            srt = t/l;
            clt = (crt+srt*m)/a;
            v = ht/ft;
            slt = v*srt/a;
        }
    }
    if( swp )
    {
        *csl = srt;
        *snl = crt;
        *csr = slt;
        *snr = clt;
    }
    else
    {
        *csl = clt;
        *snl = slt;
        *csr = crt;
        *snr = srt;
    }
    
    /*
     * Correct signs of SSMAX and SSMIN
     */
    if( pmax==1 )
    {
        tsign = bdsvd_extsignbdsqr(1, *csr, _state)*bdsvd_extsignbdsqr(1, *csl, _state)*bdsvd_extsignbdsqr(1, f, _state);
    }
    if( pmax==2 )
    {
        tsign = bdsvd_extsignbdsqr(1, *snr, _state)*bdsvd_extsignbdsqr(1, *csl, _state)*bdsvd_extsignbdsqr(1, g, _state);
    }
    if( pmax==3 )
    {
        tsign = bdsvd_extsignbdsqr(1, *snr, _state)*bdsvd_extsignbdsqr(1, *snl, _state)*bdsvd_extsignbdsqr(1, h, _state);
    }
    *ssmax = bdsvd_extsignbdsqr(*ssmax, tsign, _state);
    *ssmin = bdsvd_extsignbdsqr(*ssmin, tsign*bdsvd_extsignbdsqr(1, f, _state)*bdsvd_extsignbdsqr(1, h, _state), _state);
}


/*$ End $*/
