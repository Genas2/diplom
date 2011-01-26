/*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

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
#include "lda.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherlda(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_vector* w,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix w2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(w);
    ae_matrix_init(&w2, 0, 0, DT_REAL, _state, ae_true);

    fisherldan(xy, npoints, nvars, nclasses, info, &w2, _state);
    if( *info>0 )
    {
        ae_vector_set_length(w, nvars-1+1, _state);
        ae_v_move(&w->ptr.p_double[0], 1, &w2.ptr.pp_double[0][0], w2.stride, ae_v_len(0,nvars-1));
    }
    ae_frame_leave(_state);
}


/*************************************************************************
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherldan(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_matrix* w,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t m;
    double v;
    ae_vector c;
    ae_vector mu;
    ae_matrix muc;
    ae_vector nc;
    ae_matrix sw;
    ae_matrix st;
    ae_matrix z;
    ae_matrix z2;
    ae_matrix tm;
    ae_matrix sbroot;
    ae_matrix a;
    ae_matrix xyproj;
    ae_matrix wproj;
    ae_vector tf;
    ae_vector d;
    ae_vector d2;
    ae_vector work;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_matrix_clear(w);
    ae_vector_init(&c, 0, DT_INT, _state, ae_true);
    ae_vector_init(&mu, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&muc, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&nc, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&sw, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&st, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z2, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&tm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&sbroot, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xyproj, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&wproj, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test data
     */
    if( (npoints<0||nvars<1)||nclasses<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( ae_round(xy->ptr.pp_double[i][nvars], _state)<0||ae_round(xy->ptr.pp_double[i][nvars], _state)>=nclasses )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Special case: NPoints<=1
     * Degenerate task.
     */
    if( npoints<=1 )
    {
        *info = 2;
        ae_matrix_set_length(w, nvars-1+1, nvars-1+1, _state);
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                if( i==j )
                {
                    w->ptr.pp_double[i][j] = 1;
                }
                else
                {
                    w->ptr.pp_double[i][j] = 0;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Prepare temporaries
     */
    ae_vector_set_length(&tf, nvars-1+1, _state);
    ae_vector_set_length(&work, ae_maxint(nvars, npoints, _state)+1, _state);
    
    /*
     * Convert class labels from reals to integers (just for convenience)
     */
    ae_vector_set_length(&c, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        c.ptr.p_int[i] = ae_round(xy->ptr.pp_double[i][nvars], _state);
    }
    
    /*
     * Calculate class sizes and means
     */
    ae_vector_set_length(&mu, nvars-1+1, _state);
    ae_matrix_set_length(&muc, nclasses-1+1, nvars-1+1, _state);
    ae_vector_set_length(&nc, nclasses-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        mu.ptr.p_double[j] = 0;
    }
    for(i=0; i<=nclasses-1; i++)
    {
        nc.ptr.p_int[i] = 0;
        for(j=0; j<=nvars-1; j++)
        {
            muc.ptr.pp_double[i][j] = 0;
        }
    }
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_add(&mu.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        ae_v_add(&muc.ptr.pp_double[c.ptr.p_int[i]][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        nc.ptr.p_int[c.ptr.p_int[i]] = nc.ptr.p_int[c.ptr.p_int[i]]+1;
    }
    for(i=0; i<=nclasses-1; i++)
    {
        v = (double)1/(double)nc.ptr.p_int[i];
        ae_v_muld(&muc.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1), v);
    }
    v = (double)1/(double)npoints;
    ae_v_muld(&mu.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
    
    /*
     * Create ST matrix
     */
    ae_matrix_set_length(&st, nvars-1+1, nvars-1+1, _state);
    for(i=0; i<=nvars-1; i++)
    {
        for(j=0; j<=nvars-1; j++)
        {
            st.ptr.pp_double[i][j] = 0;
        }
    }
    for(k=0; k<=npoints-1; k++)
    {
        ae_v_move(&tf.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nvars-1));
        ae_v_sub(&tf.ptr.p_double[0], 1, &mu.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
        for(i=0; i<=nvars-1; i++)
        {
            v = tf.ptr.p_double[i];
            ae_v_addd(&st.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
        }
    }
    
    /*
     * Create SW matrix
     */
    ae_matrix_set_length(&sw, nvars-1+1, nvars-1+1, _state);
    for(i=0; i<=nvars-1; i++)
    {
        for(j=0; j<=nvars-1; j++)
        {
            sw.ptr.pp_double[i][j] = 0;
        }
    }
    for(k=0; k<=npoints-1; k++)
    {
        ae_v_move(&tf.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nvars-1));
        ae_v_sub(&tf.ptr.p_double[0], 1, &muc.ptr.pp_double[c.ptr.p_int[k]][0], 1, ae_v_len(0,nvars-1));
        for(i=0; i<=nvars-1; i++)
        {
            v = tf.ptr.p_double[i];
            ae_v_addd(&sw.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
        }
    }
    
    /*
     * Maximize ratio J=(w'*ST*w)/(w'*SW*w).
     *
     * First, make transition from w to v such that w'*ST*w becomes v'*v:
     *    v  = root(ST)*w = R*w
     *    R  = root(D)*Z'
     *    w  = (root(ST)^-1)*v = RI*v
     *    RI = Z*inv(root(D))
     *    J  = (v'*v)/(v'*(RI'*SW*RI)*v)
     *    ST = Z*D*Z'
     *
     *    so we have
     *
     *    J = (v'*v) / (v'*(inv(root(D))*Z'*SW*Z*inv(root(D)))*v)  =
     *      = (v'*v) / (v'*A*v)
     */
    if( !smatrixevd(&st, nvars, 1, ae_true, &d, &z, _state) )
    {
        *info = -4;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(w, nvars-1+1, nvars-1+1, _state);
    if( ae_fp_less_eq(d.ptr.p_double[nvars-1],0)||ae_fp_less_eq(d.ptr.p_double[0],1000*ae_machineepsilon*d.ptr.p_double[nvars-1]) )
    {
        
        /*
         * Special case: D[NVars-1]<=0
         * Degenerate task (all variables takes the same value).
         */
        if( ae_fp_less_eq(d.ptr.p_double[nvars-1],0) )
        {
            *info = 2;
            for(i=0; i<=nvars-1; i++)
            {
                for(j=0; j<=nvars-1; j++)
                {
                    if( i==j )
                    {
                        w->ptr.pp_double[i][j] = 1;
                    }
                    else
                    {
                        w->ptr.pp_double[i][j] = 0;
                    }
                }
            }
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Special case: degenerate ST matrix, multicollinearity found.
         * Since we know ST eigenvalues/vectors we can translate task to
         * non-degenerate form.
         *
         * Let WG is orthogonal basis of the non zero variance subspace
         * of the ST and let WZ is orthogonal basis of the zero variance
         * subspace.
         *
         * Projection on WG allows us to use LDA on reduced M-dimensional
         * subspace, N-M vectors of WZ allows us to update reduced LDA
         * factors to full N-dimensional subspace.
         */
        m = 0;
        for(k=0; k<=nvars-1; k++)
        {
            if( ae_fp_less_eq(d.ptr.p_double[k],1000*ae_machineepsilon*d.ptr.p_double[nvars-1]) )
            {
                m = k+1;
            }
        }
        ae_assert(m!=0, "FisherLDAN: internal error #1", _state);
        ae_matrix_set_length(&xyproj, npoints-1+1, nvars-m+1, _state);
        matrixmatrixmultiply(xy, 0, npoints-1, 0, nvars-1, ae_false, &z, 0, nvars-1, m, nvars-1, ae_false, 1.0, &xyproj, 0, npoints-1, 0, nvars-m-1, 0.0, &work, _state);
        for(i=0; i<=npoints-1; i++)
        {
            xyproj.ptr.pp_double[i][nvars-m] = xy->ptr.pp_double[i][nvars];
        }
        fisherldan(&xyproj, npoints, nvars-m, nclasses, info, &wproj, _state);
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
        matrixmatrixmultiply(&z, 0, nvars-1, m, nvars-1, ae_false, &wproj, 0, nvars-m-1, 0, nvars-m-1, ae_false, 1.0, w, 0, nvars-1, 0, nvars-m-1, 0.0, &work, _state);
        for(k=nvars-m; k<=nvars-1; k++)
        {
            ae_v_move(&w->ptr.pp_double[0][k], w->stride, &z.ptr.pp_double[0][k-(nvars-m)], z.stride, ae_v_len(0,nvars-1));
        }
        *info = 2;
    }
    else
    {
        
        /*
         * General case: no multicollinearity
         */
        ae_matrix_set_length(&tm, nvars-1+1, nvars-1+1, _state);
        ae_matrix_set_length(&a, nvars-1+1, nvars-1+1, _state);
        matrixmatrixmultiply(&sw, 0, nvars-1, 0, nvars-1, ae_false, &z, 0, nvars-1, 0, nvars-1, ae_false, 1.0, &tm, 0, nvars-1, 0, nvars-1, 0.0, &work, _state);
        matrixmatrixmultiply(&z, 0, nvars-1, 0, nvars-1, ae_true, &tm, 0, nvars-1, 0, nvars-1, ae_false, 1.0, &a, 0, nvars-1, 0, nvars-1, 0.0, &work, _state);
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                a.ptr.pp_double[i][j] = a.ptr.pp_double[i][j]/ae_sqrt(d.ptr.p_double[i]*d.ptr.p_double[j], _state);
            }
        }
        if( !smatrixevd(&a, nvars, 1, ae_true, &d2, &z2, _state) )
        {
            *info = -4;
            ae_frame_leave(_state);
            return;
        }
        for(k=0; k<=nvars-1; k++)
        {
            for(i=0; i<=nvars-1; i++)
            {
                tf.ptr.p_double[i] = z2.ptr.pp_double[i][k]/ae_sqrt(d.ptr.p_double[i], _state);
            }
            for(i=0; i<=nvars-1; i++)
            {
                v = ae_v_dotproduct(&z.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                w->ptr.pp_double[i][k] = v;
            }
        }
    }
    
    /*
     * Post-processing:
     * * normalization
     * * converting to non-negative form, if possible
     */
    for(k=0; k<=nvars-1; k++)
    {
        v = ae_v_dotproduct(&w->ptr.pp_double[0][k], w->stride, &w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1));
        v = 1/ae_sqrt(v, _state);
        ae_v_muld(&w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1), v);
        v = 0;
        for(i=0; i<=nvars-1; i++)
        {
            v = v+w->ptr.pp_double[i][k];
        }
        if( ae_fp_less(v,0) )
        {
            ae_v_muld(&w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1), -1);
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
