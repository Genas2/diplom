

#include <stdafx.h>
#include <stdio.h>
#include "testpsplineunit.h"


/*$ Declarations $*/
static void testpsplineunit_unsetp2(pspline2interpolant* p,
     ae_state *_state);
static void testpsplineunit_unsetp3(pspline3interpolant* p,
     ae_state *_state);
static void testpsplineunit_unset1d(/* Real    */ ae_vector* x,
     ae_state *_state);


/*$ Body $*/


ae_bool testpspline(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool p2errors;
    ae_bool p3errors;
    double nonstrictthreshold;
    double threshold;
    ae_int_t passcount;
    double lstep;
    double h;
    ae_int_t maxn;
    ae_int_t periodicity;
    ae_int_t skind;
    ae_int_t pkind;
    ae_bool periodic;
    double a;
    double b;
    ae_int_t n;
    ae_int_t tmpn;
    ae_int_t i;
    ae_int_t k;
    double vx;
    double vy;
    double vz;
    double vx2;
    double vy2;
    double vz2;
    double vdx;
    double vdy;
    double vdz;
    double vdx2;
    double vdy2;
    double vdz2;
    double vd2x;
    double vd2y;
    double vd2z;
    double vd2x2;
    double vd2y2;
    double vd2z2;
    double v0;
    double v1;
    ae_vector x;
    ae_vector y;
    ae_vector z;
    ae_vector t;
    ae_vector t2;
    ae_vector t3;
    ae_matrix xy;
    ae_matrix xyz;
    pspline2interpolant p2;
    pspline3interpolant p3;
    spline1dinterpolant s;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&z, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&t, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&t2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&t3, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xyz, 0, 0, DT_REAL, _state, ae_true);
    _pspline2interpolant_init(&p2, _state, ae_true);
    _pspline3interpolant_init(&p3, _state, ae_true);
    _spline1dinterpolant_init(&s, _state, ae_true);

    waserrors = ae_false;
    passcount = 20;
    lstep = 0.005;
    h = 0.00001;
    maxn = 10;
    threshold = 10000*ae_machineepsilon;
    nonstrictthreshold = 0.00001;
    p2errors = ae_false;
    p3errors = ae_false;
    
    /*
     * Test basic properties of 2- and 3-dimensional splines:
     * * PSpline2ParameterValues() properties
     * * values at nodes
     * * for periodic splines - periodicity properties
     *
     * Variables used:
     * * N              points count
     * * SKind          spline
     * * PKind          parameterization
     * * Periodicity    whether we have periodic spline or not
     */
    for(n=2; n<=maxn; n++)
    {
        for(skind=0; skind<=2; skind++)
        {
            for(pkind=0; pkind<=2; pkind++)
            {
                for(periodicity=0; periodicity<=1; periodicity++)
                {
                    periodic = periodicity==1;
                    
                    /*
                     * skip unsupported combinations of parameters
                     */
                    if( periodic&&n<3 )
                    {
                        continue;
                    }
                    if( periodic&&skind==0 )
                    {
                        continue;
                    }
                    if( n<5&&skind==0 )
                    {
                        continue;
                    }
                    
                    /*
                     * init
                     */
                    ae_matrix_set_length(&xy, n, 2, _state);
                    ae_matrix_set_length(&xyz, n, 3, _state);
                    taskgenint1dequidist(-1, 1, n, &t2, &x, _state);
                    ae_v_move(&xy.ptr.pp_double[0][0], xy.stride, &x.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    ae_v_move(&xyz.ptr.pp_double[0][0], xyz.stride, &x.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    taskgenint1dequidist(-1, 1, n, &t2, &y, _state);
                    ae_v_move(&xy.ptr.pp_double[0][1], xy.stride, &y.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    ae_v_move(&xyz.ptr.pp_double[0][1], xyz.stride, &y.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    taskgenint1dequidist(-1, 1, n, &t2, &z, _state);
                    ae_v_move(&xyz.ptr.pp_double[0][2], xyz.stride, &z.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    testpsplineunit_unsetp2(&p2, _state);
                    testpsplineunit_unsetp3(&p3, _state);
                    if( periodic )
                    {
                        pspline2buildperiodic(&xy, n, skind, pkind, &p2, _state);
                        pspline3buildperiodic(&xyz, n, skind, pkind, &p3, _state);
                    }
                    else
                    {
                        pspline2build(&xy, n, skind, pkind, &p2, _state);
                        pspline3build(&xyz, n, skind, pkind, &p3, _state);
                    }
                    
                    /*
                     * PSpline2ParameterValues() properties
                     */
                    pspline2parametervalues(&p2, &tmpn, &t2, _state);
                    if( tmpn!=n )
                    {
                        p2errors = ae_true;
                        continue;
                    }
                    pspline3parametervalues(&p3, &tmpn, &t3, _state);
                    if( tmpn!=n )
                    {
                        p3errors = ae_true;
                        continue;
                    }
                    p2errors = p2errors||ae_fp_neq(t2.ptr.p_double[0],0);
                    p3errors = p3errors||ae_fp_neq(t3.ptr.p_double[0],0);
                    for(i=1; i<=n-1; i++)
                    {
                        p2errors = p2errors||ae_fp_less_eq(t2.ptr.p_double[i],t2.ptr.p_double[i-1]);
                        p3errors = p3errors||ae_fp_less_eq(t3.ptr.p_double[i],t3.ptr.p_double[i-1]);
                    }
                    if( periodic )
                    {
                        p2errors = p2errors||ae_fp_greater_eq(t2.ptr.p_double[n-1],1);
                        p3errors = p3errors||ae_fp_greater_eq(t3.ptr.p_double[n-1],1);
                    }
                    else
                    {
                        p2errors = p2errors||ae_fp_neq(t2.ptr.p_double[n-1],1);
                        p3errors = p3errors||ae_fp_neq(t3.ptr.p_double[n-1],1);
                    }
                    
                    /*
                     * Now we have parameter values stored at T,
                     * and want to test whether the actully correspond to
                     * points
                     */
                    for(i=0; i<=n-1; i++)
                    {
                        
                        /*
                         * 2-dimensional test
                         */
                        pspline2calc(&p2, t2.ptr.p_double[i], &vx, &vy, _state);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                        
                        /*
                         * 3-dimensional test
                         */
                        pspline3calc(&p3, t3.ptr.p_double[i], &vx, &vy, &vz, _state);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-z.ptr.p_double[i], _state),threshold);
                    }
                    
                    /*
                     * Test periodicity (if needed)
                     */
                    if( periodic )
                    {
                        
                        /*
                         * periodicity at nodes
                         */
                        for(i=0; i<=n-1; i++)
                        {
                            
                            /*
                             * 2-dimensional test
                             */
                            pspline2calc(&p2, t2.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vy, _state);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            pspline2diff(&p2, t2.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vdx, &vy, &vdy, _state);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            pspline2diff2(&p2, t2.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, _state);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            
                            /*
                             * 3-dimensional test
                             */
                            pspline3calc(&p3, t3.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vy, &vz, _state);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vz-z.ptr.p_double[i], _state),threshold);
                            pspline3diff(&p3, t3.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vdx, &vy, &vdy, &vz, &vdz, _state);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vz-z.ptr.p_double[i], _state),threshold);
                            pspline3diff2(&p3, t3.ptr.p_double[i]+ae_randominteger(10, _state)-5, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, &vz, &vdz, &vd2z, _state);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vx-x.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vy-y.ptr.p_double[i], _state),threshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vz-z.ptr.p_double[i], _state),threshold);
                        }
                        
                        /*
                         * periodicity between nodes
                         */
                        v0 = ae_randomreal(_state);
                        pspline2calc(&p2, v0, &vx, &vy, _state);
                        pspline2calc(&p2, v0+ae_randominteger(10, _state)-5, &vx2, &vy2, _state);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        pspline3calc(&p3, v0, &vx, &vy, &vz, _state);
                        pspline3calc(&p3, v0+ae_randominteger(10, _state)-5, &vx2, &vy2, &vz2, _state);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
                        
                        /*
                         * near-boundary test for continuity of function values and derivatives:
                         * 2-dimensional curve
                         */
                        ae_assert(skind==1||skind==2, "TEST: unexpected spline type!", _state);
                        v0 = 100*ae_machineepsilon;
                        v1 = 1-v0;
                        pspline2calc(&p2, v0, &vx, &vy, _state);
                        pspline2calc(&p2, v1, &vx2, &vy2, _state);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        pspline2diff(&p2, v0, &vx, &vdx, &vy, &vdy, _state);
                        pspline2diff(&p2, v1, &vx2, &vdx2, &vy2, &vdy2, _state);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),nonstrictthreshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),nonstrictthreshold);
                        pspline2diff2(&p2, v0, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, _state);
                        pspline2diff2(&p2, v1, &vx2, &vdx2, &vd2x2, &vy2, &vdy2, &vd2y2, _state);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),nonstrictthreshold);
                        p2errors = p2errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),nonstrictthreshold);
                        if( skind==2 )
                        {
                            
                            /*
                             * second derivative test only for cubic splines
                             */
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vd2x-vd2x2, _state),nonstrictthreshold);
                            p2errors = p2errors||ae_fp_greater(ae_fabs(vd2y-vd2y2, _state),nonstrictthreshold);
                        }
                        
                        /*
                         * near-boundary test for continuity of function values and derivatives:
                         * 3-dimensional curve
                         */
                        ae_assert(skind==1||skind==2, "TEST: unexpected spline type!", _state);
                        v0 = 100*ae_machineepsilon;
                        v1 = 1-v0;
                        pspline3calc(&p3, v0, &vx, &vy, &vz, _state);
                        pspline3calc(&p3, v1, &vx2, &vy2, &vz2, _state);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
                        pspline3diff(&p3, v0, &vx, &vdx, &vy, &vdy, &vz, &vdz, _state);
                        pspline3diff(&p3, v1, &vx2, &vdx2, &vy2, &vdy2, &vz2, &vdz2, _state);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),nonstrictthreshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),nonstrictthreshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdz-vdz2, _state),nonstrictthreshold);
                        pspline3diff2(&p3, v0, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, &vz, &vdz, &vd2z, _state);
                        pspline3diff2(&p3, v1, &vx2, &vdx2, &vd2x2, &vy2, &vdy2, &vd2y2, &vz2, &vdz2, &vd2z2, _state);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),nonstrictthreshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),nonstrictthreshold);
                        p3errors = p3errors||ae_fp_greater(ae_fabs(vdz-vdz2, _state),nonstrictthreshold);
                        if( skind==2 )
                        {
                            
                            /*
                             * second derivative test only for cubic splines
                             */
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vd2x-vd2x2, _state),nonstrictthreshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vd2y-vd2y2, _state),nonstrictthreshold);
                            p3errors = p3errors||ae_fp_greater(ae_fabs(vd2z-vd2z2, _state),nonstrictthreshold);
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Test differentiation, tangents, calculation between nodes.
     *
     * Because differentiation is done in parameterization/spline/periodicity
     * oblivious manner, we don't have to test all possible combinations
     * of spline types and parameterizations.
     *
     * Actually we test special combination with properties which allow us
     * to easily solve this problem:
     * * 2 (3) variables
     * * first variable is sampled from equidistant grid on [0,1]
     * * other variables are random
     * * uniform parameterization is used
     * * periodicity - none
     * * spline type - any (we use cubic splines)
     * Same problem allows us to test calculation BETWEEN nodes.
     */
    for(n=2; n<=maxn; n++)
    {
        
        /*
         * init
         */
        ae_matrix_set_length(&xy, n, 2, _state);
        ae_matrix_set_length(&xyz, n, 3, _state);
        taskgenint1dequidist(0, 1, n, &t, &x, _state);
        ae_v_move(&xy.ptr.pp_double[0][0], xy.stride, &x.ptr.p_double[0], 1, ae_v_len(0,n-1));
        ae_v_move(&xyz.ptr.pp_double[0][0], xyz.stride, &x.ptr.p_double[0], 1, ae_v_len(0,n-1));
        taskgenint1dequidist(0, 1, n, &t, &y, _state);
        ae_v_move(&xy.ptr.pp_double[0][1], xy.stride, &y.ptr.p_double[0], 1, ae_v_len(0,n-1));
        ae_v_move(&xyz.ptr.pp_double[0][1], xyz.stride, &y.ptr.p_double[0], 1, ae_v_len(0,n-1));
        taskgenint1dequidist(0, 1, n, &t, &z, _state);
        ae_v_move(&xyz.ptr.pp_double[0][2], xyz.stride, &z.ptr.p_double[0], 1, ae_v_len(0,n-1));
        testpsplineunit_unsetp2(&p2, _state);
        testpsplineunit_unsetp3(&p3, _state);
        pspline2build(&xy, n, 2, 0, &p2, _state);
        pspline3build(&xyz, n, 2, 0, &p3, _state);
        
        /*
         * Test 2D/3D spline:
         * * build non-parametric cubic spline from T and X/Y
         * * calculate its value and derivatives at V0
         * * compare with Spline2Calc/Spline2Diff/Spline2Diff2
         * Because of task properties both variants should
         * return same answer.
         */
        v0 = ae_randomreal(_state);
        spline1dbuildcubic(&t, &x, n, 0, 0.0, 0, 0.0, &s, _state);
        spline1ddiff(&s, v0, &vx2, &vdx2, &vd2x2, _state);
        spline1dbuildcubic(&t, &y, n, 0, 0.0, 0, 0.0, &s, _state);
        spline1ddiff(&s, v0, &vy2, &vdy2, &vd2y2, _state);
        spline1dbuildcubic(&t, &z, n, 0, 0.0, 0, 0.0, &s, _state);
        spline1ddiff(&s, v0, &vz2, &vdz2, &vd2z2, _state);
        
        /*
         * 2D test
         */
        pspline2calc(&p2, v0, &vx, &vy, _state);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        pspline2diff(&p2, v0, &vx, &vdx, &vy, &vdy, _state);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),threshold);
        pspline2diff2(&p2, v0, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, _state);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vd2x-vd2x2, _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vd2y-vd2y2, _state),threshold);
        
        /*
         * 3D test
         */
        pspline3calc(&p3, v0, &vx, &vy, &vz, _state);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
        pspline3diff(&p3, v0, &vx, &vdx, &vy, &vdy, &vz, &vdz, _state);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdz-vdz2, _state),threshold);
        pspline3diff2(&p3, v0, &vx, &vdx, &vd2x, &vy, &vdy, &vd2y, &vz, &vdz, &vd2z, _state);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vx2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vy2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vz2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdx-vdx2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdy-vdy2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vdz-vdz2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vd2x-vd2x2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vd2y-vd2y2, _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vd2z-vd2z2, _state),threshold);
        
        /*
         * Test tangents for 2D/3D
         */
        pspline2tangent(&p2, v0, &vx, &vy, _state);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vx-vdx2/safepythag2(vdx2, vdy2, _state), _state),threshold);
        p2errors = p2errors||ae_fp_greater(ae_fabs(vy-vdy2/safepythag2(vdx2, vdy2, _state), _state),threshold);
        pspline3tangent(&p3, v0, &vx, &vy, &vz, _state);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vx-vdx2/safepythag3(vdx2, vdy2, vdz2, _state), _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vy-vdy2/safepythag3(vdx2, vdy2, vdz2, _state), _state),threshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(vz-vdz2/safepythag3(vdx2, vdy2, vdz2, _state), _state),threshold);
    }
    
    /*
     * Arc length test.
     *
     * Simple problem with easy solution (points on a straight line with
     * uniform parameterization).
     */
    for(n=2; n<=maxn; n++)
    {
        ae_matrix_set_length(&xy, n, 2, _state);
        ae_matrix_set_length(&xyz, n, 3, _state);
        for(i=0; i<=n-1; i++)
        {
            xy.ptr.pp_double[i][0] = i;
            xy.ptr.pp_double[i][1] = i;
            xyz.ptr.pp_double[i][0] = i;
            xyz.ptr.pp_double[i][1] = i;
            xyz.ptr.pp_double[i][2] = i;
        }
        pspline2build(&xy, n, 1, 0, &p2, _state);
        pspline3build(&xyz, n, 1, 0, &p3, _state);
        a = ae_randomreal(_state);
        b = ae_randomreal(_state);
        p2errors = p2errors||ae_fp_greater(ae_fabs(pspline2arclength(&p2, a, b, _state)-(b-a)*ae_sqrt(2, _state)*(n-1), _state),nonstrictthreshold);
        p3errors = p3errors||ae_fp_greater(ae_fabs(pspline3arclength(&p3, a, b, _state)-(b-a)*ae_sqrt(3, _state)*(n-1), _state),nonstrictthreshold);
    }
    
    /*
     * report
     */
    waserrors = p2errors||p3errors;
    if( !silent )
    {
        printf("TESTING SPLINE INTERPOLATION\n");
        
        /*
         * Normal tests
         */
        printf("2D TEST:                                 ");
        if( p2errors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("3D TEST:                                 ");
        if( p3errors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    
    /*
     * end
     */
    result = !waserrors;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************/
static void testpsplineunit_unsetp2(pspline2interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xy;

    ae_frame_make(_state, &_frame_block);
    _pspline2interpolant_clear(p);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);

    ae_matrix_set_length(&xy, 2, 2, _state);
    xy.ptr.pp_double[0][0] = -1;
    xy.ptr.pp_double[0][1] = -1;
    xy.ptr.pp_double[1][0] = 1;
    xy.ptr.pp_double[1][1] = 1;
    pspline2build(&xy, 2, 1, 0, p, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************/
static void testpsplineunit_unsetp3(pspline3interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xy;

    ae_frame_make(_state, &_frame_block);
    _pspline3interpolant_clear(p);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);

    ae_matrix_set_length(&xy, 2, 3, _state);
    xy.ptr.pp_double[0][0] = -1;
    xy.ptr.pp_double[0][1] = -1;
    xy.ptr.pp_double[0][2] = -1;
    xy.ptr.pp_double[1][0] = 1;
    xy.ptr.pp_double[1][1] = 1;
    xy.ptr.pp_double[1][2] = 1;
    pspline3build(&xy, 2, 1, 0, p, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Unsets real vector
*************************************************************************/
static void testpsplineunit_unset1d(/* Real    */ ae_vector* x,
     ae_state *_state)
{


    ae_vector_set_length(x, 1, _state);
    x->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
}


/*$ End $*/
