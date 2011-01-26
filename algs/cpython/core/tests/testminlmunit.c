

#include <stdafx.h>
#include <stdio.h>
#include "testminlmunit.h"


/*$ Declarations $*/
static ae_bool testminlmunit_rkindvsstatecheck(ae_int_t rkind,
     minlmstate* state,
     ae_state *_state);
static void testminlmunit_axmb(minlmstate* state,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     ae_state *_state);


/*$ Body $*/


ae_bool testminlm(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool referror;
    ae_bool lin1error;
    ae_bool lin2error;
    ae_bool eqerror;
    ae_bool converror;
    ae_bool scerror;
    ae_bool restartserror;
    ae_bool othererrors;
    ae_int_t rkind;
    ae_int_t ckind;
    double epsf;
    double epsx;
    double epsg;
    ae_int_t maxits;
    ae_int_t n;
    ae_int_t m;
    ae_vector x;
    ae_vector xe;
    ae_vector b;
    ae_vector xlast;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double v;
    double s;
    double stpmax;
    double h;
    ae_matrix a;
    double fprev;
    double xprev;
    minlmstate state;
    minlmreport rep;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xe, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xlast, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    _minlmstate_init(&state, _state, ae_true);
    _minlmreport_init(&rep, _state, ae_true);

    waserrors = ae_false;
    referror = ae_false;
    lin1error = ae_false;
    lin2error = ae_false;
    eqerror = ae_false;
    converror = ae_false;
    scerror = ae_false;
    othererrors = ae_false;
    restartserror = ae_false;
    
    /*
     * Reference problem.
     * See comments for RKindVsStateCheck() for more info about RKind.
     *
     * NOTES: we also test negative RKind's corresponding to "inexact" schemes
     * which use approximate finite difference Jacobian.
     */
    ae_vector_set_length(&x, 3, _state);
    n = 3;
    m = 3;
    h = 0.0001;
    for(rkind=-2; rkind<=5; rkind++)
    {
        x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
        x.ptr.p_double[1] = 100*ae_randomreal(_state)-50;
        x.ptr.p_double[2] = 100*ae_randomreal(_state)-50;
        if( rkind==-2 )
        {
            minlmcreatev(n, m, &x, h, &state, _state);
            minlmsetacctype(&state, 1, _state);
        }
        if( rkind==-1 )
        {
            minlmcreatev(n, m, &x, h, &state, _state);
            minlmsetacctype(&state, 0, _state);
        }
        if( rkind==0 )
        {
            minlmcreatefj(n, m, &x, &state, _state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, &x, &state, _state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, &x, &state, _state);
        }
        if( rkind==3 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 0, _state);
        }
        if( rkind==4 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 1, _state);
        }
        if( rkind==5 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 2, _state);
        }
        while(minlmiteration(&state, _state))
        {
            
            /*
             * (x-2)^2 + y^2 + (z-x)^2
             */
            if( state.needfi )
            {
                state.fi.ptr.p_double[0] = state.x.ptr.p_double[0]-2;
                state.fi.ptr.p_double[1] = state.x.ptr.p_double[1];
                state.fi.ptr.p_double[2] = state.x.ptr.p_double[2]-state.x.ptr.p_double[0];
            }
            if( state.needfij )
            {
                state.fi.ptr.p_double[0] = state.x.ptr.p_double[0]-2;
                state.fi.ptr.p_double[1] = state.x.ptr.p_double[1];
                state.fi.ptr.p_double[2] = state.x.ptr.p_double[2]-state.x.ptr.p_double[0];
                state.j.ptr.pp_double[0][0] = 1;
                state.j.ptr.pp_double[0][1] = 0;
                state.j.ptr.pp_double[0][2] = 0;
                state.j.ptr.pp_double[1][0] = 0;
                state.j.ptr.pp_double[1][1] = 1;
                state.j.ptr.pp_double[1][2] = 0;
                state.j.ptr.pp_double[2][0] = -1;
                state.j.ptr.pp_double[2][1] = 0;
                state.j.ptr.pp_double[2][2] = 1;
            }
            if( (state.needf||state.needfg)||state.needfgh )
            {
                state.f = ae_sqr(state.x.ptr.p_double[0]-2, _state)+ae_sqr(state.x.ptr.p_double[1], _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
            }
            if( state.needfg||state.needfgh )
            {
                state.g.ptr.p_double[0] = 2*(state.x.ptr.p_double[0]-2)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
                state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
                state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
            }
            if( state.needfgh )
            {
                state.h.ptr.pp_double[0][0] = 4;
                state.h.ptr.pp_double[0][1] = 0;
                state.h.ptr.pp_double[0][2] = -2;
                state.h.ptr.pp_double[1][0] = 0;
                state.h.ptr.pp_double[1][1] = 2;
                state.h.ptr.pp_double[1][2] = 0;
                state.h.ptr.pp_double[2][0] = -2;
                state.h.ptr.pp_double[2][1] = 0;
                state.h.ptr.pp_double[2][2] = 2;
            }
            scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
        }
        minlmresults(&state, &x, &rep, _state);
        referror = (((referror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-2, _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-2, _state),0.001);
    }
    
    /*
     * 1D problem #1
     *
     * NOTES: we also test negative RKind's corresponding to "inexact" schemes
     * which use approximate finite difference Jacobian.
     */
    for(rkind=-2; rkind<=5; rkind++)
    {
        ae_vector_set_length(&x, 1, _state);
        n = 1;
        m = 1;
        h = 0.00001;
        x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
        if( rkind==-2 )
        {
            minlmcreatev(n, m, &x, h, &state, _state);
            minlmsetacctype(&state, 1, _state);
        }
        if( rkind==-1 )
        {
            minlmcreatev(n, m, &x, h, &state, _state);
            minlmsetacctype(&state, 0, _state);
        }
        if( rkind==0 )
        {
            minlmcreatefj(n, m, &x, &state, _state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, &x, &state, _state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, &x, &state, _state);
        }
        if( rkind==3 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 0, _state);
        }
        if( rkind==4 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 1, _state);
        }
        if( rkind==5 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 2, _state);
        }
        while(minlmiteration(&state, _state))
        {
            if( state.needfi )
            {
                state.fi.ptr.p_double[0] = ae_sin(state.x.ptr.p_double[0], _state);
            }
            if( state.needfij )
            {
                state.fi.ptr.p_double[0] = ae_sin(state.x.ptr.p_double[0], _state);
                state.j.ptr.pp_double[0][0] = ae_cos(state.x.ptr.p_double[0], _state);
            }
            if( (state.needf||state.needfg)||state.needfgh )
            {
                state.f = ae_sqr(ae_sin(state.x.ptr.p_double[0], _state), _state);
            }
            if( state.needfg||state.needfgh )
            {
                state.g.ptr.p_double[0] = 2*ae_sin(state.x.ptr.p_double[0], _state)*ae_cos(state.x.ptr.p_double[0], _state);
            }
            if( state.needfgh )
            {
                state.h.ptr.pp_double[0][0] = 2*(ae_cos(state.x.ptr.p_double[0], _state)*ae_cos(state.x.ptr.p_double[0], _state)-ae_sin(state.x.ptr.p_double[0], _state)*ae_sin(state.x.ptr.p_double[0], _state));
            }
            scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
        }
        minlmresults(&state, &x, &rep, _state);
        lin1error = rep.terminationtype<=0||ae_fp_greater(ae_fabs(x.ptr.p_double[0]/ae_pi-ae_round(x.ptr.p_double[0]/ae_pi, _state), _state),0.001);
    }
    
    /*
     * Linear equations: test normal optimization and optimization with restarts
     */
    for(n=1; n<=10; n++)
    {
        
        /*
         * Prepare task
         */
        h = 0.00001;
        rmatrixrndcond(n, 100, &a, _state);
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xe, n, _state);
        ae_vector_set_length(&b, n, _state);
        for(i=0; i<=n-1; i++)
        {
            xe.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        for(i=0; i<=n-1; i++)
        {
            v = ae_v_dotproduct(&a.ptr.pp_double[i][0], 1, &xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
            b.ptr.p_double[i] = v;
        }
        
        /*
         * Test different RKind
         *
         * NOTES: we also test negative RKind's corresponding to "inexact" schemes
         * which use approximate finite difference Jacobian.
         */
        for(rkind=-2; rkind<=5; rkind++)
        {
            
            /*
             * Solve task (first attempt)
             */
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            if( rkind==-2 )
            {
                minlmcreatev(n, n, &x, h, &state, _state);
                minlmsetacctype(&state, 1, _state);
            }
            if( rkind==-1 )
            {
                minlmcreatev(n, n, &x, h, &state, _state);
                minlmsetacctype(&state, 0, _state);
            }
            if( rkind==0 )
            {
                minlmcreatefj(n, n, &x, &state, _state);
            }
            if( rkind==1 )
            {
                minlmcreatefgj(n, n, &x, &state, _state);
            }
            if( rkind==2 )
            {
                minlmcreatefgh(n, &x, &state, _state);
            }
            if( rkind==3 )
            {
                minlmcreatevj(n, n, &x, &state, _state);
                minlmsetacctype(&state, 0, _state);
            }
            if( rkind==4 )
            {
                minlmcreatevj(n, n, &x, &state, _state);
                minlmsetacctype(&state, 1, _state);
            }
            if( rkind==5 )
            {
                minlmcreatevj(n, n, &x, &state, _state);
                minlmsetacctype(&state, 2, _state);
            }
            while(minlmiteration(&state, _state))
            {
                testminlmunit_axmb(&state, &a, &b, n, _state);
                scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
            }
            minlmresults(&state, &x, &rep, _state);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i=0; i<=n-1; i++)
            {
                eqerror = eqerror||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-xe.ptr.p_double[i], _state),0.001);
            }
            
            /*
             * Now we try to restart algorithm from new point
             */
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            minlmrestartfrom(&state, &x, _state);
            while(minlmiteration(&state, _state))
            {
                testminlmunit_axmb(&state, &a, &b, n, _state);
                scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
            }
            minlmresults(&state, &x, &rep, _state);
            restartserror = restartserror||rep.terminationtype<=0;
            for(i=0; i<=n-1; i++)
            {
                restartserror = restartserror||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-xe.ptr.p_double[i], _state),0.001);
            }
        }
    }
    
    /*
     * Testing convergence properties using
     * different optimizer types and different conditions.
     *
     * Only limited subset of optimizers is tested because some
     * optimizers converge too quickly.
     */
    s = 100;
    for(rkind=0; rkind<=5; rkind++)
    {
        
        /*
         * Skip FGH optimizer - it converges too quickly
         */
        if( rkind==2 )
        {
            continue;
        }
        
        /*
         * Test
         */
        for(ckind=0; ckind<=3; ckind++)
        {
            epsg = 0;
            epsf = 0;
            epsx = 0;
            maxits = 0;
            if( ckind==0 )
            {
                epsf = 0.000001;
            }
            if( ckind==1 )
            {
                epsx = 0.000001;
            }
            if( ckind==2 )
            {
                maxits = 2;
            }
            if( ckind==3 )
            {
                epsg = 0.0001;
            }
            ae_vector_set_length(&x, 3, _state);
            n = 3;
            m = 3;
            for(i=0; i<=2; i++)
            {
                x.ptr.p_double[i] = 6;
            }
            if( rkind==0 )
            {
                minlmcreatefj(n, m, &x, &state, _state);
            }
            if( rkind==1 )
            {
                minlmcreatefgj(n, m, &x, &state, _state);
            }
            ae_assert(rkind!=2, "Assertion failed", _state);
            if( rkind==3 )
            {
                minlmcreatevj(n, m, &x, &state, _state);
                minlmsetacctype(&state, 0, _state);
            }
            if( rkind==4 )
            {
                minlmcreatevj(n, m, &x, &state, _state);
                minlmsetacctype(&state, 1, _state);
            }
            if( rkind==5 )
            {
                minlmcreatevj(n, m, &x, &state, _state);
                minlmsetacctype(&state, 2, _state);
            }
            minlmsetcond(&state, epsg, epsf, epsx, maxits, _state);
            while(minlmiteration(&state, _state))
            {
                if( state.needfi||state.needfij )
                {
                    state.fi.ptr.p_double[0] = s*(ae_exp(state.x.ptr.p_double[0], _state)-2);
                    state.fi.ptr.p_double[1] = ae_sqr(state.x.ptr.p_double[1], _state)+1;
                    state.fi.ptr.p_double[2] = state.x.ptr.p_double[2]-state.x.ptr.p_double[0];
                }
                if( state.needfij )
                {
                    state.j.ptr.pp_double[0][0] = s*ae_exp(state.x.ptr.p_double[0], _state);
                    state.j.ptr.pp_double[0][1] = 0;
                    state.j.ptr.pp_double[0][2] = 0;
                    state.j.ptr.pp_double[1][0] = 0;
                    state.j.ptr.pp_double[1][1] = 2*state.x.ptr.p_double[1];
                    state.j.ptr.pp_double[1][2] = 0;
                    state.j.ptr.pp_double[2][0] = -1;
                    state.j.ptr.pp_double[2][1] = 0;
                    state.j.ptr.pp_double[2][2] = 1;
                }
                if( (state.needf||state.needfg)||state.needfgh )
                {
                    state.f = s*ae_sqr(ae_exp(state.x.ptr.p_double[0], _state)-2, _state)+ae_sqr(ae_sqr(state.x.ptr.p_double[1], _state)+1, _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
                }
                if( state.needfg||state.needfgh )
                {
                    state.g.ptr.p_double[0] = s*2*(ae_exp(state.x.ptr.p_double[0], _state)-2)*ae_exp(state.x.ptr.p_double[0], _state)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
                    state.g.ptr.p_double[1] = 2*(ae_sqr(state.x.ptr.p_double[1], _state)+1)*2*state.x.ptr.p_double[1];
                    state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
                }
                if( state.needfgh )
                {
                    state.h.ptr.pp_double[0][0] = s*(4*ae_sqr(ae_exp(state.x.ptr.p_double[0], _state), _state)-4*ae_exp(state.x.ptr.p_double[0], _state))+2;
                    state.h.ptr.pp_double[0][1] = 0;
                    state.h.ptr.pp_double[0][2] = -2;
                    state.h.ptr.pp_double[1][0] = 0;
                    state.h.ptr.pp_double[1][1] = 12*ae_sqr(state.x.ptr.p_double[1], _state)+4;
                    state.h.ptr.pp_double[1][2] = 0;
                    state.h.ptr.pp_double[2][0] = -2;
                    state.h.ptr.pp_double[2][1] = 0;
                    state.h.ptr.pp_double[2][2] = 2;
                }
                scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
            }
            minlmresults(&state, &x, &rep, _state);
            if( ckind==0 )
            {
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.05);
                converror = converror||rep.terminationtype!=1;
            }
            if( ckind==1 )
            {
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.05);
                converror = converror||rep.terminationtype!=2;
            }
            if( ckind==2 )
            {
                converror = (converror||rep.terminationtype!=5)||rep.iterationscount!=maxits;
            }
            if( ckind==3 )
            {
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.05);
                converror = converror||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.05);
                converror = converror||rep.terminationtype!=4;
            }
        }
    }
    
    /*
     * Other properties:
     * 1. test reports (F should form monotone sequence)
     * 2. test maximum step
     */
    for(rkind=0; rkind<=5; rkind++)
    {
        
        /*
         * reports:
         * * check that first report is initial point
         * * check that F is monotone decreasing
         * * check that last report is final result
         */
        n = 3;
        m = 3;
        s = 100;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xlast, n, _state);
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = 6;
        }
        if( rkind==0 )
        {
            minlmcreatefj(n, m, &x, &state, _state);
        }
        if( rkind==1 )
        {
            minlmcreatefgj(n, m, &x, &state, _state);
        }
        if( rkind==2 )
        {
            minlmcreatefgh(n, &x, &state, _state);
        }
        if( rkind==3 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 0, _state);
        }
        if( rkind==4 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 1, _state);
        }
        if( rkind==5 )
        {
            minlmcreatevj(n, m, &x, &state, _state);
            minlmsetacctype(&state, 2, _state);
        }
        minlmsetcond(&state, 0, 0, 0, 4, _state);
        minlmsetxrep(&state, ae_true, _state);
        fprev = ae_maxrealnumber;
        while(minlmiteration(&state, _state))
        {
            if( state.needfi||state.needfij )
            {
                state.fi.ptr.p_double[0] = ae_sqrt(s, _state)*(ae_exp(state.x.ptr.p_double[0], _state)-2);
                state.fi.ptr.p_double[1] = state.x.ptr.p_double[1];
                state.fi.ptr.p_double[2] = state.x.ptr.p_double[2]-state.x.ptr.p_double[0];
            }
            if( state.needfij )
            {
                state.j.ptr.pp_double[0][0] = ae_sqrt(s, _state)*ae_exp(state.x.ptr.p_double[0], _state);
                state.j.ptr.pp_double[0][1] = 0;
                state.j.ptr.pp_double[0][2] = 0;
                state.j.ptr.pp_double[1][0] = 0;
                state.j.ptr.pp_double[1][1] = 1;
                state.j.ptr.pp_double[1][2] = 0;
                state.j.ptr.pp_double[2][0] = -1;
                state.j.ptr.pp_double[2][1] = 0;
                state.j.ptr.pp_double[2][2] = 1;
            }
            if( (state.needf||state.needfg)||state.needfgh )
            {
                state.f = s*ae_sqr(ae_exp(state.x.ptr.p_double[0], _state)-2, _state)+ae_sqr(state.x.ptr.p_double[1], _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
            }
            if( state.needfg||state.needfgh )
            {
                state.g.ptr.p_double[0] = s*2*(ae_exp(state.x.ptr.p_double[0], _state)-2)*ae_exp(state.x.ptr.p_double[0], _state)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
                state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
                state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
            }
            if( state.needfgh )
            {
                state.h.ptr.pp_double[0][0] = s*(4*ae_sqr(ae_exp(state.x.ptr.p_double[0], _state), _state)-4*ae_exp(state.x.ptr.p_double[0], _state))+2;
                state.h.ptr.pp_double[0][1] = 0;
                state.h.ptr.pp_double[0][2] = -2;
                state.h.ptr.pp_double[1][0] = 0;
                state.h.ptr.pp_double[1][1] = 2;
                state.h.ptr.pp_double[1][2] = 0;
                state.h.ptr.pp_double[2][0] = -2;
                state.h.ptr.pp_double[2][1] = 0;
                state.h.ptr.pp_double[2][2] = 2;
            }
            scerror = scerror||!testminlmunit_rkindvsstatecheck(rkind, &state, _state);
            if( state.xupdated )
            {
                othererrors = othererrors||ae_fp_greater(state.f,fprev);
                if( ae_fp_eq(fprev,ae_maxrealnumber) )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        othererrors = othererrors||ae_fp_neq(state.x.ptr.p_double[i],x.ptr.p_double[i]);
                    }
                }
                fprev = state.f;
                ae_v_move(&xlast.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
            }
        }
        minlmresults(&state, &x, &rep, _state);
        for(i=0; i<=n-1; i++)
        {
            othererrors = othererrors||ae_fp_neq(x.ptr.p_double[i],xlast.ptr.p_double[i]);
        }
    }
    n = 1;
    ae_vector_set_length(&x, n, _state);
    x.ptr.p_double[0] = 100;
    stpmax = 0.05+0.05*ae_randomreal(_state);
    minlmcreatefgh(n, &x, &state, _state);
    minlmsetcond(&state, 1.0E-9, 0, 0, 0, _state);
    minlmsetstpmax(&state, stpmax, _state);
    minlmsetxrep(&state, ae_true, _state);
    xprev = x.ptr.p_double[0];
    while(minlmiteration(&state, _state))
    {
        if( (state.needf||state.needfg)||state.needfgh )
        {
            state.f = ae_exp(state.x.ptr.p_double[0], _state)+ae_exp(-state.x.ptr.p_double[0], _state);
        }
        if( state.needfg||state.needfgh )
        {
            state.g.ptr.p_double[0] = ae_exp(state.x.ptr.p_double[0], _state)-ae_exp(-state.x.ptr.p_double[0], _state);
        }
        if( state.needfgh )
        {
            state.h.ptr.pp_double[0][0] = ae_exp(state.x.ptr.p_double[0], _state)+ae_exp(-state.x.ptr.p_double[0], _state);
        }
        othererrors = othererrors||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]-xprev, _state),(1+ae_sqrt(ae_machineepsilon, _state))*stpmax);
        if( state.xupdated )
        {
            xprev = state.x.ptr.p_double[0];
        }
    }
    
    /*
     * end
     */
    waserrors = ((((((referror||lin1error)||lin2error)||eqerror)||converror)||scerror)||othererrors)||restartserror;
    if( !silent )
    {
        printf("TESTING LEVENBERG-MARQUARDT OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("1-D PROBLEM #1:                           ");
        if( lin1error )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("1-D PROBLEM #2:                           ");
        if( lin2error )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LINEAR EQUATIONS:                         ");
        if( eqerror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("RESTARTS:                                 ");
        if( restartserror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("CONVERGENCE PROPERTIES:                   ");
        if( converror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("STATE FIELDS CONSISTENCY:                 ");
        if( scerror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("OTHER PROPERTIES:                         ");
        if( othererrors )
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
    result = !waserrors;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Asserts that State fields are consistent with RKind.
Returns False otherwise.

RKind is an algorithm selector:
* -2 = V, AccType=1
* -1 = V, AccType=0
*  0 = FJ
*  1 = FGJ
*  2 = FGH
*  3 = VJ, AccType=0
*  4 = VJ, AccType=1
*  5 = VJ, AccType=2

*************************************************************************/
static ae_bool testminlmunit_rkindvsstatecheck(ae_int_t rkind,
     minlmstate* state,
     ae_state *_state)
{
    ae_int_t nset;
    ae_bool result;


    nset = 0;
    if( state->needfi )
    {
        nset = nset+1;
    }
    if( state->needf )
    {
        nset = nset+1;
    }
    if( state->needfg )
    {
        nset = nset+1;
    }
    if( state->needfij )
    {
        nset = nset+1;
    }
    if( state->needfgh )
    {
        nset = nset+1;
    }
    if( state->xupdated )
    {
        nset = nset+1;
    }
    if( nset!=1 )
    {
        result = ae_false;
        return result;
    }
    if( rkind==-2 )
    {
        result = state->needfi||state->xupdated;
        return result;
    }
    if( rkind==-1 )
    {
        result = state->needfi||state->xupdated;
        return result;
    }
    if( rkind==0 )
    {
        result = (state->needf||state->needfij)||state->xupdated;
        return result;
    }
    if( rkind==1 )
    {
        result = ((state->needf||state->needfij)||state->needfg)||state->xupdated;
        return result;
    }
    if( rkind==2 )
    {
        result = ((state->needf||state->needfg)||state->needfgh)||state->xupdated;
        return result;
    }
    if( rkind==3 )
    {
        result = (state->needfi||state->needfij)||state->xupdated;
        return result;
    }
    if( rkind==4 )
    {
        result = (state->needfi||state->needfij)||state->xupdated;
        return result;
    }
    if( rkind==5 )
    {
        result = (state->needfi||state->needfij)||state->xupdated;
        return result;
    }
    result = ae_false;
    return result;
}


/*************************************************************************
Calculates FI/F/G/H for problem min(||Ax-b||)
*************************************************************************/
static void testminlmunit_axmb(minlmstate* state,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double v;


    if( (state->needf||state->needfg)||state->needfgh )
    {
        state->f = 0;
    }
    if( state->needfg||state->needfgh )
    {
        for(i=0; i<=n-1; i++)
        {
            state->g.ptr.p_double[i] = 0;
        }
    }
    if( state->needfgh )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                state->h.ptr.pp_double[i][j] = 0;
            }
        }
    }
    for(i=0; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&a->ptr.pp_double[i][0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
        if( (state->needf||state->needfg)||state->needfgh )
        {
            state->f = state->f+ae_sqr(v-b->ptr.p_double[i], _state);
        }
        if( state->needfg||state->needfgh )
        {
            for(j=0; j<=n-1; j++)
            {
                state->g.ptr.p_double[j] = state->g.ptr.p_double[j]+2*(v-b->ptr.p_double[i])*a->ptr.pp_double[i][j];
            }
        }
        if( state->needfgh )
        {
            for(j=0; j<=n-1; j++)
            {
                for(k=0; k<=n-1; k++)
                {
                    state->h.ptr.pp_double[j][k] = state->h.ptr.pp_double[j][k]+2*a->ptr.pp_double[i][j]*a->ptr.pp_double[i][k];
                }
            }
        }
        if( state->needfi )
        {
            state->fi.ptr.p_double[i] = v-b->ptr.p_double[i];
        }
        if( state->needfij )
        {
            state->fi.ptr.p_double[i] = v-b->ptr.p_double[i];
            ae_v_move(&state->j.ptr.pp_double[i][0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        }
    }
}


/*$ End $*/
