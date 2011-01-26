

#include <stdafx.h>
#include <stdio.h>
#include "testminlbfgsunit.h"


/*$ Declarations $*/
static void testminlbfgsunit_testfunc1(minlbfgsstate* state,
     ae_state *_state);
static void testminlbfgsunit_testfunc2(minlbfgsstate* state,
     ae_state *_state);
static void testminlbfgsunit_testfunc3(minlbfgsstate* state,
     ae_state *_state);
static void testminlbfgsunit_calciip2(minlbfgsstate* state,
     ae_int_t n,
     ae_state *_state);


/*$ Body $*/


ae_bool testminlbfgs(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool referror;
    ae_bool nonconverror;
    ae_bool eqerror;
    ae_bool converror;
    ae_bool crashtest;
    ae_bool othererrors;
    ae_bool restartserror;
    ae_bool precerror;
    ae_int_t n;
    ae_int_t m;
    ae_vector x;
    ae_vector xe;
    ae_vector b;
    ae_vector xlast;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t pass;
    double v;
    ae_matrix a;
    ae_matrix aa;
    ae_int_t cnt1;
    ae_int_t cnt2;
    ae_int_t cnt3;
    double epsg;
    ae_int_t maxits;
    minlbfgsstate state;
    minlbfgsreport rep;
    double fprev;
    double xprev;
    double stpmax;
    ae_int_t pkind;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xe, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xlast, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&aa, 0, 0, DT_REAL, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);
    _minlbfgsreport_init(&rep, _state, ae_true);

    waserrors = ae_false;
    precerror = ae_false;
    nonconverror = ae_false;
    restartserror = ae_false;
    eqerror = ae_false;
    converror = ae_false;
    crashtest = ae_false;
    othererrors = ae_false;
    
    /*
     * Reference problem
     */
    ae_vector_set_length(&x, 2+1, _state);
    n = 3;
    m = 2;
    x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
    x.ptr.p_double[1] = 100*ae_randomreal(_state)-50;
    x.ptr.p_double[2] = 100*ae_randomreal(_state)-50;
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0, 0, 0, _state);
    while(minlbfgsiteration(&state, _state))
    {
        state.f = ae_sqr(state.x.ptr.p_double[0]-2, _state)+ae_sqr(state.x.ptr.p_double[1], _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
        state.g.ptr.p_double[0] = 2*(state.x.ptr.p_double[0]-2)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
        state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
        state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    referror = ((rep.terminationtype<=0||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-2, _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-2, _state),0.001);
    
    /*
     * Preconditioner test.
     *
     * If
     * * P1 is default preconditioner
     * * P2 is Cholesky preconditioner with unit diagonal
     * * P3 is Cholesky preconditioner based on exact Hessian with perturbations
     * then P1 is worse than P3, P2 is worse than P3.
     *
     * We test it using f(x) = sum( ((i*i+1)*x[i])^2, i=0..N-1) and L-BFGS
     * optimizer with deliberately small M=1.
     *
     * N        - problem size
     * PKind    - zero for upper triangular preconditioner, one for lower triangular.
     * K        - number of repeated passes (should be large enough to average out random factors)
     */
    m = 1;
    k = 20;
    epsg = 1.0E-10;
    for(n=5; n<=15; n++)
    {
        for(pkind=0; pkind<=1; pkind++)
        {
            ae_matrix_set_length(&a, n, n, _state);
            ae_matrix_set_length(&aa, n, n, _state);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( i==j )
                    {
                        aa.ptr.pp_double[i][i] = 1;
                    }
                    else
                    {
                        aa.ptr.pp_double[i][j] = 0;
                    }
                    if( i==j )
                    {
                        a.ptr.pp_double[i][i] = (i*i+1)*(0.8+0.4*ae_randomreal(_state));
                    }
                    else
                    {
                        if( (pkind==0&&j>i)||(pkind==1&&j<i) )
                        {
                            a.ptr.pp_double[i][j] = 0.1*ae_randomreal(_state)-0.05;
                        }
                        else
                        {
                            a.ptr.pp_double[i][j] = _state->v_nan;
                        }
                    }
                }
            }
            ae_vector_set_length(&x, n, _state);
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 0;
            }
            minlbfgscreate(n, m, &x, &state, _state);
            
            /*
             * Test it with default preconditioner
             */
            minlbfgssetdefaultpreconditioner(&state, _state);
            cnt1 = 0;
            for(pass=0; pass<=k-1; pass++)
            {
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                }
                minlbfgsrestartfrom(&state, &x, _state);
                while(minlbfgsiteration(&state, _state))
                {
                    testminlbfgsunit_calciip2(&state, n, _state);
                }
                minlbfgsresults(&state, &x, &rep, _state);
                cnt1 = cnt1+rep.iterationscount;
            }
            
            /*
             * Test it with unit preconditioner
             */
            minlbfgssetcholeskypreconditioner(&state, &aa, pkind==0, _state);
            cnt2 = 0;
            for(pass=0; pass<=k-1; pass++)
            {
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                }
                minlbfgsrestartfrom(&state, &x, _state);
                while(minlbfgsiteration(&state, _state))
                {
                    testminlbfgsunit_calciip2(&state, n, _state);
                }
                minlbfgsresults(&state, &x, &rep, _state);
                cnt2 = cnt2+rep.iterationscount;
            }
            
            /*
             * Test it with perturbed preconditioner
             */
            minlbfgssetcholeskypreconditioner(&state, &a, pkind==0, _state);
            cnt3 = 0;
            for(pass=0; pass<=k-1; pass++)
            {
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                }
                minlbfgsrestartfrom(&state, &x, _state);
                while(minlbfgsiteration(&state, _state))
                {
                    testminlbfgsunit_calciip2(&state, n, _state);
                }
                minlbfgsresults(&state, &x, &rep, _state);
                cnt3 = cnt3+rep.iterationscount;
            }
            
            /*
             * Compare
             */
            precerror = precerror||cnt1<cnt3;
            precerror = precerror||cnt2<cnt3;
        }
    }
    
    /*
     * nonconvex problems with complex surface: we start from point with very small
     * gradient, but we need ever smaller gradient in the next step due to
     * Wolfe conditions.
     */
    ae_vector_set_length(&x, 1, _state);
    n = 1;
    m = 1;
    v = -100;
    while(ae_fp_less(v,0.1))
    {
        x.ptr.p_double[0] = v;
        minlbfgscreate(n, m, &x, &state, _state);
        minlbfgssetcond(&state, 1.0E-9, 0, 0, 0, _state);
        while(minlbfgsiteration(&state, _state))
        {
            state.f = ae_sqr(state.x.ptr.p_double[0], _state)/(1+ae_sqr(state.x.ptr.p_double[0], _state));
            state.g.ptr.p_double[0] = (2*state.x.ptr.p_double[0]*(1+ae_sqr(state.x.ptr.p_double[0], _state))-ae_sqr(state.x.ptr.p_double[0], _state)*2*state.x.ptr.p_double[0])/ae_sqr(1+ae_sqr(state.x.ptr.p_double[0], _state), _state);
        }
        minlbfgsresults(&state, &x, &rep, _state);
        nonconverror = (nonconverror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0], _state),0.001);
        v = v+0.1;
    }
    
    /*
     * F2 problem with restarts:
     * * make several iterations and restart BEFORE termination
     * * iterate and restart AFTER termination
     *
     * NOTE: step is bounded from above to avoid premature convergence
     */
    ae_vector_set_length(&x, 3, _state);
    n = 3;
    m = 2;
    x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetstpmax(&state, 0.1, _state);
    minlbfgssetcond(&state, 0.0000001, 0.0, 0.0, 0, _state);
    for(i=0; i<=10; i++)
    {
        if( !minlbfgsiteration(&state, _state) )
        {
            break;
        }
        testminlbfgsunit_testfunc2(&state, _state);
    }
    x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
    minlbfgsrestartfrom(&state, &x, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc2(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
    x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
    x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
    minlbfgsrestartfrom(&state, &x, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc2(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
    
    /*
     * Linear equations
     */
    for(n=1; n<=10; n++)
    {
        
        /*
         * Prepare task
         */
        ae_matrix_set_length(&a, n-1+1, n-1+1, _state);
        ae_vector_set_length(&x, n-1+1, _state);
        ae_vector_set_length(&xe, n-1+1, _state);
        ae_vector_set_length(&b, n-1+1, _state);
        for(i=0; i<=n-1; i++)
        {
            xe.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
            }
            a.ptr.pp_double[i][i] = a.ptr.pp_double[i][i]+3*ae_sign(a.ptr.pp_double[i][i], _state);
        }
        for(i=0; i<=n-1; i++)
        {
            v = ae_v_dotproduct(&a.ptr.pp_double[i][0], 1, &xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
            b.ptr.p_double[i] = v;
        }
        
        /*
         * Test different M
         */
        for(m=1; m<=n; m++)
        {
            
            /*
             * Solve task
             */
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            minlbfgscreate(n, m, &x, &state, _state);
            minlbfgssetcond(&state, 0, 0, 0, 0, _state);
            while(minlbfgsiteration(&state, _state))
            {
                state.f = 0;
                for(i=0; i<=n-1; i++)
                {
                    state.g.ptr.p_double[i] = 0;
                }
                for(i=0; i<=n-1; i++)
                {
                    v = ae_v_dotproduct(&a.ptr.pp_double[i][0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    state.f = state.f+ae_sqr(v-b.ptr.p_double[i], _state);
                    for(j=0; j<=n-1; j++)
                    {
                        state.g.ptr.p_double[j] = state.g.ptr.p_double[j]+2*(v-b.ptr.p_double[i])*a.ptr.pp_double[i][j];
                    }
                }
            }
            minlbfgsresults(&state, &x, &rep, _state);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i=0; i<=n-1; i++)
            {
                eqerror = eqerror||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-xe.ptr.p_double[i], _state),0.001);
            }
        }
    }
    
    /*
     * Testing convergence properties
     */
    ae_vector_set_length(&x, 2+1, _state);
    n = 3;
    m = 2;
    for(i=0; i<=2; i++)
    {
        x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0.001, 0, 0, 0, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc3(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    converror = converror||rep.terminationtype!=4;
    for(i=0; i<=2; i++)
    {
        x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0.001, 0, 0, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc3(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    converror = converror||rep.terminationtype!=1;
    for(i=0; i<=2; i++)
    {
        x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0, 0.001, 0, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc3(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    converror = converror||rep.terminationtype!=2;
    for(i=0; i<=2; i++)
    {
        x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0, 0, 10, _state);
    while(minlbfgsiteration(&state, _state))
    {
        testminlbfgsunit_testfunc3(&state, _state);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    converror = (converror||rep.terminationtype!=5)||rep.iterationscount!=10;
    
    /*
     * Crash test: too many iterations on a simple tasks
     * May fail when encounter zero step, underflow or something like that
     */
    ae_vector_set_length(&x, 2+1, _state);
    n = 3;
    m = 2;
    maxits = 10000;
    for(i=0; i<=2; i++)
    {
        x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0, 0, maxits, _state);
    while(minlbfgsiteration(&state, _state))
    {
        state.f = ae_sqr(ae_exp(state.x.ptr.p_double[0], _state)-2, _state)+ae_sqr(state.x.ptr.p_double[1], _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
        state.g.ptr.p_double[0] = 2*(ae_exp(state.x.ptr.p_double[0], _state)-2)*ae_exp(state.x.ptr.p_double[0], _state)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
        state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
        state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
    }
    minlbfgsresults(&state, &x, &rep, _state);
    crashtest = crashtest||rep.terminationtype<=0;
    
    /*
     * Other properties:
     * 1. test reports (F should form monotone sequence)
     * 2. test maximum step
     */
    n = 50;
    m = 2;
    ae_vector_set_length(&x, n, _state);
    ae_vector_set_length(&xlast, n, _state);
    for(i=0; i<=n-1; i++)
    {
        x.ptr.p_double[i] = 1;
    }
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 0, 0, 0, 100, _state);
    minlbfgssetxrep(&state, ae_true, _state);
    fprev = ae_maxrealnumber;
    while(minlbfgsiteration(&state, _state))
    {
        if( state.needfg )
        {
            state.f = 0;
            for(i=0; i<=n-1; i++)
            {
                state.f = state.f+ae_sqr((1+i)*state.x.ptr.p_double[i], _state);
                state.g.ptr.p_double[i] = 2*(1+i)*state.x.ptr.p_double[i];
            }
        }
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
    minlbfgsresults(&state, &x, &rep, _state);
    for(i=0; i<=n-1; i++)
    {
        othererrors = othererrors||ae_fp_neq(x.ptr.p_double[i],xlast.ptr.p_double[i]);
    }
    n = 1;
    m = 1;
    ae_vector_set_length(&x, n, _state);
    x.ptr.p_double[0] = 100;
    stpmax = 0.05+0.05*ae_randomreal(_state);
    minlbfgscreate(n, m, &x, &state, _state);
    minlbfgssetcond(&state, 1.0E-9, 0, 0, 0, _state);
    minlbfgssetstpmax(&state, stpmax, _state);
    minlbfgssetxrep(&state, ae_true, _state);
    xprev = x.ptr.p_double[0];
    while(minlbfgsiteration(&state, _state))
    {
        if( state.needfg )
        {
            state.f = ae_exp(state.x.ptr.p_double[0], _state)+ae_exp(-state.x.ptr.p_double[0], _state);
            state.g.ptr.p_double[0] = ae_exp(state.x.ptr.p_double[0], _state)-ae_exp(-state.x.ptr.p_double[0], _state);
            othererrors = othererrors||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]-xprev, _state),(1+ae_sqrt(ae_machineepsilon, _state))*stpmax);
        }
        if( state.xupdated )
        {
            othererrors = othererrors||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]-xprev, _state),(1+ae_sqrt(ae_machineepsilon, _state))*stpmax);
            xprev = state.x.ptr.p_double[0];
        }
    }
    
    /*
     * end
     */
    waserrors = ((((((referror||nonconverror)||eqerror)||converror)||crashtest)||othererrors)||restartserror)||precerror;
    if( !silent )
    {
        printf("TESTING L-BFGS OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("NON-CONVEX PROBLEM:                       ");
        if( nonconverror )
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
        printf("PRECONDITIONER:                           ");
        if( precerror )
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
        printf("CRASH TEST:                               ");
        if( crashtest )
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
Calculate test function #1

It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
constraint.
*************************************************************************/
static void testminlbfgsunit_testfunc1(minlbfgsstate* state,
     ae_state *_state)
{


    if( ae_fp_less(state->x.ptr.p_double[0],100) )
    {
        state->f = ae_sqr(ae_exp(state->x.ptr.p_double[0], _state)-2, _state)+ae_sqr(state->x.ptr.p_double[1], _state)+ae_sqr(state->x.ptr.p_double[2]-state->x.ptr.p_double[0], _state);
        state->g.ptr.p_double[0] = 2*(ae_exp(state->x.ptr.p_double[0], _state)-2)*ae_exp(state->x.ptr.p_double[0], _state)+2*(state->x.ptr.p_double[0]-state->x.ptr.p_double[2]);
        state->g.ptr.p_double[1] = 2*state->x.ptr.p_double[1];
        state->g.ptr.p_double[2] = 2*(state->x.ptr.p_double[2]-state->x.ptr.p_double[0]);
    }
    else
    {
        state->f = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[0] = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[1] = 0;
        state->g.ptr.p_double[2] = 0;
    }
}


/*************************************************************************
Calculate test function #2

Simple variation of #1, much more nonlinear, which makes unlikely premature
convergence of algorithm .
*************************************************************************/
static void testminlbfgsunit_testfunc2(minlbfgsstate* state,
     ae_state *_state)
{


    if( ae_fp_less(state->x.ptr.p_double[0],100) )
    {
        state->f = ae_sqr(ae_exp(state->x.ptr.p_double[0], _state)-2, _state)+ae_sqr(ae_sqr(state->x.ptr.p_double[1], _state), _state)+ae_sqr(state->x.ptr.p_double[2]-state->x.ptr.p_double[0], _state);
        state->g.ptr.p_double[0] = 2*(ae_exp(state->x.ptr.p_double[0], _state)-2)*ae_exp(state->x.ptr.p_double[0], _state)+2*(state->x.ptr.p_double[0]-state->x.ptr.p_double[2]);
        state->g.ptr.p_double[1] = 4*state->x.ptr.p_double[1]*ae_sqr(state->x.ptr.p_double[1], _state);
        state->g.ptr.p_double[2] = 2*(state->x.ptr.p_double[2]-state->x.ptr.p_double[0]);
    }
    else
    {
        state->f = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[0] = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[1] = 0;
        state->g.ptr.p_double[2] = 0;
    }
}


/*************************************************************************
Calculate test function #3

Simple variation of #1, much more nonlinear, with non-zero value at minimum.
It achieve two goals:
* makes unlikely premature convergence of algorithm .
* solves some issues with EpsF stopping condition which arise when
  F(minimum) is zero

*************************************************************************/
static void testminlbfgsunit_testfunc3(minlbfgsstate* state,
     ae_state *_state)
{
    double s;


    s = 0.001;
    if( ae_fp_less(state->x.ptr.p_double[0],100) )
    {
        state->f = ae_sqr(ae_exp(state->x.ptr.p_double[0], _state)-2, _state)+ae_sqr(ae_sqr(state->x.ptr.p_double[1], _state)+s, _state)+ae_sqr(state->x.ptr.p_double[2]-state->x.ptr.p_double[0], _state);
        state->g.ptr.p_double[0] = 2*(ae_exp(state->x.ptr.p_double[0], _state)-2)*ae_exp(state->x.ptr.p_double[0], _state)+2*(state->x.ptr.p_double[0]-state->x.ptr.p_double[2]);
        state->g.ptr.p_double[1] = 2*(ae_sqr(state->x.ptr.p_double[1], _state)+s)*2*state->x.ptr.p_double[1];
        state->g.ptr.p_double[2] = 2*(state->x.ptr.p_double[2]-state->x.ptr.p_double[0]);
    }
    else
    {
        state->f = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[0] = ae_sqrt(ae_maxrealnumber, _state);
        state->g.ptr.p_double[1] = 0;
        state->g.ptr.p_double[2] = 0;
    }
}


/*************************************************************************
Calculate test function IIP2

f(x) = sum( ((i*i+1)*x[i])^2, i=0..N-1)

It has high condition number which makes unlikely fast convergence without
good preconditioner.

*************************************************************************/
static void testminlbfgsunit_calciip2(minlbfgsstate* state,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( state->needfg )
    {
        state->f = 0;
    }
    for(i=0; i<=n-1; i++)
    {
        if( state->needfg )
        {
            state->f = state->f+ae_sqr(i*i+1, _state)*ae_sqr(state->x.ptr.p_double[i], _state);
            state->g.ptr.p_double[i] = ae_sqr(i*i+1, _state)*2*state->x.ptr.p_double[i];
        }
    }
}


/*$ End $*/
