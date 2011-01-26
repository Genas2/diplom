

#include <stdafx.h>
#include <stdio.h>
#include "testminasaunit.h"


/*$ Declarations $*/
static void testminasaunit_testfunc1(minasastate* state, ae_state *_state);
static void testminasaunit_testfunc2(minasastate* state, ae_state *_state);
static void testminasaunit_testfunc3(minasastate* state, ae_state *_state);
static void testminasaunit_checkbounds(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_int_t n,
     ae_bool* err,
     ae_state *_state);
static double testminasaunit_asaboundval(double x,
     double b1,
     double b2,
     ae_state *_state);


/*$ Body $*/


ae_bool testminasa(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool referror;
    ae_bool restartserror;
    ae_bool converror;
    ae_bool othererrors;
    ae_int_t n;
    ae_vector x;
    ae_vector xe;
    ae_vector c;
    ae_vector bndl;
    ae_vector bndu;
    ae_vector xlast;
    double fprev;
    double xprev;
    double stpmax;
    ae_int_t i;
    ae_int_t j;
    double v;
    double s;
    double tol;
    ae_int_t algotype;
    ae_matrix a;
    minasastate state;
    minasareport rep;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xe, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bndl, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bndu, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xlast, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    _minasastate_init(&state, _state, ae_true);
    _minasareport_init(&rep, _state, ae_true);

    waserrors = ae_false;
    referror = ae_false;
    converror = ae_false;
    othererrors = ae_false;
    restartserror = ae_false;
    
    /*
     * Different algorithms
     */
    for(algotype=-1; algotype<=1; algotype++)
    {
        
        /*
         * reference problem, simple convex optimization
         */
        for(n=1; n<=5; n++)
        {
            
            /*
             * min(x'*diag(c)*x) on a random box
             */
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&xe, n, _state);
            ae_vector_set_length(&c, n, _state);
            ae_vector_set_length(&bndl, n, _state);
            ae_vector_set_length(&bndu, n, _state);
            for(i=0; i<=n-1; i++)
            {
                c.ptr.p_double[i] = 1+ae_randomreal(_state);
                xe.ptr.p_double[i] = 4*ae_randomreal(_state)-2;
                bndl.ptr.p_double[i] = -ae_maxreal(ae_randomreal(_state), 0.2, _state);
                bndu.ptr.p_double[i] = ae_maxreal(ae_randomreal(_state), 0.2, _state);
                x.ptr.p_double[i] = 0.5*(bndl.ptr.p_double[i]+bndu.ptr.p_double[i]);
            }
            tol = 0.001;
            minasacreate(n, &x, &bndl, &bndu, &state, _state);
            minasasetcond(&state, tol, 0.0, 0.0, 0, _state);
            minasasetalgorithm(&state, algotype, _state);
            while(minasaiteration(&state, _state))
            {
                testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
                state.f = 0;
                for(i=0; i<=n-1; i++)
                {
                    state.f = state.f+c.ptr.p_double[i]*ae_sqr(state.x.ptr.p_double[i]-xe.ptr.p_double[i], _state);
                    state.g.ptr.p_double[i] = 2*c.ptr.p_double[i]*(state.x.ptr.p_double[i]-xe.ptr.p_double[i]);
                }
            }
            minasaresults(&state, &x, &rep, _state);
            referror = referror||rep.terminationtype<=0;
            for(i=0; i<=n-1; i++)
            {
                referror = referror||ae_fp_greater(ae_fabs(testminasaunit_asaboundval(xe.ptr.p_double[i], bndl.ptr.p_double[i], bndu.ptr.p_double[i], _state)-x.ptr.p_double[i], _state),0.01);
            }
        }
        
        /*
         * F2 problem with restarts:
         * * make several iterations and restart BEFORE termination
         * * iterate and restart AFTER termination
         *
         * NOTE: step is bounded from above to avoid premature convergence
         */
        ae_vector_set_length(&x, 3, _state);
        ae_vector_set_length(&bndl, 3, _state);
        ae_vector_set_length(&bndu, 3, _state);
        n = 3;
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        bndl.ptr.p_double[0] = -10000;
        bndl.ptr.p_double[1] = -10000;
        bndl.ptr.p_double[2] = -10000;
        bndu.ptr.p_double[0] = 10000;
        bndu.ptr.p_double[1] = 10000;
        bndu.ptr.p_double[2] = 10000;
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetalgorithm(&state, algotype, _state);
        minasasetstpmax(&state, 0.1, _state);
        minasasetcond(&state, 0.0000001, 0.0, 0.0, 0, _state);
        for(i=0; i<=10; i++)
        {
            if( !minasaiteration(&state, _state) )
            {
                break;
            }
            testminasaunit_testfunc2(&state, _state);
        }
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        minasarestartfrom(&state, &x, &bndl, &bndu, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_testfunc2(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        minasarestartfrom(&state, &x, &bndl, &bndu, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_testfunc2(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
        
        /*
         * reference problem 2: non-convex optimization on [-2,2] x [1,2]
         *
         * A saddle function is minimized:
         * * stationary point [0,0] (non-feasible)
         * * constrained minimum [-2,2].
         * * starting point [+2,2]
         *
         * Path from start to end may be very complex, with multiple changes
         * in active constraints, so it is interesting task for our method.
         *
         * Scale parameter is used to make optimization more interesting
         * during GPA runs.
         */
        ae_vector_set_length(&x, 2, _state);
        ae_vector_set_length(&bndl, 2, _state);
        ae_vector_set_length(&bndu, 2, _state);
        bndl.ptr.p_double[0] = -2;
        bndu.ptr.p_double[0] = 2;
        x.ptr.p_double[0] = 2;
        bndl.ptr.p_double[1] = 1;
        bndu.ptr.p_double[1] = 2;
        x.ptr.p_double[1] = 2;
        tol = 0.001;
        s = 0.01;
        minasacreate(2, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, tol, 0.0, 0.0, 0, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, 2, &othererrors, _state);
            state.f = s*(ae_sqr(state.x.ptr.p_double[0]+state.x.ptr.p_double[1], _state)-ae_sqr(state.x.ptr.p_double[0]-state.x.ptr.p_double[1], _state));
            state.g.ptr.p_double[0] = s*(2*(state.x.ptr.p_double[0]+state.x.ptr.p_double[1])-2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[1]));
            state.g.ptr.p_double[1] = s*(2*(state.x.ptr.p_double[0]+state.x.ptr.p_double[1])+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[1]));
        }
        minasaresults(&state, &x, &rep, _state);
        referror = ((referror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]+2, _state),0.01))||ae_fp_greater(ae_fabs(state.x.ptr.p_double[1]-2, _state),0.01);
        
        /*
         * function #1 with 'x[0]>=ln(2)' constraint.
         * may show very interesting behavior.
         */
        ae_vector_set_length(&x, 3, _state);
        ae_vector_set_length(&bndl, 3, _state);
        ae_vector_set_length(&bndu, 3, _state);
        n = 3;
        for(i=0; i<=2; i++)
        {
            bndl.ptr.p_double[i] = -10000;
            bndu.ptr.p_double[i] = 10000;
        }
        bndl.ptr.p_double[0] = ae_log(2, _state);
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 3*ae_randomreal(_state)+3;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0.0000001, 0.0, 0.0, 0, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
            testminasaunit_testfunc1(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        referror = referror||rep.terminationtype<=0;
        referror = referror||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.05);
        referror = referror||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.05);
        referror = referror||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.05);
        
        /*
         * Testing convergence properties
         */
        ae_vector_set_length(&x, 3, _state);
        ae_vector_set_length(&bndl, 3, _state);
        ae_vector_set_length(&bndu, 3, _state);
        n = 3;
        for(i=0; i<=2; i++)
        {
            bndl.ptr.p_double[i] = -10000;
            bndu.ptr.p_double[i] = 10000;
        }
        bndl.ptr.p_double[0] = ae_log(2, _state);
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 3*ae_randomreal(_state)+3;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0.001, 0.0, 0.0, 0, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
            testminasaunit_testfunc3(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=4;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 3*ae_randomreal(_state)+3;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0.0, 0.001, 0.0, 0, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
            testminasaunit_testfunc3(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=1;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 3*ae_randomreal(_state)+3;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0.0, 0.0, 0.001, 0, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
            testminasaunit_testfunc3(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=2;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 3*ae_randomreal(_state)+3;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0.0, 0.0, 0.0, 3, _state);
        minasasetalgorithm(&state, algotype, _state);
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
            testminasaunit_testfunc3(&state, _state);
        }
        minasaresults(&state, &x, &rep, _state);
        converror = converror||!((rep.terminationtype==5&&rep.iterationscount==3)||rep.terminationtype==7);
        
        /*
         * Other properties
         *
         *
         * Other properties:
         * 1. test reports (F should form monotone sequence)
         * 2. test maximum step
         */
        n = 50;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xlast, n, _state);
        ae_vector_set_length(&bndl, n, _state);
        ae_vector_set_length(&bndu, n, _state);
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = 1;
            xlast.ptr.p_double[i] = ae_randomreal(_state);
            bndl.ptr.p_double[i] = -100000;
            bndu.ptr.p_double[i] = 100000;
        }
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 0, 0, 0, 100, _state);
        minasasetxrep(&state, ae_true, _state);
        fprev = ae_maxrealnumber;
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
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
        minasaresults(&state, &x, &rep, _state);
        for(i=0; i<=n-1; i++)
        {
            othererrors = othererrors||ae_fp_neq(x.ptr.p_double[i],xlast.ptr.p_double[i]);
        }
        n = 1;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&bndl, n, _state);
        ae_vector_set_length(&bndu, n, _state);
        x.ptr.p_double[0] = 100;
        bndl.ptr.p_double[0] = -1000000;
        bndu.ptr.p_double[0] = 1000000;
        stpmax = 0.05+0.05*ae_randomreal(_state);
        minasacreate(n, &x, &bndl, &bndu, &state, _state);
        minasasetcond(&state, 1.0E-9, 0, 0, 0, _state);
        minasasetstpmax(&state, stpmax, _state);
        minasasetxrep(&state, ae_true, _state);
        xprev = x.ptr.p_double[0];
        while(minasaiteration(&state, _state))
        {
            testminasaunit_checkbounds(&state.x, &bndl, &bndu, n, &othererrors, _state);
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
    }
    
    /*
     * end
     */
    waserrors = ((referror||converror)||othererrors)||restartserror;
    if( !silent )
    {
        printf("TESTING ASA OPTIMIZATION\n");
        printf("REFERENCE PROBLEMS:                       ");
        if( referror )
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
static void testminasaunit_testfunc1(minasastate* state, ae_state *_state)
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
static void testminasaunit_testfunc2(minasastate* state, ae_state *_state)
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
static void testminasaunit_testfunc3(minasastate* state, ae_state *_state)
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
Checks that X is bounded with respect to BndL/BndU.

If it is not, True is assigned to the Err variable (which is not changed
otherwise).
*************************************************************************/
static void testminasaunit_checkbounds(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_int_t n,
     ae_bool* err,
     ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_less(x->ptr.p_double[i],bndl->ptr.p_double[i])||ae_fp_greater(x->ptr.p_double[i],bndu->ptr.p_double[i]) )
        {
            *err = ae_true;
        }
    }
}


/*************************************************************************
'bound' value: map X to [B1,B2]
*************************************************************************/
static double testminasaunit_asaboundval(double x,
     double b1,
     double b2,
     ae_state *_state)
{
    double result;


    if( ae_fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ae_fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


/*$ End $*/
