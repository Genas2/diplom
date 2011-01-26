

#include <stdafx.h>
#include <stdio.h>
#include "testmincgunit.h"


/*$ Declarations $*/
static void testmincgunit_testfunc1(mincgstate* state, ae_state *_state);
static void testmincgunit_testfunc2(mincgstate* state, ae_state *_state);
static void testmincgunit_testfunc3(mincgstate* state, ae_state *_state);


/*$ Body $*/


ae_bool testmincg(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool referror;
    ae_bool eqerror;
    ae_bool linerror1;
    ae_bool linerror2;
    ae_bool restartserror;
    ae_bool converror;
    ae_bool othererrors;
    ae_int_t n;
    ae_vector x;
    ae_vector xe;
    ae_vector b;
    ae_vector xlast;
    double fprev;
    double xprev;
    double stpmax;
    ae_int_t i;
    ae_int_t j;
    double v;
    ae_matrix a;
    mincgstate state;
    mincgreport rep;
    ae_int_t cgtype;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xe, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xlast, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    _mincgstate_init(&state, _state, ae_true);
    _mincgreport_init(&rep, _state, ae_true);

    waserrors = ae_false;
    referror = ae_false;
    linerror1 = ae_false;
    linerror2 = ae_false;
    eqerror = ae_false;
    converror = ae_false;
    restartserror = ae_false;
    othererrors = ae_false;
    for(cgtype=0; cgtype<=1; cgtype++)
    {
        
        /*
         * Reference problem
         */
        ae_vector_set_length(&x, 2+1, _state);
        n = 3;
        x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
        x.ptr.p_double[1] = 100*ae_randomreal(_state)-50;
        x.ptr.p_double[2] = 100*ae_randomreal(_state)-50;
        mincgcreate(n, &x, &state, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            state.f = ae_sqr(state.x.ptr.p_double[0]-2, _state)+ae_sqr(state.x.ptr.p_double[1], _state)+ae_sqr(state.x.ptr.p_double[2]-state.x.ptr.p_double[0], _state);
            state.g.ptr.p_double[0] = 2*(state.x.ptr.p_double[0]-2)+2*(state.x.ptr.p_double[0]-state.x.ptr.p_double[2]);
            state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
            state.g.ptr.p_double[2] = 2*(state.x.ptr.p_double[2]-state.x.ptr.p_double[0]);
        }
        mincgresults(&state, &x, &rep, _state);
        referror = (((referror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-2, _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.001))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-2, _state),0.001);
        
        /*
         * F2 problem with restarts:
         * * make several iterations and restart BEFORE termination
         * * iterate and restart AFTER termination
         *
         * NOTE: step is bounded from above to avoid premature convergence
         */
        ae_vector_set_length(&x, 3, _state);
        n = 3;
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        mincgcreate(n, &x, &state, _state);
        mincgsetcgtype(&state, cgtype, _state);
        mincgsetstpmax(&state, 0.1, _state);
        mincgsetcond(&state, 0.0000001, 0.0, 0.0, 0, _state);
        for(i=0; i<=10; i++)
        {
            if( !mincgiteration(&state, _state) )
            {
                break;
            }
            testmincgunit_testfunc2(&state, _state);
        }
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        mincgrestartfrom(&state, &x, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc2(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
        x.ptr.p_double[0] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[1] = 10+10*ae_randomreal(_state);
        x.ptr.p_double[2] = 10+10*ae_randomreal(_state);
        mincgrestartfrom(&state, &x, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc2(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        restartserror = (((restartserror||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-ae_log(2, _state), _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),0.01))||ae_fp_greater(ae_fabs(x.ptr.p_double[2]-ae_log(2, _state), _state),0.01);
        
        /*
         * 1D problem #1
         */
        ae_vector_set_length(&x, 0+1, _state);
        n = 1;
        x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
        mincgcreate(n, &x, &state, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            state.f = -ae_cos(state.x.ptr.p_double[0], _state);
            state.g.ptr.p_double[0] = ae_sin(state.x.ptr.p_double[0], _state);
        }
        mincgresults(&state, &x, &rep, _state);
        linerror1 = (linerror1||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0]/ae_pi-ae_round(x.ptr.p_double[0]/ae_pi, _state), _state),0.001);
        
        /*
         * 1D problem #2
         */
        ae_vector_set_length(&x, 0+1, _state);
        n = 1;
        x.ptr.p_double[0] = 100*ae_randomreal(_state)-50;
        mincgcreate(n, &x, &state, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            state.f = ae_sqr(state.x.ptr.p_double[0], _state)/(1+ae_sqr(state.x.ptr.p_double[0], _state));
            state.g.ptr.p_double[0] = (2*state.x.ptr.p_double[0]*(1+ae_sqr(state.x.ptr.p_double[0], _state))-ae_sqr(state.x.ptr.p_double[0], _state)*2*state.x.ptr.p_double[0])/ae_sqr(1+ae_sqr(state.x.ptr.p_double[0], _state), _state);
        }
        mincgresults(&state, &x, &rep, _state);
        linerror2 = (linerror2||rep.terminationtype<=0)||ae_fp_greater(ae_fabs(x.ptr.p_double[0], _state),0.001);
        
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
             * Solve task
             */
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            mincgcreate(n, &x, &state, _state);
            mincgsetcgtype(&state, cgtype, _state);
            while(mincgiteration(&state, _state))
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
            mincgresults(&state, &x, &rep, _state);
            eqerror = eqerror||rep.terminationtype<=0;
            for(i=0; i<=n-1; i++)
            {
                eqerror = eqerror||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-xe.ptr.p_double[i], _state),0.001);
            }
        }
        
        /*
         * Testing convergence properties
         */
        ae_vector_set_length(&x, 2+1, _state);
        n = 3;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
        }
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 0.001, 0.0, 0.0, 0, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc3(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=4;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
        }
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 0.0, 0.001, 0.0, 0, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc3(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=1;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 6*ae_randomreal(_state)-3;
        }
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 0.0, 0.0, 0.001, 0, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc3(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        converror = converror||rep.terminationtype!=2;
        for(i=0; i<=2; i++)
        {
            x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 0.0, 0.0, 0.0, 10, _state);
        mincgsetcgtype(&state, cgtype, _state);
        while(mincgiteration(&state, _state))
        {
            testmincgunit_testfunc3(&state, _state);
        }
        mincgresults(&state, &x, &rep, _state);
        converror = converror||!((rep.terminationtype==5&&rep.iterationscount==10)||rep.terminationtype==7);
        
        /*
         * Other properties:
         * 1. test reports (F should form monotone sequence)
         * 2. test maximum step
         */
        n = 50;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xlast, n, _state);
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = 1;
        }
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 0, 0, 0, 100, _state);
        mincgsetxrep(&state, ae_true, _state);
        fprev = ae_maxrealnumber;
        while(mincgiteration(&state, _state))
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
        mincgresults(&state, &x, &rep, _state);
        for(i=0; i<=n-1; i++)
        {
            othererrors = othererrors||ae_fp_neq(x.ptr.p_double[i],xlast.ptr.p_double[i]);
        }
        n = 1;
        ae_vector_set_length(&x, n, _state);
        x.ptr.p_double[0] = 100;
        stpmax = 0.05+0.05*ae_randomreal(_state);
        mincgcreate(n, &x, &state, _state);
        mincgsetcond(&state, 1.0E-9, 0, 0, 0, _state);
        mincgsetstpmax(&state, stpmax, _state);
        mincgsetxrep(&state, ae_true, _state);
        xprev = x.ptr.p_double[0];
        while(mincgiteration(&state, _state))
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
    }
    
    /*
     * end
     */
    waserrors = (((((referror||eqerror)||linerror1)||linerror2)||converror)||othererrors)||restartserror;
    if( !silent )
    {
        printf("TESTING CG OPTIMIZATION\n");
        printf("REFERENCE PROBLEM:                        ");
        if( referror )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LIN-1 PROBLEM:                            ");
        if( linerror1 )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LIN-2 PROBLEM:                            ");
        if( linerror2 )
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
*************************************************************************/
static void testmincgunit_testfunc1(mincgstate* state, ae_state *_state)
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
static void testmincgunit_testfunc2(mincgstate* state, ae_state *_state)
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
static void testmincgunit_testfunc3(mincgstate* state, ae_state *_state)
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


/*$ End $*/
