

#include <stdafx.h>
#include <stdio.h>
#include "testminbleicunit.h"


/*$ Declarations $*/
static void testminbleicunit_checkbounds(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_int_t n,
     ae_bool* err,
     ae_state *_state);
static void testminbleicunit_testfeasibility(ae_bool* err,
     ae_state *_state);
static void testminbleicunit_testother(ae_bool* err, ae_state *_state);
static void testminbleicunit_testconv(ae_bool* err, ae_state *_state);


/*$ Body $*/


ae_bool testminbleic(ae_bool silent, ae_state *_state)
{
    ae_bool waserrors;
    ae_bool feasibilityerrors;
    ae_bool othererrors;
    ae_bool converrors;
    ae_bool result;


    waserrors = ae_false;
    feasibilityerrors = ae_false;
    othererrors = ae_false;
    converrors = ae_false;
    testminbleicunit_testfeasibility(&feasibilityerrors, _state);
    testminbleicunit_testother(&othererrors, _state);
    testminbleicunit_testconv(&converrors, _state);
    
    /*
     * end
     */
    waserrors = (feasibilityerrors||othererrors)||converrors;
    if( !silent )
    {
        printf("TESTING BLEIC OPTIMIZATION\n");
        printf("FEASIBILITY PROPERTIES:                   ");
        if( feasibilityerrors )
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
        printf("CONVERGENCE PROPERTIES:                   ");
        if( converrors )
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
    return result;
}


/*************************************************************************
Checks that X is bounded with respect to BndL/BndU.

If it is not, True is assigned to the Err variable (which is not changed
otherwise).
*************************************************************************/
static void testminbleicunit_checkbounds(/* Real    */ ae_vector* x,
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
This function test feasibility properties.
It launches a sequence of problems and examines their solutions.
Most of the attention is directed towards feasibility properties,
although we make some quick checks to ensure that actual solution is found.

On failure sets Err to True (leaves it unchanged otherwise)
*************************************************************************/
static void testminbleicunit_testfeasibility(ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t pkind;
    ae_int_t passcount;
    ae_int_t pass;
    ae_int_t n;
    ae_int_t nmax;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t p;
    double v;
    double v2;
    double v3;
    double vv;
    double muinit;
    ae_vector bl;
    ae_vector bu;
    ae_vector x;
    ae_vector g;
    ae_vector x0;
    ae_vector xs;
    ae_matrix c;
    ae_vector ct;
    minbleicstate state;
    double epsc;
    double epsg;
    minbleicreport rep;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&bl, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bu, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&g, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xs, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);
    _minbleicstate_init(&state, _state, ae_true);
    _minbleicreport_init(&rep, _state, ae_true);

    nmax = 5;
    epsc = 1.0E-4;
    epsg = 1.0E-8;
    muinit = 1.0E-3;
    passcount = 10;
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Test problem 1:
         * * no boundary and inequality constraints
         * * randomly generated plane as equality constraint
         * * random point (not necessarily on the plane)
         * * f = |x|^P, P = {2, 4} is used as target function
         * * we check that after work is over we are on the plane
         */
        for(pkind=1; pkind<=2; pkind++)
        {
            for(n=1; n<=nmax; n++)
            {
                
                /*
                 * Generate X, BL, BU, CT and left part of C.
                 *
                 * Right part of C is generated using somewhat complex algo:
                 * * we generate random vector and multiply it by C.
                 * * result is used as the right part.
                 * * calculations are done on the fly, vector itself is not stored
                 * We use such algo to be sure that our system is consistent.
                 */
                p = 2*pkind;
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&g, n, _state);
                ae_matrix_set_length(&c, 1, n+1, _state);
                ae_vector_set_length(&ct, 1, _state);
                c.ptr.pp_double[0][n] = 0;
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    c.ptr.pp_double[0][i] = 2*ae_randomreal(_state)-1;
                    v = 2*ae_randomreal(_state)-1;
                    c.ptr.pp_double[0][n] = c.ptr.pp_double[0][n]+c.ptr.pp_double[0][i]*v;
                }
                ct.ptr.p_int[0] = 0;
                
                /*
                 * Create and optimize
                 */
                minbleiccreate(n, &x, &state, _state);
                minbleicsetbarrierwidth(&state, muinit, _state);
                minbleicsetlc(&state, &c, &ct, 1, _state);
                minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                minbleicsetoutercond(&state, epsc, epsc, _state);
                while(minbleiciteration(&state, _state))
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+ae_pow(state.x.ptr.p_double[i], p, _state);
                            state.g.ptr.p_double[i] = p*ae_pow(state.x.ptr.p_double[i], p-1, _state);
                        }
                        continue;
                    }
                    
                    /*
                     * Unknown protocol specified
                     */
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                minbleicresults(&state, &x, &rep, _state);
                if( rep.terminationtype<=0 )
                {
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                
                /*
                 * Test feasibility of solution
                 */
                v = ae_v_dotproduct(&c.ptr.pp_double[0][0], 1, &x.ptr.p_double[0], 1, ae_v_len(0,n-1));
                *err = *err||ae_fp_greater(ae_fabs(v-c.ptr.pp_double[0][n], _state),epsc);
                
                /*
                 * if C is nonzero, test that result is
                 * a stationary point of constrained F.
                 *
                 * NOTE: this check is done only if C is nonzero
                 */
                vv = ae_v_dotproduct(&c.ptr.pp_double[0][0], 1, &c.ptr.pp_double[0][0], 1, ae_v_len(0,n-1));
                if( ae_fp_neq(vv,0) )
                {
                    
                    /*
                     * Calculate gradient at the result
                     * Project gradient into C
                     * Check projected norm
                     */
                    for(i=0; i<=n-1; i++)
                    {
                        g.ptr.p_double[i] = p*ae_pow(x.ptr.p_double[i], p-1, _state);
                    }
                    v2 = ae_v_dotproduct(&c.ptr.pp_double[0][0], 1, &c.ptr.pp_double[0][0], 1, ae_v_len(0,n-1));
                    v = ae_v_dotproduct(&c.ptr.pp_double[0][0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    vv = v/v2;
                    ae_v_subd(&g.ptr.p_double[0], 1, &c.ptr.pp_double[0][0], 1, ae_v_len(0,n-1), vv);
                    v3 = ae_v_dotproduct(&g.ptr.p_double[0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    *err = *err||ae_fp_greater(ae_sqrt(v3, _state),0.001);
                }
            }
        }
        
        /*
         * Test problem 2 (multiple equality constraints):
         * * 1<=N<=NMax, 1<=K<=N
         * * no boundary constraints
         * * N-dimensional space
         * * randomly generated point xs
         * * K randomly generated hyperplanes which all pass through xs
         *   define K equality constraints: (a[k],x)=b[k]
         * *f(x) = |x-x0|^2, x0 = xs+a[0]
         * * extremum of f(x) is exactly xs because:
         *   * xs is the closest point in the plane defined by (a[0],x)=b[0]
         *   * xs is feasible by definition
         */
        for(n=2; n<=nmax; n++)
        {
            for(k=1; k<=n; k++)
            {
                
                /*
                 * Generate X, X0, XS, BL, BU, CT and left part of C.
                 *
                 * Right part of C is generated using somewhat complex algo:
                 * * we generate random vector and multiply it by C.
                 * * result is used as the right part.
                 * * calculations are done on the fly, vector itself is not stored
                 * We use such algo to be sure that our system is consistent.
                 */
                p = 2*pkind;
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&x0, n, _state);
                ae_vector_set_length(&xs, n, _state);
                ae_vector_set_length(&g, n, _state);
                ae_matrix_set_length(&c, k, n+1, _state);
                ae_vector_set_length(&ct, k, _state);
                c.ptr.pp_double[0][n] = 0;
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    xs.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                }
                for(i=0; i<=k-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                    }
                    v = ae_v_dotproduct(&c.ptr.pp_double[i][0], 1, &xs.ptr.p_double[0], 1, ae_v_len(0,n-1));
                    c.ptr.pp_double[i][n] = v;
                    ct.ptr.p_int[i] = 0;
                }
                ae_v_move(&x0.ptr.p_double[0], 1, &xs.ptr.p_double[0], 1, ae_v_len(0,n-1));
                ae_v_add(&x0.ptr.p_double[0], 1, &c.ptr.pp_double[0][0], 1, ae_v_len(0,n-1));
                
                /*
                 * Create and optimize
                 */
                minbleiccreate(n, &x, &state, _state);
                minbleicsetbarrierwidth(&state, muinit, _state);
                minbleicsetlc(&state, &c, &ct, k, _state);
                minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                minbleicsetoutercond(&state, epsc, epsc, _state);
                while(minbleiciteration(&state, _state))
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+ae_sqr(state.x.ptr.p_double[i]-x0.ptr.p_double[i], _state);
                            state.g.ptr.p_double[i] = 2*(state.x.ptr.p_double[i]-x0.ptr.p_double[i]);
                        }
                        continue;
                    }
                    
                    /*
                     * Unknown protocol specified
                     */
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                minbleicresults(&state, &x, &rep, _state);
                if( rep.terminationtype<=0 )
                {
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                
                /*
                 * Compare with XS
                 */
                v = 0;
                for(i=0; i<=n-1; i++)
                {
                    v = v+ae_sqr(x.ptr.p_double[i]-xs.ptr.p_double[i], _state);
                }
                v = ae_sqrt(v, _state);
                *err = *err||ae_fp_greater(ae_fabs(v, _state),0.001);
            }
        }
        
        /*
         * Another simple problem:
         * * bound constraints 0 <= x[i] <= 1
         * * no linear constraints
         * * F(x) = |x-x0|^P, where P={2,4} and x0 is randomly selected from [-1,+2]^N
         * * with such simple boundaries and function it is easy to find
         *   analytic form of solution: S[i] = bound(x0[i], 0, 1).
         * * however, we can't guarantee that solution is strictly feasible
         *   with respect to nonlinearity constraint, so we check
         *   for approximate feasibility.
         */
        for(pkind=1; pkind<=2; pkind++)
        {
            for(n=1; n<=nmax; n++)
            {
                
                /*
                 * Generate X, BL, BU.
                 */
                p = 2*pkind;
                ae_vector_set_length(&bl, n, _state);
                ae_vector_set_length(&bu, n, _state);
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&x0, n, _state);
                for(i=0; i<=n-1; i++)
                {
                    bl.ptr.p_double[i] = 0;
                    bu.ptr.p_double[i] = 1;
                    x.ptr.p_double[i] = ae_randomreal(_state);
                    x0.ptr.p_double[i] = 3*ae_randomreal(_state)-1;
                }
                
                /*
                 * Create and optimize
                 */
                minbleiccreate(n, &x, &state, _state);
                minbleicsetbarrierwidth(&state, muinit, _state);
                minbleicsetbc(&state, &bl, &bu, _state);
                minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                minbleicsetoutercond(&state, epsc, epsc, _state);
                while(minbleiciteration(&state, _state))
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p, _state);
                            state.g.ptr.p_double[i] = p*ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p-1, _state);
                        }
                        continue;
                    }
                    
                    /*
                     * Unknown protocol specified
                     */
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                minbleicresults(&state, &x, &rep, _state);
                if( rep.terminationtype<=0 )
                {
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                
                /*
                 * * compare solution with analytic one
                 * * check feasibility
                 */
                for(i=0; i<=n-1; i++)
                {
                    *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-boundval(x0.ptr.p_double[i], 0.0, 1.0, _state), _state),0.01);
                    *err = *err||ae_fp_less(x.ptr.p_double[i],0.0-epsc);
                    *err = *err||ae_fp_greater(x.ptr.p_double[i],1.0+epsc);
                }
            }
        }
        
        /*
         * Same as previous one, but with bound constraints posed
         * as general linear ones:
         * * no bound constraints
         * * 2*N linear constraints 0 <= x[i] <= 1
         * * F(x) = |x-x0|^P, where P={2,4} and x0 is randomly selected from [-1,+2]^N
         * * with such simple constraints and function it is easy to find
         *   analytic form of solution: S[i] = bound(x0[i], 0, 1).
         * * however, we can't guarantee that solution is strictly feasible
         *   with respect to nonlinearity constraint, so we check
         *   for approximate feasibility.
         */
        for(pkind=1; pkind<=2; pkind++)
        {
            for(n=1; n<=nmax; n++)
            {
                
                /*
                 * Generate X, BL, BU.
                 */
                p = 2*pkind;
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&x0, n, _state);
                ae_matrix_set_length(&c, 2*n, n+1, _state);
                ae_vector_set_length(&ct, 2*n, _state);
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = ae_randomreal(_state);
                    x0.ptr.p_double[i] = 3*ae_randomreal(_state)-1;
                    for(j=0; j<=n; j++)
                    {
                        c.ptr.pp_double[2*i+0][j] = 0;
                        c.ptr.pp_double[2*i+1][j] = 0;
                    }
                    c.ptr.pp_double[2*i+0][i] = 1;
                    c.ptr.pp_double[2*i+0][n] = 0;
                    ct.ptr.p_int[2*i+0] = 1;
                    c.ptr.pp_double[2*i+1][i] = 1;
                    c.ptr.pp_double[2*i+1][n] = 1;
                    ct.ptr.p_int[2*i+1] = -1;
                }
                
                /*
                 * Create and optimize
                 */
                minbleiccreate(n, &x, &state, _state);
                minbleicsetbarrierwidth(&state, muinit, _state);
                minbleicsetlc(&state, &c, &ct, 2*n, _state);
                minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                minbleicsetoutercond(&state, epsc, epsc, _state);
                while(minbleiciteration(&state, _state))
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p, _state);
                            state.g.ptr.p_double[i] = p*ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p-1, _state);
                        }
                        continue;
                    }
                    
                    /*
                     * Unknown protocol specified
                     */
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                minbleicresults(&state, &x, &rep, _state);
                if( rep.terminationtype<=0 )
                {
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                
                /*
                 * * compare solution with analytic one
                 * * check feasibility
                 */
                for(i=0; i<=n-1; i++)
                {
                    *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[i]-boundval(x0.ptr.p_double[i], 0.0, 1.0, _state), _state),0.01);
                    *err = *err||ae_fp_less(x.ptr.p_double[i],0.0-epsc);
                    *err = *err||ae_fp_greater(x.ptr.p_double[i],1.0+epsc);
                }
            }
        }
        
        /*
         * Infeasible problem:
         * * all bound constraints are 0 <= x[i] <= 1 except for one
         * * that one is 0 >= x[i] >= 1
         * * no linear constraints
         * * F(x) = |x-x0|^P, where P={2,4} and x0 is randomly selected from [-1,+2]^N
         * * algorithm must return correct error code on such problem
         */
        for(pkind=1; pkind<=2; pkind++)
        {
            for(n=1; n<=nmax; n++)
            {
                
                /*
                 * Generate X, BL, BU.
                 */
                p = 2*pkind;
                ae_vector_set_length(&bl, n, _state);
                ae_vector_set_length(&bu, n, _state);
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&x0, n, _state);
                for(i=0; i<=n-1; i++)
                {
                    bl.ptr.p_double[i] = 0;
                    bu.ptr.p_double[i] = 1;
                    x.ptr.p_double[i] = ae_randomreal(_state);
                    x0.ptr.p_double[i] = 3*ae_randomreal(_state)-1;
                }
                i = ae_randominteger(n, _state);
                bl.ptr.p_double[i] = 1;
                bu.ptr.p_double[i] = 0;
                
                /*
                 * Create and optimize
                 */
                minbleiccreate(n, &x, &state, _state);
                minbleicsetbarrierwidth(&state, muinit, _state);
                minbleicsetbc(&state, &bl, &bu, _state);
                minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                minbleicsetoutercond(&state, epsc, epsc, _state);
                while(minbleiciteration(&state, _state))
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p, _state);
                            state.g.ptr.p_double[i] = p*ae_pow(state.x.ptr.p_double[i]-x0.ptr.p_double[i], p-1, _state);
                        }
                        continue;
                    }
                    
                    /*
                     * Unknown protocol specified
                     */
                    *err = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                minbleicresults(&state, &x, &rep, _state);
                *err = *err||rep.terminationtype!=-3;
            }
        }
        
        /*
         * Infeasible problem (2):
         * * no bound and inequality constraints
         * * 1<=K<=N arbitrary equality constraints
         * * (K+1)th constraint which is equal to the first constraint a*x=c,
         *   but with c:=c+1. I.e. we have both a*x=c and a*x=c+1, which can't
         *   be true (other constraints may be inconsistent too, but we don't
         *   have to check it).
         * * F(x) = |x|^P, where P={2,4}
         * * algorithm must return correct error code on such problem
         */
        for(pkind=1; pkind<=2; pkind++)
        {
            for(n=1; n<=nmax; n++)
            {
                for(k=1; k<=n; k++)
                {
                    
                    /*
                     * Generate X, BL, BU.
                     */
                    p = 2*pkind;
                    ae_vector_set_length(&x, n, _state);
                    ae_matrix_set_length(&c, k+1, n+1, _state);
                    ae_vector_set_length(&ct, k+1, _state);
                    for(i=0; i<=n-1; i++)
                    {
                        x.ptr.p_double[i] = ae_randomreal(_state);
                    }
                    for(i=0; i<=k-1; i++)
                    {
                        for(j=0; j<=n; j++)
                        {
                            c.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                        }
                        ct.ptr.p_int[i] = 0;
                    }
                    ct.ptr.p_int[k] = 0;
                    ae_v_move(&c.ptr.pp_double[k][0], 1, &c.ptr.pp_double[0][0], 1, ae_v_len(0,n-1));
                    c.ptr.pp_double[k][n] = c.ptr.pp_double[0][n]+1;
                    
                    /*
                     * Create and optimize
                     */
                    minbleiccreate(n, &x, &state, _state);
                    minbleicsetbarrierwidth(&state, muinit, _state);
                    minbleicsetlc(&state, &c, &ct, k+1, _state);
                    minbleicsetinnercond(&state, epsg, 0.0, 0.0, _state);
                    minbleicsetoutercond(&state, epsc, epsc, _state);
                    while(minbleiciteration(&state, _state))
                    {
                        if( state.needfg )
                        {
                            state.f = 0;
                            for(i=0; i<=n-1; i++)
                            {
                                state.f = state.f+ae_pow(state.x.ptr.p_double[i], p, _state);
                                state.g.ptr.p_double[i] = p*ae_pow(state.x.ptr.p_double[i], p-1, _state);
                            }
                            continue;
                        }
                        
                        /*
                         * Unknown protocol specified
                         */
                        *err = ae_true;
                        ae_frame_leave(_state);
                        return;
                    }
                    minbleicresults(&state, &x, &rep, _state);
                    *err = *err||rep.terminationtype!=-3;
                }
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function additional properties.

On failure sets Err to True (leaves it unchanged otherwise)
*************************************************************************/
static void testminbleicunit_testother(ae_bool* err, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t passcount;
    ae_int_t pass;
    ae_int_t n;
    ae_int_t nmax;
    ae_int_t i;
    ae_vector bl;
    ae_vector bu;
    ae_vector x;
    ae_vector xf;
    ae_vector xlast;
    ae_matrix c;
    ae_vector ct;
    double fprev;
    double xprev;
    double stpmax;
    double mu;
    double muinit;
    minbleicstate state;
    double epsc;
    double epsg;
    minbleicreport rep;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&bl, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bu, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xlast, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);
    _minbleicstate_init(&state, _state, ae_true);
    _minbleicreport_init(&rep, _state, ae_true);

    nmax = 5;
    epsc = 1.0E-4;
    epsg = 1.0E-8;
    muinit = 1.0E-3;
    passcount = 10;
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Test reports:
         * * first value must be starting point
         * * last value must be last point
         */
        n = 50;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xlast, n, _state);
        ae_vector_set_length(&bl, n, _state);
        ae_vector_set_length(&bu, n, _state);
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = 10;
            bl.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            bu.ptr.p_double[i] = _state->v_posinf;
        }
        minbleiccreate(n, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetinnercond(&state, 0, 0, 0, _state);
        minbleicsetmaxits(&state, 10, _state);
        minbleicsetoutercond(&state, 1.0E-64, 1.0E-64, _state);
        minbleicsetxrep(&state, ae_true, _state);
        fprev = ae_maxrealnumber;
        while(minbleiciteration(&state, _state))
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
                if( ae_fp_eq(fprev,ae_maxrealnumber) )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        *err = *err||ae_fp_neq(state.x.ptr.p_double[i],x.ptr.p_double[i]);
                    }
                }
                fprev = state.f;
                ae_v_move(&xlast.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,n-1));
            }
        }
        minbleicresults(&state, &x, &rep, _state);
        for(i=0; i<=n-1; i++)
        {
            *err = *err||ae_fp_neq(x.ptr.p_double[i],xlast.ptr.p_double[i]);
        }
        
        /*
         * Test stpmax
         */
        n = 1;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&bl, n, _state);
        ae_vector_set_length(&bu, n, _state);
        x.ptr.p_double[0] = 100;
        bl.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        bu.ptr.p_double[0] = _state->v_posinf;
        stpmax = 0.05+0.05*ae_randomreal(_state);
        minbleiccreate(n, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, epsc, epsc, _state);
        minbleicsetxrep(&state, ae_true, _state);
        minbleicsetstpmax(&state, stpmax, _state);
        xprev = x.ptr.p_double[0];
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = ae_exp(state.x.ptr.p_double[0], _state)+ae_exp(-state.x.ptr.p_double[0], _state);
                state.g.ptr.p_double[0] = ae_exp(state.x.ptr.p_double[0], _state)-ae_exp(-state.x.ptr.p_double[0], _state);
                *err = *err||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]-xprev, _state),(1+ae_sqrt(ae_machineepsilon, _state))*stpmax);
            }
            if( state.xupdated )
            {
                *err = *err||ae_fp_greater(ae_fabs(state.x.ptr.p_double[0]-xprev, _state),(1+ae_sqrt(ae_machineepsilon, _state))*stpmax);
                xprev = state.x.ptr.p_double[0];
            }
        }
        
        /*
         * Ability to solve problems with function which is unbounded from below
         */
        n = 1;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&bl, n, _state);
        ae_vector_set_length(&bu, n, _state);
        bl.ptr.p_double[0] = 4*ae_randomreal(_state)+1;
        bu.ptr.p_double[0] = bl.ptr.p_double[0]+1;
        x.ptr.p_double[0] = 0.5*(bl.ptr.p_double[0]+bu.ptr.p_double[0]);
        minbleiccreate(n, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, epsc, epsc, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = -1.0E8*ae_sqr(state.x.ptr.p_double[0], _state);
                state.g.ptr.p_double[0] = -2.0E8*state.x.ptr.p_double[0];
            }
        }
        minbleicresults(&state, &x, &rep, _state);
        *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[0]-bu.ptr.p_double[0], _state),epsc);
        
        /*
         * * constant barrier width is not changed during outer iteration
         * * decaying barrier width is decreased at least once (if very stringent
         *   stopping conditions are specified)
         */
        n = 1;
        mu = 1.0E-2;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&xf, n, _state);
        ae_vector_set_length(&bl, n, _state);
        ae_vector_set_length(&bu, n, _state);
        bl.ptr.p_double[0] = 4*ae_randomreal(_state)+1;
        bu.ptr.p_double[0] = bl.ptr.p_double[0]+1;
        x.ptr.p_double[0] = 0.5*(bl.ptr.p_double[0]+bu.ptr.p_double[0]);
        minbleiccreate(n, &x, &state, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, 1.0E-10, 1.0E-10, _state);
        minbleicsetbarrierwidth(&state, mu, _state);
        minbleicsetbarrierdecay(&state, 1.0, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = 1.0E3*ae_sqr(state.x.ptr.p_double[0], _state);
                state.g.ptr.p_double[0] = 2.0E3*state.x.ptr.p_double[0];
            }
        }
        minbleicresults(&state, &xf, &rep, _state);
        *err = *err||rep.terminationtype<=0;
        *err = *err||ae_fp_neq(state.mu,mu);
        minbleicrestartfrom(&state, &x, _state);
        minbleicsetbarrierdecay(&state, 0.5, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = 1.0E6*ae_sqr(state.x.ptr.p_double[0], _state);
                state.g.ptr.p_double[0] = 2.0E6*state.x.ptr.p_double[0];
            }
        }
        minbleicresults(&state, &xf, &rep, _state);
        *err = *err||rep.terminationtype<=0;
        *err = *err||ae_fp_greater_eq(state.mu,mu);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function tests convergence properties.
We solve several simple problems with different combinations of constraints

On failure sets Err to True (leaves it unchanged otherwise)
*************************************************************************/
static void testminbleicunit_testconv(ae_bool* err, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t passcount;
    ae_int_t pass;
    ae_int_t n;
    ae_vector bl;
    ae_vector bu;
    ae_vector x;
    ae_matrix c;
    ae_vector ct;
    double muinit;
    minbleicstate state;
    double epsc;
    double epsg;
    double tol;
    minbleicreport rep;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&bl, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bu, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);
    _minbleicstate_init(&state, _state, ae_true);
    _minbleicreport_init(&rep, _state, ae_true);

    epsc = 1.0E-4;
    epsg = 1.0E-8;
    muinit = 1.0E-3;
    tol = 0.001;
    passcount = 10;
    
    /*
     * Three closely connected problems:
     * * 2-dimensional space
     * * octagonal area bounded by:
     *   * -1<=x<=+1
     *   * -1<=y<=+1
     *   * x+y<=1.5
     *   * x-y<=1.5
     *   * -x+y<=1.5
     *   * -x-y<=1.5
     * * several target functions:
     *   * f0=x+0.001*y, minimum at x=-1, y=-0.5
     *   * f1=(x+10)^2+y^2, minimum at x=-1, y=0
     *   * f2=(x+10)^2+(y-0.6)^2, minimum at x=-1, y=0.5
     */
    ae_vector_set_length(&x, 2, _state);
    ae_vector_set_length(&bl, 2, _state);
    ae_vector_set_length(&bu, 2, _state);
    ae_matrix_set_length(&c, 4, 3, _state);
    ae_vector_set_length(&ct, 4, _state);
    bl.ptr.p_double[0] = -1;
    bl.ptr.p_double[1] = -1;
    bu.ptr.p_double[0] = 1;
    bu.ptr.p_double[1] = 1;
    c.ptr.pp_double[0][0] = 1;
    c.ptr.pp_double[0][1] = 1;
    c.ptr.pp_double[0][2] = 1.5;
    ct.ptr.p_int[0] = -1;
    c.ptr.pp_double[1][0] = 1;
    c.ptr.pp_double[1][1] = -1;
    c.ptr.pp_double[1][2] = 1.5;
    ct.ptr.p_int[1] = -1;
    c.ptr.pp_double[2][0] = -1;
    c.ptr.pp_double[2][1] = 1;
    c.ptr.pp_double[2][2] = 1.5;
    ct.ptr.p_int[2] = -1;
    c.ptr.pp_double[3][0] = -1;
    c.ptr.pp_double[3][1] = -1;
    c.ptr.pp_double[3][2] = 1.5;
    ct.ptr.p_int[3] = -1;
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * f0
         */
        x.ptr.p_double[0] = 0.2*ae_randomreal(_state)-0.1;
        x.ptr.p_double[1] = 0.2*ae_randomreal(_state)-0.1;
        minbleiccreate(2, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetlc(&state, &c, &ct, 4, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, epsc, epsc, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = state.x.ptr.p_double[0]+0.001*state.x.ptr.p_double[1];
                state.g.ptr.p_double[0] = 1;
                state.g.ptr.p_double[1] = 0.001;
            }
        }
        minbleicresults(&state, &x, &rep, _state);
        if( rep.terminationtype>0 )
        {
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[0]+1, _state),tol);
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[1]+0.5, _state),tol);
        }
        else
        {
            *err = ae_true;
        }
        
        /*
         * f1
         */
        x.ptr.p_double[0] = 0.2*ae_randomreal(_state)-0.1;
        x.ptr.p_double[1] = 0.2*ae_randomreal(_state)-0.1;
        minbleiccreate(2, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetlc(&state, &c, &ct, 4, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, epsc, epsc, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = ae_sqr(state.x.ptr.p_double[0]+10, _state)+ae_sqr(state.x.ptr.p_double[1], _state);
                state.g.ptr.p_double[0] = 2*(state.x.ptr.p_double[0]+10);
                state.g.ptr.p_double[1] = 2*state.x.ptr.p_double[1];
            }
        }
        minbleicresults(&state, &x, &rep, _state);
        if( rep.terminationtype>0 )
        {
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[0]+1, _state),tol);
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[1], _state),tol);
        }
        else
        {
            *err = ae_true;
        }
        
        /*
         * f2
         */
        x.ptr.p_double[0] = 0.2*ae_randomreal(_state)-0.1;
        x.ptr.p_double[1] = 0.2*ae_randomreal(_state)-0.1;
        minbleiccreate(2, &x, &state, _state);
        minbleicsetbarrierwidth(&state, muinit, _state);
        minbleicsetbc(&state, &bl, &bu, _state);
        minbleicsetlc(&state, &c, &ct, 4, _state);
        minbleicsetinnercond(&state, epsg, 0, 0, _state);
        minbleicsetoutercond(&state, epsc, epsc, _state);
        while(minbleiciteration(&state, _state))
        {
            if( state.needfg )
            {
                state.f = ae_sqr(state.x.ptr.p_double[0]+10, _state)+ae_sqr(state.x.ptr.p_double[1]-0.6, _state);
                state.g.ptr.p_double[0] = 2*(state.x.ptr.p_double[0]+10);
                state.g.ptr.p_double[1] = 2*(state.x.ptr.p_double[1]-0.6);
            }
        }
        minbleicresults(&state, &x, &rep, _state);
        if( rep.terminationtype>0 )
        {
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[0]+1, _state),tol);
            *err = *err||ae_fp_greater(ae_fabs(x.ptr.p_double[1]-0.5, _state),tol);
        }
        else
        {
            *err = ae_true;
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
