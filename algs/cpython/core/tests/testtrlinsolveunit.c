

#include <stdafx.h>
#include <stdio.h>
#include "testtrlinsolveunit.h"


/*$ Declarations $*/
static void testtrlinsolveunit_makeacopy(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Main unittest subroutine
*************************************************************************/
ae_bool testtrlinsolve(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t maxmn;
    ae_int_t passcount;
    double threshold;
    ae_matrix aeffective;
    ae_matrix aparam;
    ae_vector xe;
    ae_vector b;
    ae_int_t n;
    ae_int_t pass;
    ae_int_t i;
    ae_int_t j;
    ae_int_t cnts;
    ae_int_t cntu;
    ae_int_t cntt;
    ae_int_t cntm;
    ae_bool waserrors;
    ae_bool isupper;
    ae_bool istrans;
    ae_bool isunit;
    double v;
    double s;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&aeffective, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&aparam, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xe, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);

    waserrors = ae_false;
    maxmn = 15;
    passcount = 15;
    threshold = 1000*ae_machineepsilon;
    
    /*
     * Different problems
     */
    for(n=1; n<=maxmn; n++)
    {
        ae_matrix_set_length(&aeffective, n-1+1, n-1+1, _state);
        ae_matrix_set_length(&aparam, n-1+1, n-1+1, _state);
        ae_vector_set_length(&xe, n-1+1, _state);
        ae_vector_set_length(&b, n-1+1, _state);
        for(pass=1; pass<=passcount; pass++)
        {
            for(cnts=0; cnts<=1; cnts++)
            {
                for(cntu=0; cntu<=1; cntu++)
                {
                    for(cntt=0; cntt<=1; cntt++)
                    {
                        for(cntm=0; cntm<=2; cntm++)
                        {
                            isupper = cnts==0;
                            isunit = cntu==0;
                            istrans = cntt==0;
                            
                            /*
                             * Skip meaningless combinations of parameters:
                             * (matrix is singular) AND (matrix is unit diagonal)
                             */
                            if( cntm==2&&isunit )
                            {
                                continue;
                            }
                            
                            /*
                             * Clear matrices
                             */
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    aeffective.ptr.pp_double[i][j] = 0;
                                    aparam.ptr.pp_double[i][j] = 0;
                                }
                            }
                            
                            /*
                             * Prepare matrices
                             */
                            if( isupper )
                            {
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=i; j<=n-1; j++)
                                    {
                                        aeffective.ptr.pp_double[i][j] = 0.9*(2*ae_randomreal(_state)-1);
                                        aparam.ptr.pp_double[i][j] = aeffective.ptr.pp_double[i][j];
                                    }
                                    aeffective.ptr.pp_double[i][i] = (2*ae_randominteger(2, _state)-1)*(0.8+ae_randomreal(_state));
                                    aparam.ptr.pp_double[i][i] = aeffective.ptr.pp_double[i][i];
                                }
                            }
                            else
                            {
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=i; j++)
                                    {
                                        aeffective.ptr.pp_double[i][j] = 0.9*(2*ae_randomreal(_state)-1);
                                        aparam.ptr.pp_double[i][j] = aeffective.ptr.pp_double[i][j];
                                    }
                                    aeffective.ptr.pp_double[i][i] = (2*ae_randominteger(2, _state)-1)*(0.8+ae_randomreal(_state));
                                    aparam.ptr.pp_double[i][i] = aeffective.ptr.pp_double[i][i];
                                }
                            }
                            if( isunit )
                            {
                                for(i=0; i<=n-1; i++)
                                {
                                    aeffective.ptr.pp_double[i][i] = 1;
                                    aparam.ptr.pp_double[i][i] = 0;
                                }
                            }
                            if( istrans )
                            {
                                if( isupper )
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        for(j=i+1; j<=n-1; j++)
                                        {
                                            aeffective.ptr.pp_double[j][i] = aeffective.ptr.pp_double[i][j];
                                            aeffective.ptr.pp_double[i][j] = 0;
                                        }
                                    }
                                }
                                else
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        for(j=i+1; j<=n-1; j++)
                                        {
                                            aeffective.ptr.pp_double[i][j] = aeffective.ptr.pp_double[j][i];
                                            aeffective.ptr.pp_double[j][i] = 0;
                                        }
                                    }
                                }
                            }
                            
                            /*
                             * Prepare task, solve, compare
                             */
                            for(i=0; i<=n-1; i++)
                            {
                                xe.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                v = ae_v_dotproduct(&aeffective.ptr.pp_double[i][0], 1, &xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
                                b.ptr.p_double[i] = v;
                            }
                            rmatrixtrsafesolve(&aparam, n, &b, &s, isupper, istrans, isunit, _state);
                            ae_v_muld(&xe.ptr.p_double[0], 1, ae_v_len(0,n-1), s);
                            ae_v_sub(&xe.ptr.p_double[0], 1, &b.ptr.p_double[0], 1, ae_v_len(0,n-1));
                            v = ae_v_dotproduct(&xe.ptr.p_double[0], 1, &xe.ptr.p_double[0], 1, ae_v_len(0,n-1));
                            v = ae_sqrt(v, _state);
                            waserrors = waserrors||ae_fp_greater(v,threshold);
                        }
                    }
                }
            }
        }
    }
    
    /*
     * report
     */
    if( !silent )
    {
        printf("TESTING RMatrixTRSafeSolve\n");
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
Copy
*************************************************************************/
static void testtrlinsolveunit_makeacopy(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(b);

    ae_matrix_set_length(b, m-1+1, n-1+1, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            b->ptr.pp_double[i][j] = a->ptr.pp_double[i][j];
        }
    }
}


/*$ End $*/
