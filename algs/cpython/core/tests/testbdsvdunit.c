

#include <stdafx.h>
#include <stdio.h>
#include "testbdsvdunit.h"


/*$ Declarations $*/
static void testbdsvdunit_fillidentity(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_state *_state);
static void testbdsvdunit_fillsparsede(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     double sparcity,
     ae_state *_state);
static void testbdsvdunit_getbdsvderror(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* u,
     /* Real    */ ae_matrix* c,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* vt,
     double* materr,
     double* orterr,
     ae_bool* wsorted,
     ae_state *_state);
static void testbdsvdunit_testbdsvdproblem(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     double* materr,
     double* orterr,
     ae_bool* wsorted,
     ae_bool* wfailed,
     ae_int_t* failcount,
     ae_int_t* succcount,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
ae_bool testbdsvd(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector d;
    ae_vector e;
    ae_matrix mempty;
    ae_int_t n;
    ae_int_t maxn;
    ae_int_t i;
    ae_int_t j;
    ae_int_t gpass;
    ae_int_t pass;
    ae_bool waserrors;
    ae_bool wsorted;
    ae_bool wfailed;
    ae_bool failcase;
    double materr;
    double orterr;
    double threshold;
    double failthreshold;
    double failr;
    ae_int_t failcount;
    ae_int_t succcount;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&e, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&mempty, 0, 0, DT_REAL, _state, ae_true);

    failcount = 0;
    succcount = 0;
    materr = 0;
    orterr = 0;
    wsorted = ae_true;
    wfailed = ae_false;
    waserrors = ae_false;
    maxn = 15;
    threshold = 5*100*ae_machineepsilon;
    failthreshold = 1.0E-2;
    ae_vector_set_length(&d, maxn-1+1, _state);
    ae_vector_set_length(&e, maxn-2+1, _state);
    
    /*
     * special case: fail matrix
     */
    n = 5;
    d.ptr.p_double[0] = -8.27448347422711894000e-01;
    d.ptr.p_double[1] = -8.16705832087160854600e-01;
    d.ptr.p_double[2] = -2.53974358904729382800e-17;
    d.ptr.p_double[3] = -1.24626684881972815700e+00;
    d.ptr.p_double[4] = -4.64744131545637651000e-01;
    e.ptr.p_double[0] = -3.25785088656270038800e-01;
    e.ptr.p_double[1] = -1.03732413708914436580e-01;
    e.ptr.p_double[2] = -9.57365642262031357700e-02;
    e.ptr.p_double[3] = -2.71564153973817390400e-01;
    failcase = rmatrixbdsvd(&d, &e, n, ae_true, ae_false, &mempty, 0, &mempty, 0, &mempty, 0, _state);
    
    /*
     * special case: zero divide matrix
     * unfixed LAPACK routine should fail on this problem
     */
    n = 7;
    d.ptr.p_double[0] = -6.96462904751731892700e-01;
    d.ptr.p_double[1] = 0.00000000000000000000e+00;
    d.ptr.p_double[2] = -5.73827770385971991400e-01;
    d.ptr.p_double[3] = -6.62562624399371191700e-01;
    d.ptr.p_double[4] = 5.82737148001782223600e-01;
    d.ptr.p_double[5] = 3.84825263580925003300e-01;
    d.ptr.p_double[6] = 9.84087420830525472200e-01;
    e.ptr.p_double[0] = -7.30307931760612871800e-02;
    e.ptr.p_double[1] = -2.30079042939542843800e-01;
    e.ptr.p_double[2] = -6.87824621739351216300e-01;
    e.ptr.p_double[3] = -1.77306437707837570600e-02;
    e.ptr.p_double[4] = 1.78285126526551632000e-15;
    e.ptr.p_double[5] = -4.89434737751289969400e-02;
    rmatrixbdsvd(&d, &e, n, ae_true, ae_false, &mempty, 0, &mempty, 0, &mempty, 0, _state);
    
    /*
     * zero matrix, several cases
     */
    for(i=0; i<=maxn-1; i++)
    {
        d.ptr.p_double[i] = 0;
    }
    for(i=0; i<=maxn-2; i++)
    {
        e.ptr.p_double[i] = 0;
    }
    for(n=1; n<=maxn; n++)
    {
        testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
    }
    
    /*
     * Dense matrix
     */
    for(n=1; n<=maxn; n++)
    {
        for(pass=1; pass<=10; pass++)
        {
            for(i=0; i<=maxn-1; i++)
            {
                d.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            for(i=0; i<=maxn-2; i++)
            {
                e.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
        }
    }
    
    /*
     * Sparse matrices, very sparse matrices, incredible sparse matrices
     */
    for(n=1; n<=maxn; n++)
    {
        for(pass=1; pass<=10; pass++)
        {
            testbdsvdunit_fillsparsede(&d, &e, n, 0.5, _state);
            testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
            testbdsvdunit_fillsparsede(&d, &e, n, 0.8, _state);
            testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
            testbdsvdunit_fillsparsede(&d, &e, n, 0.9, _state);
            testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
            testbdsvdunit_fillsparsede(&d, &e, n, 0.95, _state);
            testbdsvdunit_testbdsvdproblem(&d, &e, n, &materr, &orterr, &wsorted, &wfailed, &failcount, &succcount, _state);
        }
    }
    
    /*
     * report
     */
    failr = (double)failcount/(double)(succcount+failcount);
    waserrors = ((ae_fp_greater(materr,threshold)||ae_fp_greater(orterr,threshold))||!wsorted)||ae_fp_greater(failr,failthreshold);
    if( !silent )
    {
        printf("TESTING BIDIAGONAL SVD DECOMPOSITION\n");
        printf("SVD decomposition error:                 %5.3e\n",
            (double)(materr));
        printf("SVD orthogonality error:                 %5.3e\n",
            (double)(orterr));
        printf("Singular values order:                   ");
        if( wsorted )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("Always converged:                        ");
        if( !wfailed )
        {
            printf("YES\n");
        }
        else
        {
            printf("NO\n");
            printf("Fail ratio:                              %5.3f\n",
                (double)(failr));
        }
        printf("Fail matrix test:                        ");
        if( !failcase )
        {
            printf("AS EXPECTED\n");
        }
        else
        {
            printf("CONVERGED (UNEXPECTED)\n");
        }
        printf("Threshold:                               %5.3e\n",
            (double)(threshold));
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


static void testbdsvdunit_fillidentity(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    ae_matrix_set_length(a, n-1+1, n-1+1, _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( i==j )
            {
                a->ptr.pp_double[i][j] = 1;
            }
            else
            {
                a->ptr.pp_double[i][j] = 0;
            }
        }
    }
}


static void testbdsvdunit_fillsparsede(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     double sparcity,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    ae_vector_set_length(d, n-1+1, _state);
    ae_vector_set_length(e, ae_maxint(0, n-2, _state)+1, _state);
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_greater_eq(ae_randomreal(_state),sparcity) )
        {
            d->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        else
        {
            d->ptr.p_double[i] = 0;
        }
    }
    for(i=0; i<=n-2; i++)
    {
        if( ae_fp_greater_eq(ae_randomreal(_state),sparcity) )
        {
            e->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        else
        {
            e->ptr.p_double[i] = 0;
        }
    }
}


static void testbdsvdunit_getbdsvderror(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* u,
     /* Real    */ ae_matrix* c,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* vt,
     double* materr,
     double* orterr,
     ae_bool* wsorted,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double locerr;
    double sm;


    
    /*
     * decomposition error
     */
    locerr = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            sm = 0;
            for(k=0; k<=n-1; k++)
            {
                sm = sm+w->ptr.p_double[k]*u->ptr.pp_double[i][k]*vt->ptr.pp_double[k][j];
            }
            if( isupper )
            {
                if( i==j )
                {
                    locerr = ae_maxreal(locerr, ae_fabs(d->ptr.p_double[i]-sm, _state), _state);
                }
                else
                {
                    if( i==j-1 )
                    {
                        locerr = ae_maxreal(locerr, ae_fabs(e->ptr.p_double[i]-sm, _state), _state);
                    }
                    else
                    {
                        locerr = ae_maxreal(locerr, ae_fabs(sm, _state), _state);
                    }
                }
            }
            else
            {
                if( i==j )
                {
                    locerr = ae_maxreal(locerr, ae_fabs(d->ptr.p_double[i]-sm, _state), _state);
                }
                else
                {
                    if( i-1==j )
                    {
                        locerr = ae_maxreal(locerr, ae_fabs(e->ptr.p_double[j]-sm, _state), _state);
                    }
                    else
                    {
                        locerr = ae_maxreal(locerr, ae_fabs(sm, _state), _state);
                    }
                }
            }
        }
    }
    *materr = ae_maxreal(*materr, locerr, _state);
    
    /*
     * check for C = U'
     * we consider it as decomposition error
     */
    locerr = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            locerr = ae_maxreal(locerr, ae_fabs(u->ptr.pp_double[i][j]-c->ptr.pp_double[j][i], _state), _state);
        }
    }
    *materr = ae_maxreal(*materr, locerr, _state);
    
    /*
     * orthogonality error
     */
    locerr = 0;
    for(i=0; i<=n-1; i++)
    {
        for(j=i; j<=n-1; j++)
        {
            sm = ae_v_dotproduct(&u->ptr.pp_double[0][i], u->stride, &u->ptr.pp_double[0][j], u->stride, ae_v_len(0,n-1));
            if( i!=j )
            {
                locerr = ae_maxreal(locerr, ae_fabs(sm, _state), _state);
            }
            else
            {
                locerr = ae_maxreal(locerr, ae_fabs(sm-1, _state), _state);
            }
            sm = ae_v_dotproduct(&vt->ptr.pp_double[i][0], 1, &vt->ptr.pp_double[j][0], 1, ae_v_len(0,n-1));
            if( i!=j )
            {
                locerr = ae_maxreal(locerr, ae_fabs(sm, _state), _state);
            }
            else
            {
                locerr = ae_maxreal(locerr, ae_fabs(sm-1, _state), _state);
            }
        }
    }
    *orterr = ae_maxreal(*orterr, locerr, _state);
    
    /*
     * values order error
     */
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_greater(w->ptr.p_double[i],w->ptr.p_double[i-1]) )
        {
            *wsorted = ae_false;
        }
    }
}


static void testbdsvdunit_testbdsvdproblem(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     double* materr,
     double* orterr,
     ae_bool* wsorted,
     ae_bool* wfailed,
     ae_int_t* failcount,
     ae_int_t* succcount,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix u;
    ae_matrix vt;
    ae_matrix c;
    ae_vector w;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double v;
    double mx;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);

    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_greater(ae_fabs(d->ptr.p_double[i], _state),mx) )
        {
            mx = ae_fabs(d->ptr.p_double[i], _state);
        }
    }
    for(i=0; i<=n-2; i++)
    {
        if( ae_fp_greater(ae_fabs(e->ptr.p_double[i], _state),mx) )
        {
            mx = ae_fabs(e->ptr.p_double[i], _state);
        }
    }
    if( ae_fp_eq(mx,0) )
    {
        mx = 1;
    }
    
    /*
     * Upper BDSVD tests
     */
    ae_vector_set_length(&w, n-1+1, _state);
    testbdsvdunit_fillidentity(&u, n, _state);
    testbdsvdunit_fillidentity(&vt, n, _state);
    testbdsvdunit_fillidentity(&c, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = d->ptr.p_double[i];
    }
    if( !rmatrixbdsvd(&w, e, n, ae_true, ae_false, &u, n, &c, n, &vt, n, _state) )
    {
        *failcount = *failcount+1;
        *wfailed = ae_true;
        ae_frame_leave(_state);
        return;
    }
    testbdsvdunit_getbdsvderror(d, e, n, ae_true, &u, &c, &w, &vt, materr, orterr, wsorted, _state);
    testbdsvdunit_fillidentity(&u, n, _state);
    testbdsvdunit_fillidentity(&vt, n, _state);
    testbdsvdunit_fillidentity(&c, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = d->ptr.p_double[i];
    }
    if( !rmatrixbdsvd(&w, e, n, ae_true, ae_true, &u, n, &c, n, &vt, n, _state) )
    {
        *failcount = *failcount+1;
        *wfailed = ae_true;
        ae_frame_leave(_state);
        return;
    }
    testbdsvdunit_getbdsvderror(d, e, n, ae_true, &u, &c, &w, &vt, materr, orterr, wsorted, _state);
    
    /*
     * Lower BDSVD tests
     */
    ae_vector_set_length(&w, n-1+1, _state);
    testbdsvdunit_fillidentity(&u, n, _state);
    testbdsvdunit_fillidentity(&vt, n, _state);
    testbdsvdunit_fillidentity(&c, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = d->ptr.p_double[i];
    }
    if( !rmatrixbdsvd(&w, e, n, ae_false, ae_false, &u, n, &c, n, &vt, n, _state) )
    {
        *failcount = *failcount+1;
        *wfailed = ae_true;
        ae_frame_leave(_state);
        return;
    }
    testbdsvdunit_getbdsvderror(d, e, n, ae_false, &u, &c, &w, &vt, materr, orterr, wsorted, _state);
    testbdsvdunit_fillidentity(&u, n, _state);
    testbdsvdunit_fillidentity(&vt, n, _state);
    testbdsvdunit_fillidentity(&c, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = d->ptr.p_double[i];
    }
    if( !rmatrixbdsvd(&w, e, n, ae_false, ae_true, &u, n, &c, n, &vt, n, _state) )
    {
        *failcount = *failcount+1;
        *wfailed = ae_true;
        ae_frame_leave(_state);
        return;
    }
    testbdsvdunit_getbdsvderror(d, e, n, ae_false, &u, &c, &w, &vt, materr, orterr, wsorted, _state);
    
    /*
     * update counter
     */
    *succcount = *succcount+1;
    ae_frame_leave(_state);
}


/*$ End $*/
