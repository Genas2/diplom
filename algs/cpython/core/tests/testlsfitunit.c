

#include <stdafx.h>
#include <stdio.h>
#include "testlsfitunit.h"


/*$ Declarations $*/
static void testlsfitunit_testpolynomialfitting(ae_bool* fiterrors,
     ae_state *_state);
static void testlsfitunit_testrationalfitting(ae_bool* fiterrors,
     ae_state *_state);
static void testlsfitunit_testsplinefitting(ae_bool* fiterrors,
     ae_state *_state);
static void testlsfitunit_testgeneralfitting(ae_bool* llserrors,
     ae_bool* nlserrors,
     ae_state *_state);
static ae_bool testlsfitunit_isglssolution(ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     /* Real    */ ae_vector* c,
     ae_state *_state);
static double testlsfitunit_getglserror(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_vector* c,
     ae_state *_state);
static void testlsfitunit_fitlinearnonlinear(ae_int_t m,
     ae_int_t deravailable,
     /* Real    */ ae_matrix* xy,
     lsfitstate* state,
     ae_bool* nlserrors,
     ae_state *_state);


/*$ Body $*/


ae_bool testlsfit(ae_bool silent, ae_state *_state)
{
    ae_bool waserrors;
    ae_bool llserrors;
    ae_bool nlserrors;
    ae_bool polfiterrors;
    ae_bool ratfiterrors;
    ae_bool splfiterrors;
    ae_bool result;


    waserrors = ae_false;
    testlsfitunit_testpolynomialfitting(&polfiterrors, _state);
    testlsfitunit_testrationalfitting(&ratfiterrors, _state);
    testlsfitunit_testsplinefitting(&splfiterrors, _state);
    testlsfitunit_testgeneralfitting(&llserrors, &nlserrors, _state);
    
    /*
     * report
     */
    waserrors = (((llserrors||nlserrors)||polfiterrors)||ratfiterrors)||splfiterrors;
    if( !silent )
    {
        printf("TESTING LEAST SQUARES\n");
        printf("POLYNOMIAL LEAST SQUARES:                ");
        if( polfiterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("RATIONAL LEAST SQUARES:                  ");
        if( ratfiterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("SPLINE LEAST SQUARES:                    ");
        if( splfiterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("LINEAR LEAST SQUARES:                    ");
        if( llserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("NON-LINEAR LEAST SQUARES:                ");
        if( nlserrors )
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
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
static void testlsfitunit_testpolynomialfitting(ae_bool* fiterrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    double threshold;
    ae_vector x;
    ae_vector y;
    ae_vector w;
    ae_vector x2;
    ae_vector y2;
    ae_vector w2;
    ae_vector xfull;
    ae_vector yfull;
    double a;
    double b;
    double t;
    ae_int_t i;
    ae_int_t k;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;
    ae_int_t info;
    ae_int_t info2;
    double v;
    double v0;
    double v1;
    double v2;
    double s;
    double xmin;
    double xmax;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    barycentricinterpolant p;
    barycentricinterpolant p1;
    barycentricinterpolant p2;
    polynomialfitreport rep;
    polynomialfitreport rep2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t maxn;
    ae_int_t pass;
    ae_int_t passcount;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xfull, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yfull, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);
    _barycentricinterpolant_init(&p, _state, ae_true);
    _barycentricinterpolant_init(&p1, _state, ae_true);
    _barycentricinterpolant_init(&p2, _state, ae_true);
    _polynomialfitreport_init(&rep, _state, ae_true);
    _polynomialfitreport_init(&rep2, _state, ae_true);

    *fiterrors = ae_false;
    maxn = 5;
    passcount = 20;
    threshold = 1.0E8*ae_machineepsilon;
    
    /*
     * Test polunomial fitting
     */
    for(pass=1; pass<=passcount; pass++)
    {
        for(n=1; n<=maxn; n++)
        {
            
            /*
             * N=M+K fitting (i.e. interpolation)
             */
            for(k=0; k<=n-1; k++)
            {
                taskgenint1d(-1, 1, n, &xfull, &yfull, _state);
                ae_vector_set_length(&x, n-k, _state);
                ae_vector_set_length(&y, n-k, _state);
                ae_vector_set_length(&w, n-k, _state);
                if( k>0 )
                {
                    ae_vector_set_length(&xc, k, _state);
                    ae_vector_set_length(&yc, k, _state);
                    ae_vector_set_length(&dc, k, _state);
                }
                for(i=0; i<=n-k-1; i++)
                {
                    x.ptr.p_double[i] = xfull.ptr.p_double[i];
                    y.ptr.p_double[i] = yfull.ptr.p_double[i];
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                }
                for(i=0; i<=k-1; i++)
                {
                    xc.ptr.p_double[i] = xfull.ptr.p_double[n-k+i];
                    yc.ptr.p_double[i] = yfull.ptr.p_double[n-k+i];
                    dc.ptr.p_int[i] = 0;
                }
                polynomialfitwc(&x, &y, &w, n-k, &xc, &yc, &dc, k, n, &info, &p1, &rep, _state);
                if( info<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    for(i=0; i<=n-k-1; i++)
                    {
                        *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(barycentriccalc(&p1, x.ptr.p_double[i], _state)-y.ptr.p_double[i], _state),threshold);
                    }
                    for(i=0; i<=k-1; i++)
                    {
                        *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(barycentriccalc(&p1, xc.ptr.p_double[i], _state)-yc.ptr.p_double[i], _state),threshold);
                    }
                }
            }
            
            /*
             * Testing constraints on derivatives.
             * Special tasks which will always have solution:
             * 1. P(0)=YC[0]
             * 2. P(0)=YC[0], P'(0)=YC[1]
             */
            if( n>1 )
            {
                for(m=3; m<=5; m++)
                {
                    for(k=1; k<=2; k++)
                    {
                        taskgenint1d(-1, 1, n, &x, &y, _state);
                        ae_vector_set_length(&w, n, _state);
                        ae_vector_set_length(&xc, 2, _state);
                        ae_vector_set_length(&yc, 2, _state);
                        ae_vector_set_length(&dc, 2, _state);
                        for(i=0; i<=n-1; i++)
                        {
                            w.ptr.p_double[i] = 1+ae_randomreal(_state);
                        }
                        xc.ptr.p_double[0] = 0;
                        yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
                        dc.ptr.p_int[0] = 0;
                        xc.ptr.p_double[1] = 0;
                        yc.ptr.p_double[1] = 2*ae_randomreal(_state)-1;
                        dc.ptr.p_int[1] = 1;
                        polynomialfitwc(&x, &y, &w, n, &xc, &yc, &dc, k, m, &info, &p1, &rep, _state);
                        if( info<=0 )
                        {
                            *fiterrors = ae_true;
                        }
                        else
                        {
                            barycentricdiff1(&p1, 0.0, &v0, &v1, _state);
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v0-yc.ptr.p_double[0], _state),threshold);
                            if( k==2 )
                            {
                                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v1-yc.ptr.p_double[1], _state),threshold);
                            }
                        }
                    }
                }
            }
        }
    }
    for(m=2; m<=8; m++)
    {
        for(pass=1; pass<=passcount; pass++)
        {
            
            /*
             * General fitting
             *
             * interpolating function through M nodes should have
             * greater RMS error than fitting it through the same M nodes
             */
            n = 100;
            ae_vector_set_length(&x2, n, _state);
            ae_vector_set_length(&y2, n, _state);
            ae_vector_set_length(&w2, n, _state);
            xmin = 0;
            xmax = 2*ae_pi;
            for(i=0; i<=n-1; i++)
            {
                x2.ptr.p_double[i] = 2*ae_pi*ae_randomreal(_state);
                y2.ptr.p_double[i] = ae_sin(x2.ptr.p_double[i], _state);
                w2.ptr.p_double[i] = 1;
            }
            ae_vector_set_length(&x, m, _state);
            ae_vector_set_length(&y, m, _state);
            for(i=0; i<=m-1; i++)
            {
                x.ptr.p_double[i] = xmin+(xmax-xmin)*i/(m-1);
                y.ptr.p_double[i] = ae_sin(x.ptr.p_double[i], _state);
            }
            polynomialbuild(&x, &y, m, &p1, _state);
            polynomialfitwc(&x2, &y2, &w2, n, &xc, &yc, &dc, 0, m, &info, &p2, &rep, _state);
            if( info<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * calculate P1 (interpolant) RMS error, compare with P2 error
                 */
                v1 = 0;
                v2 = 0;
                for(i=0; i<=n-1; i++)
                {
                    v1 = v1+ae_sqr(barycentriccalc(&p1, x2.ptr.p_double[i], _state)-y2.ptr.p_double[i], _state);
                    v2 = v2+ae_sqr(barycentriccalc(&p2, x2.ptr.p_double[i], _state)-y2.ptr.p_double[i], _state);
                }
                v1 = ae_sqrt(v1/n, _state);
                v2 = ae_sqrt(v2/n, _state);
                *fiterrors = *fiterrors||ae_fp_greater(v2,v1);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v2-rep.rmserror, _state),threshold);
            }
            
            /*
             * compare weighted and non-weighted
             */
            n = 20;
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                w.ptr.p_double[i] = 1;
            }
            polynomialfitwc(&x, &y, &w, n, &xc, &yc, &dc, 0, m, &info, &p1, &rep, _state);
            polynomialfit(&x, &y, n, m, &info2, &p2, &rep2, _state);
            if( info<=0||info2<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * calculate P1 (interpolant), compare with P2 error
                 * compare RMS errors
                 */
                t = 2*ae_randomreal(_state)-1;
                v1 = barycentriccalc(&p1, t, _state);
                v2 = barycentriccalc(&p2, t, _state);
                *fiterrors = *fiterrors||ae_fp_neq(v2,v1);
                *fiterrors = *fiterrors||ae_fp_neq(rep.rmserror,rep2.rmserror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.avgerror,rep2.avgerror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.avgrelerror,rep2.avgrelerror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.maxerror,rep2.maxerror);
            }
        }
    }
    for(m=1; m<=maxn; m++)
    {
        for(pass=1; pass<=passcount; pass++)
        {
            ae_assert(passcount>=2, "PassCount should be 2 or greater!", _state);
            
            /*
             * solve simple task (all X[] are the same, Y[] are specially
             * calculated to ensure simple form of all types of errors)
             * and check correctness of the errors calculated by subroutines
             *
             * First pass is done with zero Y[], other passes - with random Y[].
             * It should test both ability to correctly calculate errors and
             * ability to not fail while working with zeros :)
             */
            n = 4*maxn;
            if( pass==1 )
            {
                v1 = 0;
                v2 = 0;
                v = 0;
            }
            else
            {
                v1 = ae_randomreal(_state);
                v2 = ae_randomreal(_state);
                v = 1+ae_randomreal(_state);
            }
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            for(i=0; i<=maxn-1; i++)
            {
                x.ptr.p_double[4*i+0] = i;
                y.ptr.p_double[4*i+0] = v-v2;
                w.ptr.p_double[4*i+0] = 1;
                x.ptr.p_double[4*i+1] = i;
                y.ptr.p_double[4*i+1] = v-v1;
                w.ptr.p_double[4*i+1] = 1;
                x.ptr.p_double[4*i+2] = i;
                y.ptr.p_double[4*i+2] = v+v1;
                w.ptr.p_double[4*i+2] = 1;
                x.ptr.p_double[4*i+3] = i;
                y.ptr.p_double[4*i+3] = v+v2;
                w.ptr.p_double[4*i+3] = 1;
            }
            refrms = ae_sqrt((ae_sqr(v1, _state)+ae_sqr(v2, _state))/2, _state);
            refavg = (ae_fabs(v1, _state)+ae_fabs(v2, _state))/2;
            if( pass==1 )
            {
                refavgrel = 0;
            }
            else
            {
                refavgrel = 0.25*(ae_fabs(v2, _state)/ae_fabs(v-v2, _state)+ae_fabs(v1, _state)/ae_fabs(v-v1, _state)+ae_fabs(v1, _state)/ae_fabs(v+v1, _state)+ae_fabs(v2, _state)/ae_fabs(v+v2, _state));
            }
            refmax = ae_maxreal(v1, v2, _state);
            
            /*
             * Test errors correctness
             */
            polynomialfit(&x, &y, n, m, &info, &p, &rep, _state);
            if( info<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                s = barycentriccalc(&p, 0, _state);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-v, _state),threshold);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
            }
        }
    }
    ae_frame_leave(_state);
}


static void testlsfitunit_testrationalfitting(ae_bool* fiterrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    double threshold;
    ae_int_t maxn;
    ae_int_t passcount;
    barycentricinterpolant b1;
    barycentricinterpolant b2;
    ae_vector x;
    ae_vector x2;
    ae_vector y;
    ae_vector y2;
    ae_vector w;
    ae_vector w2;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;
    double h;
    double s1;
    double s2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t n2;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t d;
    ae_int_t pass;
    double err;
    double maxerr;
    double t;
    double a;
    double b;
    double s;
    double v;
    double v0;
    double v1;
    double v2;
    double v3;
    double d0;
    double d1;
    double d2;
    ae_int_t info;
    ae_int_t info2;
    double xmin;
    double xmax;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    barycentricfitreport rep;
    barycentricfitreport rep2;

    ae_frame_make(_state, &_frame_block);
    _barycentricinterpolant_init(&b1, _state, ae_true);
    _barycentricinterpolant_init(&b2, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);
    _barycentricfitreport_init(&rep, _state, ae_true);
    _barycentricfitreport_init(&rep2, _state, ae_true);

    *fiterrors = ae_false;
    
    /*
     * PassCount        number of repeated passes
     * Threshold        error tolerance
     * LipschitzTol     Lipschitz constant increase allowed
     *                  when calculating constant on a twice denser grid
     */
    passcount = 5;
    maxn = 15;
    threshold = 1000000*ae_machineepsilon;
    
    /*
     * Test rational fitting:
     */
    for(pass=1; pass<=passcount; pass++)
    {
        for(n=2; n<=maxn; n++)
        {
            
            /*
             * N=M+K fitting (i.e. interpolation)
             */
            for(k=0; k<=n-1; k++)
            {
                ae_vector_set_length(&x, n-k, _state);
                ae_vector_set_length(&y, n-k, _state);
                ae_vector_set_length(&w, n-k, _state);
                if( k>0 )
                {
                    ae_vector_set_length(&xc, k, _state);
                    ae_vector_set_length(&yc, k, _state);
                    ae_vector_set_length(&dc, k, _state);
                }
                for(i=0; i<=n-k-1; i++)
                {
                    x.ptr.p_double[i] = (double)i/(double)(n-1);
                    y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                }
                for(i=0; i<=k-1; i++)
                {
                    xc.ptr.p_double[i] = (double)(n-k+i)/(double)(n-1);
                    yc.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    dc.ptr.p_int[i] = 0;
                }
                barycentricfitfloaterhormannwc(&x, &y, &w, n-k, &xc, &yc, &dc, k, n, &info, &b1, &rep, _state);
                if( info<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    for(i=0; i<=n-k-1; i++)
                    {
                        *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(barycentriccalc(&b1, x.ptr.p_double[i], _state)-y.ptr.p_double[i], _state),threshold);
                    }
                    for(i=0; i<=k-1; i++)
                    {
                        *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(barycentriccalc(&b1, xc.ptr.p_double[i], _state)-yc.ptr.p_double[i], _state),threshold);
                    }
                }
            }
            
            /*
             * Testing constraints on derivatives:
             * * several M's are tried
             * * several K's are tried - 1, 2.
             * * constraints at the ends of the interval
             */
            for(m=3; m<=5; m++)
            {
                for(k=1; k<=2; k++)
                {
                    ae_vector_set_length(&x, n, _state);
                    ae_vector_set_length(&y, n, _state);
                    ae_vector_set_length(&w, n, _state);
                    ae_vector_set_length(&xc, 2, _state);
                    ae_vector_set_length(&yc, 2, _state);
                    ae_vector_set_length(&dc, 2, _state);
                    for(i=0; i<=n-1; i++)
                    {
                        x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                        y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                        w.ptr.p_double[i] = 1+ae_randomreal(_state);
                    }
                    xc.ptr.p_double[0] = -1;
                    yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
                    dc.ptr.p_int[0] = 0;
                    xc.ptr.p_double[1] = 1;
                    yc.ptr.p_double[1] = 2*ae_randomreal(_state)-1;
                    dc.ptr.p_int[1] = 0;
                    barycentricfitfloaterhormannwc(&x, &y, &w, n, &xc, &yc, &dc, k, m, &info, &b1, &rep, _state);
                    if( info<=0 )
                    {
                        *fiterrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=k-1; i++)
                        {
                            barycentricdiff1(&b1, xc.ptr.p_double[i], &v0, &v1, _state);
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v0-yc.ptr.p_double[i], _state),threshold);
                        }
                    }
                }
            }
        }
    }
    for(m=2; m<=8; m++)
    {
        for(pass=1; pass<=passcount; pass++)
        {
            
            /*
             * General fitting
             *
             * interpolating function through M nodes should have
             * greater RMS error than fitting it through the same M nodes
             */
            n = 100;
            ae_vector_set_length(&x2, n, _state);
            ae_vector_set_length(&y2, n, _state);
            ae_vector_set_length(&w2, n, _state);
            xmin = ae_maxrealnumber;
            xmax = -ae_maxrealnumber;
            for(i=0; i<=n-1; i++)
            {
                x2.ptr.p_double[i] = 2*ae_pi*ae_randomreal(_state);
                y2.ptr.p_double[i] = ae_sin(x2.ptr.p_double[i], _state);
                w2.ptr.p_double[i] = 1;
                xmin = ae_minreal(xmin, x2.ptr.p_double[i], _state);
                xmax = ae_maxreal(xmax, x2.ptr.p_double[i], _state);
            }
            ae_vector_set_length(&x, m, _state);
            ae_vector_set_length(&y, m, _state);
            for(i=0; i<=m-1; i++)
            {
                x.ptr.p_double[i] = xmin+(xmax-xmin)*i/(m-1);
                y.ptr.p_double[i] = ae_sin(x.ptr.p_double[i], _state);
            }
            barycentricbuildfloaterhormann(&x, &y, m, 3, &b1, _state);
            barycentricfitfloaterhormannwc(&x2, &y2, &w2, n, &xc, &yc, &dc, 0, m, &info, &b2, &rep, _state);
            if( info<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * calculate B1 (interpolant) RMS error, compare with B2 error
                 */
                v1 = 0;
                v2 = 0;
                for(i=0; i<=n-1; i++)
                {
                    v1 = v1+ae_sqr(barycentriccalc(&b1, x2.ptr.p_double[i], _state)-y2.ptr.p_double[i], _state);
                    v2 = v2+ae_sqr(barycentriccalc(&b2, x2.ptr.p_double[i], _state)-y2.ptr.p_double[i], _state);
                }
                v1 = ae_sqrt(v1/n, _state);
                v2 = ae_sqrt(v2/n, _state);
                *fiterrors = *fiterrors||ae_fp_greater(v2,v1);
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v2-rep.rmserror, _state),threshold);
            }
            
            /*
             * compare weighted and non-weighted
             */
            n = 20;
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                w.ptr.p_double[i] = 1;
            }
            barycentricfitfloaterhormannwc(&x, &y, &w, n, &xc, &yc, &dc, 0, m, &info, &b1, &rep, _state);
            barycentricfitfloaterhormann(&x, &y, n, m, &info2, &b2, &rep2, _state);
            if( info<=0||info2<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * calculate B1 (interpolant), compare with B2
                 * compare RMS errors
                 */
                t = 2*ae_randomreal(_state)-1;
                v1 = barycentriccalc(&b1, t, _state);
                v2 = barycentriccalc(&b2, t, _state);
                *fiterrors = *fiterrors||ae_fp_neq(v2,v1);
                *fiterrors = *fiterrors||ae_fp_neq(rep.rmserror,rep2.rmserror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.avgerror,rep2.avgerror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.avgrelerror,rep2.avgrelerror);
                *fiterrors = *fiterrors||ae_fp_neq(rep.maxerror,rep2.maxerror);
            }
        }
    }
    for(pass=1; pass<=passcount; pass++)
    {
        ae_assert(passcount>=2, "PassCount should be 2 or greater!", _state);
        
        /*
         * solve simple task (all X[] are the same, Y[] are specially
         * calculated to ensure simple form of all types of errors)
         * and check correctness of the errors calculated by subroutines
         *
         * First pass is done with zero Y[], other passes - with random Y[].
         * It should test both ability to correctly calculate errors and
         * ability to not fail while working with zeros :)
         */
        n = 4;
        if( pass==1 )
        {
            v1 = 0;
            v2 = 0;
            v = 0;
        }
        else
        {
            v1 = ae_randomreal(_state);
            v2 = ae_randomreal(_state);
            v = 1+ae_randomreal(_state);
        }
        ae_vector_set_length(&x, 4, _state);
        ae_vector_set_length(&y, 4, _state);
        ae_vector_set_length(&w, 4, _state);
        x.ptr.p_double[0] = 0;
        y.ptr.p_double[0] = v-v2;
        w.ptr.p_double[0] = 1;
        x.ptr.p_double[1] = 0;
        y.ptr.p_double[1] = v-v1;
        w.ptr.p_double[1] = 1;
        x.ptr.p_double[2] = 0;
        y.ptr.p_double[2] = v+v1;
        w.ptr.p_double[2] = 1;
        x.ptr.p_double[3] = 0;
        y.ptr.p_double[3] = v+v2;
        w.ptr.p_double[3] = 1;
        refrms = ae_sqrt((ae_sqr(v1, _state)+ae_sqr(v2, _state))/2, _state);
        refavg = (ae_fabs(v1, _state)+ae_fabs(v2, _state))/2;
        if( pass==1 )
        {
            refavgrel = 0;
        }
        else
        {
            refavgrel = 0.25*(ae_fabs(v2, _state)/ae_fabs(v-v2, _state)+ae_fabs(v1, _state)/ae_fabs(v-v1, _state)+ae_fabs(v1, _state)/ae_fabs(v+v1, _state)+ae_fabs(v2, _state)/ae_fabs(v+v2, _state));
        }
        refmax = ae_maxreal(v1, v2, _state);
        
        /*
         * Test errors correctness
         */
        barycentricfitfloaterhormann(&x, &y, 4, 2, &info, &b1, &rep, _state);
        if( info<=0 )
        {
            *fiterrors = ae_true;
        }
        else
        {
            s = barycentriccalc(&b1, 0, _state);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-v, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
    }
    ae_frame_leave(_state);
}


static void testlsfitunit_testsplinefitting(ae_bool* fiterrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    double threshold;
    double nonstrictthreshold;
    ae_int_t passcount;
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    ae_int_t k;
    ae_int_t pass;
    ae_vector x;
    ae_vector y;
    ae_vector w;
    ae_vector w2;
    ae_vector xc;
    ae_vector yc;
    ae_vector d;
    ae_vector dc;
    double sa;
    double sb;
    ae_int_t info;
    ae_int_t info1;
    ae_int_t info2;
    spline1dinterpolant c;
    spline1dinterpolant c2;
    spline1dfitreport rep;
    spline1dfitreport rep2;
    double s;
    double ds;
    double d2s;
    ae_int_t stype;
    double t;
    double v;
    double v1;
    double v2;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    double rho;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);
    _spline1dinterpolant_init(&c, _state, ae_true);
    _spline1dinterpolant_init(&c2, _state, ae_true);
    _spline1dfitreport_init(&rep, _state, ae_true);
    _spline1dfitreport_init(&rep2, _state, ae_true);

    
    /*
     * Valyes:
     * * pass count
     * * threshold - for tests which must be satisfied exactly
     * * nonstrictthreshold - for approximate tests
     */
    passcount = 20;
    threshold = 10000*ae_machineepsilon;
    nonstrictthreshold = 1.0E-4;
    *fiterrors = ae_false;
    
    /*
     * Test fitting by Cubic and Hermite splines (obsolete, but still supported)
     */
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Cubic splines
         * Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
         */
        for(m=4; m<=8; m++)
        {
            for(k=1; k<=4; k++)
            {
                if( k>=m )
                {
                    continue;
                }
                n = 100;
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&y, n, _state);
                ae_vector_set_length(&w, n, _state);
                ae_vector_set_length(&xc, 4, _state);
                ae_vector_set_length(&yc, 4, _state);
                ae_vector_set_length(&dc, 4, _state);
                sa = 1+ae_randomreal(_state);
                sb = 2*ae_randomreal(_state)-1;
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                    y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                }
                xc.ptr.p_double[0] = sb;
                yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[0] = 0;
                xc.ptr.p_double[1] = sb;
                yc.ptr.p_double[1] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[1] = 1;
                xc.ptr.p_double[2] = sa+sb;
                yc.ptr.p_double[2] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[2] = 0;
                xc.ptr.p_double[3] = sa+sb;
                yc.ptr.p_double[3] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[3] = 1;
                spline1dfitcubicwc(&x, &y, &w, n, &xc, &yc, &dc, k, m, &info, &c, &rep, _state);
                if( info<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    
                    /*
                     * Check that constraints are satisfied
                     */
                    for(i=0; i<=k-1; i++)
                    {
                        spline1ddiff(&c, xc.ptr.p_double[i], &s, &ds, &d2s, _state);
                        if( dc.ptr.p_int[i]==0 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-yc.ptr.p_double[i], _state),threshold);
                        }
                        if( dc.ptr.p_int[i]==1 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(ds-yc.ptr.p_double[i], _state),threshold);
                        }
                        if( dc.ptr.p_int[i]==2 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(d2s-yc.ptr.p_double[i], _state),threshold);
                        }
                    }
                }
            }
        }
        
        /*
         * Cubic splines
         * Ability to handle one internal constraint
         */
        for(m=4; m<=8; m++)
        {
            n = 100;
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            ae_vector_set_length(&xc, 1, _state);
            ae_vector_set_length(&yc, 1, _state);
            ae_vector_set_length(&dc, 1, _state);
            sa = 1+ae_randomreal(_state);
            sb = 2*ae_randomreal(_state)-1;
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                w.ptr.p_double[i] = 1+ae_randomreal(_state);
            }
            xc.ptr.p_double[0] = sa*ae_randomreal(_state)+sb;
            yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
            dc.ptr.p_int[0] = ae_randominteger(2, _state);
            spline1dfitcubicwc(&x, &y, &w, n, &xc, &yc, &dc, 1, m, &info, &c, &rep, _state);
            if( info<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * Check that constraints are satisfied
                 */
                spline1ddiff(&c, xc.ptr.p_double[0], &s, &ds, &d2s, _state);
                if( dc.ptr.p_int[0]==0 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-yc.ptr.p_double[0], _state),threshold);
                }
                if( dc.ptr.p_int[0]==1 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(ds-yc.ptr.p_double[0], _state),threshold);
                }
                if( dc.ptr.p_int[0]==2 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(d2s-yc.ptr.p_double[0], _state),threshold);
                }
            }
        }
        
        /*
         * Hermite splines
         * Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
         */
        for(m=4; m<=8; m++)
        {
            for(k=1; k<=4; k++)
            {
                if( k>=m )
                {
                    continue;
                }
                if( m%2!=0 )
                {
                    continue;
                }
                n = 100;
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&y, n, _state);
                ae_vector_set_length(&w, n, _state);
                ae_vector_set_length(&xc, 4, _state);
                ae_vector_set_length(&yc, 4, _state);
                ae_vector_set_length(&dc, 4, _state);
                sa = 1+ae_randomreal(_state);
                sb = 2*ae_randomreal(_state)-1;
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                    y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                }
                xc.ptr.p_double[0] = sb;
                yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[0] = 0;
                xc.ptr.p_double[1] = sb;
                yc.ptr.p_double[1] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[1] = 1;
                xc.ptr.p_double[2] = sa+sb;
                yc.ptr.p_double[2] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[2] = 0;
                xc.ptr.p_double[3] = sa+sb;
                yc.ptr.p_double[3] = 2*ae_randomreal(_state)-1;
                dc.ptr.p_int[3] = 1;
                spline1dfithermitewc(&x, &y, &w, n, &xc, &yc, &dc, k, m, &info, &c, &rep, _state);
                if( info<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    
                    /*
                     * Check that constraints are satisfied
                     */
                    for(i=0; i<=k-1; i++)
                    {
                        spline1ddiff(&c, xc.ptr.p_double[i], &s, &ds, &d2s, _state);
                        if( dc.ptr.p_int[i]==0 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-yc.ptr.p_double[i], _state),threshold);
                        }
                        if( dc.ptr.p_int[i]==1 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(ds-yc.ptr.p_double[i], _state),threshold);
                        }
                        if( dc.ptr.p_int[i]==2 )
                        {
                            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(d2s-yc.ptr.p_double[i], _state),threshold);
                        }
                    }
                }
            }
        }
        
        /*
         * Hermite splines
         * Ability to handle one internal constraint
         */
        for(m=4; m<=8; m++)
        {
            if( m%2!=0 )
            {
                continue;
            }
            n = 100;
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            ae_vector_set_length(&xc, 1, _state);
            ae_vector_set_length(&yc, 1, _state);
            ae_vector_set_length(&dc, 1, _state);
            sa = 1+ae_randomreal(_state);
            sb = 2*ae_randomreal(_state)-1;
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                w.ptr.p_double[i] = 1+ae_randomreal(_state);
            }
            xc.ptr.p_double[0] = sa*ae_randomreal(_state)+sb;
            yc.ptr.p_double[0] = 2*ae_randomreal(_state)-1;
            dc.ptr.p_int[0] = ae_randominteger(2, _state);
            spline1dfithermitewc(&x, &y, &w, n, &xc, &yc, &dc, 1, m, &info, &c, &rep, _state);
            if( info<=0 )
            {
                *fiterrors = ae_true;
            }
            else
            {
                
                /*
                 * Check that constraints are satisfied
                 */
                spline1ddiff(&c, xc.ptr.p_double[0], &s, &ds, &d2s, _state);
                if( dc.ptr.p_int[0]==0 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-yc.ptr.p_double[0], _state),threshold);
                }
                if( dc.ptr.p_int[0]==1 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(ds-yc.ptr.p_double[0], _state),threshold);
                }
                if( dc.ptr.p_int[0]==2 )
                {
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(d2s-yc.ptr.p_double[0], _state),threshold);
                }
            }
        }
    }
    for(m=4; m<=8; m++)
    {
        for(stype=0; stype<=1; stype++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                if( stype==1&&m%2!=0 )
                {
                    continue;
                }
                
                /*
                 * cubic/Hermite spline fitting:
                 * * generate "template spline" C2
                 * * generate 2*N points from C2, such that result of
                 *   ideal fit should be equal to C2
                 * * fit, store in C
                 * * compare C and C2
                 */
                sa = 1+ae_randomreal(_state);
                sb = 2*ae_randomreal(_state)-1;
                if( stype==0 )
                {
                    ae_vector_set_length(&x, m-2, _state);
                    ae_vector_set_length(&y, m-2, _state);
                    for(i=0; i<=m-2-1; i++)
                    {
                        x.ptr.p_double[i] = sa*i/(m-2-1)+sb;
                        y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    spline1dbuildcubic(&x, &y, m-2, 1, 2*ae_randomreal(_state)-1, 1, 2*ae_randomreal(_state)-1, &c2, _state);
                }
                if( stype==1 )
                {
                    ae_vector_set_length(&x, m/2, _state);
                    ae_vector_set_length(&y, m/2, _state);
                    ae_vector_set_length(&d, m/2, _state);
                    for(i=0; i<=m/2-1; i++)
                    {
                        x.ptr.p_double[i] = sa*i/(m/2-1)+sb;
                        y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                        d.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    spline1dbuildhermite(&x, &y, &d, m/2, &c2, _state);
                }
                n = 50;
                ae_vector_set_length(&x, 2*n, _state);
                ae_vector_set_length(&y, 2*n, _state);
                ae_vector_set_length(&w, 2*n, _state);
                for(i=0; i<=n-1; i++)
                {
                    
                    /*
                     * "if i=0" and "if i=1" are needed to
                     * synchronize interval size for C2 and
                     * spline being fitted (i.e. C).
                     */
                    t = ae_randomreal(_state);
                    x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                    if( i==0 )
                    {
                        x.ptr.p_double[i] = sb;
                    }
                    if( i==1 )
                    {
                        x.ptr.p_double[i] = sa+sb;
                    }
                    v = spline1dcalc(&c2, x.ptr.p_double[i], _state);
                    y.ptr.p_double[i] = v+t;
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                    x.ptr.p_double[n+i] = x.ptr.p_double[i];
                    y.ptr.p_double[n+i] = v-t;
                    w.ptr.p_double[n+i] = w.ptr.p_double[i];
                }
                if( stype==0 )
                {
                    spline1dfitcubicwc(&x, &y, &w, 2*n, &xc, &yc, &dc, 0, m, &info, &c, &rep, _state);
                }
                if( stype==1 )
                {
                    spline1dfithermitewc(&x, &y, &w, 2*n, &xc, &yc, &dc, 0, m, &info, &c, &rep, _state);
                }
                if( info<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    for(i=0; i<=n-1; i++)
                    {
                        v = sa*ae_randomreal(_state)+sb;
                        *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(spline1dcalc(&c, v, _state)-spline1dcalc(&c2, v, _state), _state),threshold);
                    }
                }
            }
        }
    }
    for(m=4; m<=8; m++)
    {
        for(pass=1; pass<=passcount; pass++)
        {
            
            /*
             * prepare points/weights
             */
            sa = 1+ae_randomreal(_state);
            sb = 2*ae_randomreal(_state)-1;
            n = 10+ae_randominteger(10, _state);
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
            ae_vector_set_length(&w, n, _state);
            for(i=0; i<=n-1; i++)
            {
                x.ptr.p_double[i] = sa*ae_randomreal(_state)+sb;
                y.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                w.ptr.p_double[i] = 1;
            }
            
            /*
             * Fit cubic with unity weights, without weights, then compare
             */
            if( m>=4 )
            {
                spline1dfitcubicwc(&x, &y, &w, n, &xc, &yc, &dc, 0, m, &info1, &c, &rep, _state);
                spline1dfitcubic(&x, &y, n, m, &info2, &c2, &rep2, _state);
                if( info1<=0||info2<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    for(i=0; i<=n-1; i++)
                    {
                        v = sa*ae_randomreal(_state)+sb;
                        *fiterrors = *fiterrors||ae_fp_neq(spline1dcalc(&c, v, _state),spline1dcalc(&c2, v, _state));
                        *fiterrors = *fiterrors||ae_fp_neq(rep.taskrcond,rep2.taskrcond);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.rmserror,rep2.rmserror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.avgerror,rep2.avgerror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.avgrelerror,rep2.avgrelerror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.maxerror,rep2.maxerror);
                    }
                }
            }
            
            /*
             * Fit Hermite with unity weights, without weights, then compare
             */
            if( m>=4&&m%2==0 )
            {
                spline1dfithermitewc(&x, &y, &w, n, &xc, &yc, &dc, 0, m, &info1, &c, &rep, _state);
                spline1dfithermite(&x, &y, n, m, &info2, &c2, &rep2, _state);
                if( info1<=0||info2<=0 )
                {
                    *fiterrors = ae_true;
                }
                else
                {
                    for(i=0; i<=n-1; i++)
                    {
                        v = sa*ae_randomreal(_state)+sb;
                        *fiterrors = *fiterrors||ae_fp_neq(spline1dcalc(&c, v, _state),spline1dcalc(&c2, v, _state));
                        *fiterrors = *fiterrors||ae_fp_neq(rep.taskrcond,rep2.taskrcond);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.rmserror,rep2.rmserror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.avgerror,rep2.avgerror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.avgrelerror,rep2.avgrelerror);
                        *fiterrors = *fiterrors||ae_fp_neq(rep.maxerror,rep2.maxerror);
                    }
                }
            }
        }
    }
    
    /*
     * check basic properties of penalized splines which are
     * preserved independently of Rho parameter.
     */
    for(m=4; m<=10; m++)
    {
        for(k=-5; k<=5; k++)
        {
            rho = k;
            
            /*
             * when we have two points (even with different weights),
             * resulting spline must be equal to the straight line
             */
            ae_vector_set_length(&x, 2, _state);
            ae_vector_set_length(&y, 2, _state);
            ae_vector_set_length(&w, 2, _state);
            x.ptr.p_double[0] = -0.5-ae_randomreal(_state);
            y.ptr.p_double[0] = 0.5+ae_randomreal(_state);
            w.ptr.p_double[0] = 1+ae_randomreal(_state);
            x.ptr.p_double[1] = 0.5+ae_randomreal(_state);
            y.ptr.p_double[1] = 0.5+ae_randomreal(_state);
            w.ptr.p_double[1] = 1+ae_randomreal(_state);
            spline1dfitpenalized(&x, &y, 2, m, rho, &info, &c, &rep, _state);
            if( info>0 )
            {
                v = 2*ae_randomreal(_state)-1;
                v1 = (v-x.ptr.p_double[0])/(x.ptr.p_double[1]-x.ptr.p_double[0])*y.ptr.p_double[1]+(v-x.ptr.p_double[1])/(x.ptr.p_double[0]-x.ptr.p_double[1])*y.ptr.p_double[0];
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v1-spline1dcalc(&c, v, _state), _state),nonstrictthreshold);
            }
            else
            {
                *fiterrors = ae_true;
            }
            spline1dfitpenalizedw(&x, &y, &w, 2, m, rho, &info, &c, &rep, _state);
            if( info>0 )
            {
                v = 2*ae_randomreal(_state)-1;
                v1 = (v-x.ptr.p_double[0])/(x.ptr.p_double[1]-x.ptr.p_double[0])*y.ptr.p_double[1]+(v-x.ptr.p_double[1])/(x.ptr.p_double[0]-x.ptr.p_double[1])*y.ptr.p_double[0];
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v1-spline1dcalc(&c, v, _state), _state),nonstrictthreshold);
            }
            else
            {
                *fiterrors = ae_true;
            }
            
            /*
             * spline fitting is invariant with respect to
             * scaling of weights (of course, ANY fitting algorithm
             * must be invariant, but we want to test this property
             * just to be sure that it is correctly implemented)
             */
            for(n=2; n<=2*m; n++)
            {
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&y, n, _state);
                ae_vector_set_length(&w, n, _state);
                ae_vector_set_length(&w2, n, _state);
                s = 1+ae_exp(10*ae_randomreal(_state), _state);
                for(i=0; i<=n-1; i++)
                {
                    x.ptr.p_double[i] = (double)i/(double)(n-1);
                    y.ptr.p_double[i] = ae_randomreal(_state);
                    w.ptr.p_double[i] = 0.1+ae_randomreal(_state);
                    w2.ptr.p_double[i] = w.ptr.p_double[i]*s;
                }
                spline1dfitpenalizedw(&x, &y, &w, n, m, rho, &info, &c, &rep, _state);
                spline1dfitpenalizedw(&x, &y, &w2, n, m, rho, &info2, &c2, &rep2, _state);
                if( info>0&&info2>0 )
                {
                    v = ae_randomreal(_state);
                    v1 = spline1dcalc(&c, v, _state);
                    v2 = spline1dcalc(&c2, v, _state);
                    *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(v1-v2, _state),nonstrictthreshold);
                }
                else
                {
                    *fiterrors = ae_true;
                }
            }
        }
    }
    
    /*
     * Advanced proprties:
     * * penalized spline with M about 5*N and sufficiently small Rho
     *   must pass through all points on equidistant grid
     */
    for(n=2; n<=10; n++)
    {
        m = 5*n;
        rho = -5;
        ae_vector_set_length(&x, n, _state);
        ae_vector_set_length(&y, n, _state);
        ae_vector_set_length(&w, n, _state);
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = (double)i/(double)(n-1);
            y.ptr.p_double[i] = ae_randomreal(_state);
            w.ptr.p_double[i] = 0.1+ae_randomreal(_state);
        }
        spline1dfitpenalized(&x, &y, n, m, rho, &info, &c, &rep, _state);
        if( info>0 )
        {
            for(i=0; i<=n-1; i++)
            {
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(y.ptr.p_double[i]-spline1dcalc(&c, x.ptr.p_double[i], _state), _state),nonstrictthreshold);
            }
        }
        else
        {
            *fiterrors = ae_true;
        }
        spline1dfitpenalizedw(&x, &y, &w, n, m, rho, &info, &c, &rep, _state);
        if( info>0 )
        {
            for(i=0; i<=n-1; i++)
            {
                *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(y.ptr.p_double[i]-spline1dcalc(&c, x.ptr.p_double[i], _state), _state),nonstrictthreshold);
            }
        }
        else
        {
            *fiterrors = ae_true;
        }
    }
    
    /*
     * Check correctness of error reports
     */
    for(pass=1; pass<=passcount; pass++)
    {
        ae_assert(passcount>=2, "PassCount should be 2 or greater!", _state);
        
        /*
         * solve simple task (all X[] are the same, Y[] are specially
         * calculated to ensure simple form of all types of errors)
         * and check correctness of the errors calculated by subroutines
         *
         * First pass is done with zero Y[], other passes - with random Y[].
         * It should test both ability to correctly calculate errors and
         * ability to not fail while working with zeros :)
         */
        n = 4;
        if( pass==1 )
        {
            v1 = 0;
            v2 = 0;
            v = 0;
        }
        else
        {
            v1 = ae_randomreal(_state);
            v2 = ae_randomreal(_state);
            v = 1+ae_randomreal(_state);
        }
        ae_vector_set_length(&x, 4, _state);
        ae_vector_set_length(&y, 4, _state);
        ae_vector_set_length(&w, 4, _state);
        x.ptr.p_double[0] = 0;
        y.ptr.p_double[0] = v-v2;
        w.ptr.p_double[0] = 1;
        x.ptr.p_double[1] = 0;
        y.ptr.p_double[1] = v-v1;
        w.ptr.p_double[1] = 1;
        x.ptr.p_double[2] = 0;
        y.ptr.p_double[2] = v+v1;
        w.ptr.p_double[2] = 1;
        x.ptr.p_double[3] = 0;
        y.ptr.p_double[3] = v+v2;
        w.ptr.p_double[3] = 1;
        refrms = ae_sqrt((ae_sqr(v1, _state)+ae_sqr(v2, _state))/2, _state);
        refavg = (ae_fabs(v1, _state)+ae_fabs(v2, _state))/2;
        if( pass==1 )
        {
            refavgrel = 0;
        }
        else
        {
            refavgrel = 0.25*(ae_fabs(v2, _state)/ae_fabs(v-v2, _state)+ae_fabs(v1, _state)/ae_fabs(v-v1, _state)+ae_fabs(v1, _state)/ae_fabs(v+v1, _state)+ae_fabs(v2, _state)/ae_fabs(v+v2, _state));
        }
        refmax = ae_maxreal(v1, v2, _state);
        
        /*
         * Test penalized spline
         */
        spline1dfitpenalizedw(&x, &y, &w, 4, 4, 0.0, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *fiterrors = ae_true;
        }
        else
        {
            s = spline1dcalc(&c, 0, _state);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-v, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
        
        /*
         * Test cubic fitting
         */
        spline1dfitcubic(&x, &y, 4, 4, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *fiterrors = ae_true;
        }
        else
        {
            s = spline1dcalc(&c, 0, _state);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-v, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
        
        /*
         * Test Hermite fitting
         */
        spline1dfithermite(&x, &y, 4, 4, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *fiterrors = ae_true;
        }
        else
        {
            s = spline1dcalc(&c, 0, _state);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(s-v, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *fiterrors = *fiterrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
    }
    ae_frame_leave(_state);
}


static void testlsfitunit_testgeneralfitting(ae_bool* llserrors,
     ae_bool* nlserrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    double threshold;
    double nlthreshold;
    ae_int_t maxn;
    ae_int_t maxm;
    ae_int_t passcount;
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t pass;
    double xscale;
    double diffstep;
    ae_vector x;
    ae_vector y;
    ae_vector w;
    ae_vector w2;
    ae_vector c;
    ae_vector c2;
    ae_matrix a;
    ae_matrix a2;
    ae_matrix cm;
    double v;
    double v1;
    double v2;
    lsfitreport rep;
    lsfitreport rep2;
    ae_int_t info;
    ae_int_t info2;
    double refrms;
    double refavg;
    double refavgrel;
    double refmax;
    lsfitstate state;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c2, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a2, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cm, 0, 0, DT_REAL, _state, ae_true);
    _lsfitreport_init(&rep, _state, ae_true);
    _lsfitreport_init(&rep2, _state, ae_true);
    _lsfitstate_init(&state, _state, ae_true);

    *llserrors = ae_false;
    *nlserrors = ae_false;
    threshold = 10000*ae_machineepsilon;
    nlthreshold = 0.00001;
    diffstep = 0.0001;
    maxn = 6;
    maxm = 6;
    passcount = 4;
    
    /*
     * Testing unconstrained least squares (linear/nonlinear)
     */
    for(n=1; n<=maxn; n++)
    {
        for(m=1; m<=maxm; m++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                /*
                 * Solve non-degenerate linear least squares task
                 * Use Chebyshev basis. Its condition number is very good.
                 */
                ae_matrix_set_length(&a, n, m, _state);
                ae_vector_set_length(&x, n, _state);
                ae_vector_set_length(&y, n, _state);
                ae_vector_set_length(&w, n, _state);
                xscale = 0.9+0.1*ae_randomreal(_state);
                for(i=0; i<=n-1; i++)
                {
                    if( n==1 )
                    {
                        x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    else
                    {
                        x.ptr.p_double[i] = xscale*((double)(2*i)/(double)(n-1)-1);
                    }
                    y.ptr.p_double[i] = 3*x.ptr.p_double[i]+ae_exp(x.ptr.p_double[i], _state);
                    w.ptr.p_double[i] = 1+ae_randomreal(_state);
                    a.ptr.pp_double[i][0] = 1;
                    if( m>1 )
                    {
                        a.ptr.pp_double[i][1] = x.ptr.p_double[i];
                    }
                    for(j=2; j<=m-1; j++)
                    {
                        a.ptr.pp_double[i][j] = 2*x.ptr.p_double[i]*a.ptr.pp_double[i][j-1]-a.ptr.pp_double[i][j-2];
                    }
                }
                
                /*
                 * 1. test weighted fitting (optimality)
                 * 2. Solve degenerate least squares task built on the basis
                 *    of previous task
                 */
                lsfitlinearw(&y, &w, &a, n, m, &info, &c, &rep, _state);
                if( info<=0 )
                {
                    *llserrors = ae_true;
                }
                else
                {
                    *llserrors = *llserrors||!testlsfitunit_isglssolution(n, m, 0, &y, &w, &a, &cm, &c, _state);
                }
                ae_matrix_set_length(&a2, n, 2*m, _state);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        a2.ptr.pp_double[i][2*j+0] = a.ptr.pp_double[i][j];
                        a2.ptr.pp_double[i][2*j+1] = a.ptr.pp_double[i][j];
                    }
                }
                lsfitlinearw(&y, &w, &a2, n, 2*m, &info, &c2, &rep, _state);
                if( info<=0 )
                {
                    *llserrors = ae_true;
                }
                else
                {
                    
                    /*
                     * test answer correctness using design matrix properties
                     * and previous task solution
                     */
                    for(j=0; j<=m-1; j++)
                    {
                        *llserrors = *llserrors||ae_fp_greater(ae_fabs(c2.ptr.p_double[2*j+0]+c2.ptr.p_double[2*j+1]-c.ptr.p_double[j], _state),threshold);
                    }
                }
                
                /*
                 * test non-weighted fitting
                 */
                ae_vector_set_length(&w2, n, _state);
                for(i=0; i<=n-1; i++)
                {
                    w2.ptr.p_double[i] = 1;
                }
                lsfitlinearw(&y, &w2, &a, n, m, &info, &c, &rep, _state);
                lsfitlinear(&y, &a, n, m, &info2, &c2, &rep2, _state);
                if( info<=0||info2<=0 )
                {
                    *llserrors = ae_true;
                }
                else
                {
                    
                    /*
                     * test answer correctness
                     */
                    for(j=0; j<=m-1; j++)
                    {
                        *llserrors = *llserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[j]-c2.ptr.p_double[j], _state),threshold);
                    }
                    *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.taskrcond-rep2.taskrcond, _state),threshold);
                }
                
                /*
                 * test nonlinear fitting on the linear task
                 * (only non-degenerate tasks are tested)
                 * and compare with answer from linear fitting subroutine
                 */
                if( n>=m )
                {
                    ae_vector_set_length(&c2, m, _state);
                    
                    /*
                     * test function/gradient/Hessian-based weighted fitting
                     */
                    lsfitlinearw(&y, &w, &a, n, m, &info, &c, &rep, _state);
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatewf(&a, &y, &w, &c2, n, m, m, diffstep, &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 0, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatewfg(&a, &y, &w, &c2, n, m, m, ae_fp_greater(ae_randomreal(_state),0.5), &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 1, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatewfgh(&a, &y, &w, &c2, n, m, m, &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 2, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                    
                    /*
                     * test gradient-only or Hessian-based fitting without weights
                     */
                    lsfitlinear(&y, &a, n, m, &info, &c, &rep, _state);
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatef(&a, &y, &c2, n, m, m, diffstep, &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 0, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatefg(&a, &y, &c2, n, m, m, ae_fp_greater(ae_randomreal(_state),0.5), &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 1, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        c2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                    }
                    lsfitcreatefgh(&a, &y, &c2, n, m, m, &state, _state);
                    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
                    testlsfitunit_fitlinearnonlinear(m, 2, &a, &state, nlserrors, _state);
                    lsfitresults(&state, &info, &c2, &rep2, _state);
                    if( info<=0 )
                    {
                        *nlserrors = ae_true;
                    }
                    else
                    {
                        for(i=0; i<=m-1; i++)
                        {
                            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[i]-c2.ptr.p_double[i], _state),100*nlthreshold);
                        }
                    }
                }
            }
        }
        
        /*
         * test correctness of the RCond field
         */
        ae_matrix_set_length(&a, n-1+1, n-1+1, _state);
        ae_vector_set_length(&x, n-1+1, _state);
        ae_vector_set_length(&y, n-1+1, _state);
        ae_vector_set_length(&w, n-1+1, _state);
        v1 = ae_maxrealnumber;
        v2 = ae_minrealnumber;
        for(i=0; i<=n-1; i++)
        {
            x.ptr.p_double[i] = 0.1+0.9*ae_randomreal(_state);
            y.ptr.p_double[i] = 0.1+0.9*ae_randomreal(_state);
            w.ptr.p_double[i] = 1;
            for(j=0; j<=n-1; j++)
            {
                if( i==j )
                {
                    a.ptr.pp_double[i][i] = 0.1+0.9*ae_randomreal(_state);
                    v1 = ae_minreal(v1, a.ptr.pp_double[i][i], _state);
                    v2 = ae_maxreal(v2, a.ptr.pp_double[i][i], _state);
                }
                else
                {
                    a.ptr.pp_double[i][j] = 0;
                }
            }
        }
        lsfitlinearw(&y, &w, &a, n, n, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *llserrors = ae_true;
        }
        else
        {
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.taskrcond-v1/v2, _state),threshold);
        }
    }
    
    /*
     * Test constrained least squares
     */
    for(pass=1; pass<=passcount; pass++)
    {
        for(n=1; n<=maxn; n++)
        {
            for(m=1; m<=maxm; m++)
            {
                
                /*
                 * test for K<>0
                 */
                for(k=1; k<=m-1; k++)
                {
                    
                    /*
                     * Prepare Chebyshev basis. Its condition number is very good.
                     * Prepare constraints (random numbers)
                     */
                    ae_matrix_set_length(&a, n, m, _state);
                    ae_vector_set_length(&x, n, _state);
                    ae_vector_set_length(&y, n, _state);
                    ae_vector_set_length(&w, n, _state);
                    xscale = 0.9+0.1*ae_randomreal(_state);
                    for(i=0; i<=n-1; i++)
                    {
                        if( n==1 )
                        {
                            x.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                        }
                        else
                        {
                            x.ptr.p_double[i] = xscale*((double)(2*i)/(double)(n-1)-1);
                        }
                        y.ptr.p_double[i] = 3*x.ptr.p_double[i]+ae_exp(x.ptr.p_double[i], _state);
                        w.ptr.p_double[i] = 1+ae_randomreal(_state);
                        a.ptr.pp_double[i][0] = 1;
                        if( m>1 )
                        {
                            a.ptr.pp_double[i][1] = x.ptr.p_double[i];
                        }
                        for(j=2; j<=m-1; j++)
                        {
                            a.ptr.pp_double[i][j] = 2*x.ptr.p_double[i]*a.ptr.pp_double[i][j-1]-a.ptr.pp_double[i][j-2];
                        }
                    }
                    ae_matrix_set_length(&cm, k, m+1, _state);
                    for(i=0; i<=k-1; i++)
                    {
                        for(j=0; j<=m; j++)
                        {
                            cm.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                        }
                    }
                    
                    /*
                     * Solve constrained task
                     */
                    lsfitlinearwc(&y, &w, &a, &cm, n, m, k, &info, &c, &rep, _state);
                    if( info<=0 )
                    {
                        *llserrors = ae_true;
                    }
                    else
                    {
                        *llserrors = *llserrors||!testlsfitunit_isglssolution(n, m, k, &y, &w, &a, &cm, &c, _state);
                    }
                    
                    /*
                     * test non-weighted fitting
                     */
                    ae_vector_set_length(&w2, n, _state);
                    for(i=0; i<=n-1; i++)
                    {
                        w2.ptr.p_double[i] = 1;
                    }
                    lsfitlinearwc(&y, &w2, &a, &cm, n, m, k, &info, &c, &rep, _state);
                    lsfitlinearc(&y, &a, &cm, n, m, k, &info2, &c2, &rep2, _state);
                    if( info<=0||info2<=0 )
                    {
                        *llserrors = ae_true;
                    }
                    else
                    {
                        
                        /*
                         * test answer correctness
                         */
                        for(j=0; j<=m-1; j++)
                        {
                            *llserrors = *llserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[j]-c2.ptr.p_double[j], _state),threshold);
                        }
                        *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.taskrcond-rep2.taskrcond, _state),threshold);
                    }
                }
            }
        }
    }
    
    /*
     * nonlinear task for nonlinear fitting:
     *
     *     f(X,C) = 1/(1+C*X^2),
     *     C(true) = 2.
     */
    n = 100;
    ae_vector_set_length(&c, 1, _state);
    c.ptr.p_double[0] = 1+2*ae_randomreal(_state);
    ae_matrix_set_length(&a, n, 1, _state);
    ae_vector_set_length(&y, n, _state);
    for(i=0; i<=n-1; i++)
    {
        a.ptr.pp_double[i][0] = 4*ae_randomreal(_state)-2;
        y.ptr.p_double[i] = 1/(1+2*ae_sqr(a.ptr.pp_double[i][0], _state));
    }
    lsfitcreatefg(&a, &y, &c, n, 1, 1, ae_true, &state, _state);
    lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
    while(lsfititeration(&state, _state))
    {
        if( state.needf )
        {
            state.f = 1/(1+state.c.ptr.p_double[0]*ae_sqr(state.x.ptr.p_double[0], _state));
        }
        if( state.needfg )
        {
            state.f = 1/(1+state.c.ptr.p_double[0]*ae_sqr(state.x.ptr.p_double[0], _state));
            state.g.ptr.p_double[0] = -ae_sqr(state.x.ptr.p_double[0], _state)/ae_sqr(1+state.c.ptr.p_double[0]*ae_sqr(state.x.ptr.p_double[0], _state), _state);
        }
    }
    lsfitresults(&state, &info, &c, &rep, _state);
    if( info<=0 )
    {
        *nlserrors = ae_true;
    }
    else
    {
        *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[0]-2, _state),100*nlthreshold);
    }
    
    /*
     * solve simple task (fitting by constant function) and check
     * correctness of the errors calculated by subroutines
     */
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * test on task with non-zero Yi
         */
        n = 4;
        v1 = ae_randomreal(_state);
        v2 = ae_randomreal(_state);
        v = 1+ae_randomreal(_state);
        ae_vector_set_length(&c, 1, _state);
        c.ptr.p_double[0] = 1+2*ae_randomreal(_state);
        ae_matrix_set_length(&a, 4, 1, _state);
        ae_vector_set_length(&y, 4, _state);
        a.ptr.pp_double[0][0] = 1;
        y.ptr.p_double[0] = v-v2;
        a.ptr.pp_double[1][0] = 1;
        y.ptr.p_double[1] = v-v1;
        a.ptr.pp_double[2][0] = 1;
        y.ptr.p_double[2] = v+v1;
        a.ptr.pp_double[3][0] = 1;
        y.ptr.p_double[3] = v+v2;
        refrms = ae_sqrt((ae_sqr(v1, _state)+ae_sqr(v2, _state))/2, _state);
        refavg = (ae_fabs(v1, _state)+ae_fabs(v2, _state))/2;
        refavgrel = 0.25*(ae_fabs(v2, _state)/ae_fabs(v-v2, _state)+ae_fabs(v1, _state)/ae_fabs(v-v1, _state)+ae_fabs(v1, _state)/ae_fabs(v+v1, _state)+ae_fabs(v2, _state)/ae_fabs(v+v2, _state));
        refmax = ae_maxreal(v1, v2, _state);
        
        /*
         * Test LLS
         */
        lsfitlinear(&y, &a, 4, 1, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *llserrors = ae_true;
        }
        else
        {
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[0]-v, _state),threshold);
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *llserrors = *llserrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
        
        /*
         * Test NLS
         */
        lsfitcreatefg(&a, &y, &c, 4, 1, 1, ae_true, &state, _state);
        lsfitsetcond(&state, 0.0, nlthreshold, 0, _state);
        while(lsfititeration(&state, _state))
        {
            if( state.needf )
            {
                state.f = state.c.ptr.p_double[0];
            }
            if( state.needfg )
            {
                state.f = state.c.ptr.p_double[0];
                state.g.ptr.p_double[0] = 1;
            }
        }
        lsfitresults(&state, &info, &c, &rep, _state);
        if( info<=0 )
        {
            *nlserrors = ae_true;
        }
        else
        {
            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(c.ptr.p_double[0]-v, _state),threshold);
            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(rep.rmserror-refrms, _state),threshold);
            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(rep.avgerror-refavg, _state),threshold);
            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(rep.avgrelerror-refavgrel, _state),threshold);
            *nlserrors = *nlserrors||ae_fp_greater(ae_fabs(rep.maxerror-refmax, _state),threshold);
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Tests whether C is solution of (possibly) constrained LLS problem
*************************************************************************/
static ae_bool testlsfitunit_isglssolution(ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     /* Real    */ ae_vector* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _c;
    ae_int_t i;
    ae_int_t j;
    ae_vector c2;
    ae_vector sv;
    ae_vector deltac;
    ae_vector deltaproj;
    ae_matrix u;
    ae_matrix vt;
    double v;
    double s1;
    double s2;
    double s3;
    double delta;
    double threshold;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_c, c, _state, ae_true);
    c = &_c;
    ae_vector_init(&c2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sv, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&deltac, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&deltaproj, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * Setup.
     * Threshold is small because CMatrix may be ill-conditioned
     */
    delta = 0.001;
    threshold = ae_sqrt(ae_machineepsilon, _state);
    ae_vector_set_length(&c2, m, _state);
    ae_vector_set_length(&deltac, m, _state);
    ae_vector_set_length(&deltaproj, m, _state);
    
    /*
     * test whether C is feasible point or not (projC must be close to C)
     */
    for(i=0; i<=k-1; i++)
    {
        v = ae_v_dotproduct(&cmatrix->ptr.pp_double[i][0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        if( ae_fp_greater(ae_fabs(v-cmatrix->ptr.pp_double[i][m], _state),threshold) )
        {
            result = ae_false;
            ae_frame_leave(_state);
            return result;
        }
    }
    
    /*
     * find orthogonal basis of Null(CMatrix) (stored in rows from K to M-1)
     */
    if( k>0 )
    {
        rmatrixsvd(cmatrix, k, m, 0, 2, 2, &sv, &u, &vt, _state);
    }
    
    /*
     * Test result
     */
    result = ae_true;
    s1 = testlsfitunit_getglserror(n, m, y, w, fmatrix, c, _state);
    for(j=0; j<=m-1; j++)
    {
        
        /*
         * prepare modification of C which leave us in the feasible set.
         *
         * let deltaC be increment on Jth coordinate, then project
         * deltaC in the Null(CMatrix) and store result in DeltaProj
         */
        ae_v_move(&c2.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        for(i=0; i<=m-1; i++)
        {
            if( i==j )
            {
                deltac.ptr.p_double[i] = delta;
            }
            else
            {
                deltac.ptr.p_double[i] = 0;
            }
        }
        if( k==0 )
        {
            ae_v_move(&deltaproj.ptr.p_double[0], 1, &deltac.ptr.p_double[0], 1, ae_v_len(0,m-1));
        }
        else
        {
            for(i=0; i<=m-1; i++)
            {
                deltaproj.ptr.p_double[i] = 0;
            }
            for(i=k; i<=m-1; i++)
            {
                v = ae_v_dotproduct(&vt.ptr.pp_double[i][0], 1, &deltac.ptr.p_double[0], 1, ae_v_len(0,m-1));
                ae_v_addd(&deltaproj.ptr.p_double[0], 1, &vt.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
            }
        }
        
        /*
         * now we have DeltaProj such that if C is feasible,
         * then C+DeltaProj is feasible too
         */
        ae_v_move(&c2.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        ae_v_add(&c2.ptr.p_double[0], 1, &deltaproj.ptr.p_double[0], 1, ae_v_len(0,m-1));
        s2 = testlsfitunit_getglserror(n, m, y, w, fmatrix, &c2, _state);
        ae_v_move(&c2.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        ae_v_sub(&c2.ptr.p_double[0], 1, &deltaproj.ptr.p_double[0], 1, ae_v_len(0,m-1));
        s3 = testlsfitunit_getglserror(n, m, y, w, fmatrix, &c2, _state);
        result = (result&&ae_fp_greater_eq(s2,s1/(1+threshold)))&&ae_fp_greater_eq(s3,s1/(1+threshold));
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Tests whether C is solution of LLS problem
*************************************************************************/
static double testlsfitunit_getglserror(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_vector* c,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    double result;


    result = 0;
    for(i=0; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&fmatrix->ptr.pp_double[i][0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        result = result+ae_sqr(w->ptr.p_double[i]*(v-y->ptr.p_double[i]), _state);
    }
    return result;
}


/*************************************************************************
Subroutine for nonlinear fitting of linear problem

DerAvailable:
* 0     when only function value should be used
* 1     when we can provide gradient/function
* 2     when we can provide Hessian/gradient/function

When something which is not permitted by DerAvailable is requested,
this function sets NLSErrors to True.
*************************************************************************/
static void testlsfitunit_fitlinearnonlinear(ae_int_t m,
     ae_int_t deravailable,
     /* Real    */ ae_matrix* xy,
     lsfitstate* state,
     ae_bool* nlserrors,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;


    while(lsfititeration(state, _state))
    {
        
        /*
         * assume that one and only one of flags is set
         * test that we didn't request hessian in hessian-free setting
         */
        if( deravailable<1&&state->needfg )
        {
            *nlserrors = ae_true;
        }
        if( deravailable<2&&state->needfgh )
        {
            *nlserrors = ae_true;
        }
        i = 0;
        if( state->needf )
        {
            i = i+1;
        }
        if( state->needfg )
        {
            i = i+1;
        }
        if( state->needfgh )
        {
            i = i+1;
        }
        if( i!=1 )
        {
            *nlserrors = ae_true;
        }
        
        /*
         * test that PointIndex is consistent with actual point passed
         */
        for(i=0; i<=m-1; i++)
        {
            *nlserrors = *nlserrors||ae_fp_neq(xy->ptr.pp_double[state->pointindex][i],state->x.ptr.p_double[i]);
        }
        
        /*
         * calculate
         */
        if( state->needf )
        {
            v = ae_v_dotproduct(&state->x.ptr.p_double[0], 1, &state->c.ptr.p_double[0], 1, ae_v_len(0,m-1));
            state->f = v;
            continue;
        }
        if( state->needfg )
        {
            v = ae_v_dotproduct(&state->x.ptr.p_double[0], 1, &state->c.ptr.p_double[0], 1, ae_v_len(0,m-1));
            state->f = v;
            ae_v_move(&state->g.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,m-1));
            continue;
        }
        if( state->needfgh )
        {
            v = ae_v_dotproduct(&state->x.ptr.p_double[0], 1, &state->c.ptr.p_double[0], 1, ae_v_len(0,m-1));
            state->f = v;
            ae_v_move(&state->g.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,m-1));
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    state->h.ptr.pp_double[i][j] = 0;
                }
            }
            continue;
        }
    }
}


/*$ End $*/
