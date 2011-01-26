

#include <stdafx.h>
#include <stdio.h>
#include "testkmeansunit.h"


/*$ Declarations $*/
static void testkmeansunit_simpletest1(ae_int_t nvars,
     ae_int_t nc,
     ae_int_t passcount,
     ae_bool* converrors,
     ae_bool* othererrors,
     ae_bool* simpleerrors,
     ae_state *_state);
static void testkmeansunit_restartstest(ae_bool* converrors,
     ae_bool* restartserrors,
     ae_state *_state);
static double testkmeansunit_rnormal(ae_state *_state);
static void testkmeansunit_rsphere(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t i,
     ae_state *_state);


/*$ Body $*/


ae_bool testkmeans(ae_bool silent, ae_state *_state)
{
    ae_int_t nf;
    ae_int_t maxnf;
    ae_int_t nc;
    ae_int_t maxnc;
    ae_int_t passcount;
    ae_bool waserrors;
    ae_bool converrors;
    ae_bool simpleerrors;
    ae_bool complexerrors;
    ae_bool othererrors;
    ae_bool restartserrors;
    ae_bool result;


    
    /*
     * Primary settings
     */
    maxnf = 5;
    maxnc = 5;
    passcount = 10;
    waserrors = ae_false;
    converrors = ae_false;
    othererrors = ae_false;
    simpleerrors = ae_false;
    complexerrors = ae_false;
    restartserrors = ae_false;
    
    /*
     *
     */
    for(nf=1; nf<=maxnf; nf++)
    {
        for(nc=1; nc<=maxnc; nc++)
        {
            testkmeansunit_simpletest1(nf, nc, passcount, &converrors, &othererrors, &simpleerrors, _state);
        }
    }
    testkmeansunit_restartstest(&converrors, &restartserrors, _state);
    
    /*
     * Final report
     */
    waserrors = (((converrors||othererrors)||simpleerrors)||complexerrors)||restartserrors;
    if( !silent )
    {
        printf("K-MEANS TEST\n");
        printf("TOTAL RESULTS:                           ");
        if( !waserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* CONVERGENCE:                           ");
        if( !converrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* SIMPLE TASKS:                          ");
        if( !simpleerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* COMPLEX TASKS:                         ");
        if( !complexerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* OTHER PROPERTIES:                      ");
        if( !othererrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* RESTARTS PROPERTIES:                   ");
        if( !restartserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST SUMMARY: FAILED\n");
        }
        else
        {
            printf("TEST SUMMARY: PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Simple test 1: ellipsoid in NF-dimensional space.
compare k-means centers with random centers
*************************************************************************/
static void testkmeansunit_simpletest1(ae_int_t nvars,
     ae_int_t nc,
     ae_int_t passcount,
     ae_bool* converrors,
     ae_bool* othererrors,
     ae_bool* simpleerrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t npoints;
    ae_int_t majoraxis;
    ae_matrix xy;
    ae_vector tmp;
    double v;
    ae_int_t i;
    ae_int_t j;
    ae_int_t info;
    ae_matrix c;
    ae_vector xyc;
    ae_int_t pass;
    ae_int_t restarts;
    double ekmeans;
    double erandom;
    double dclosest;
    ae_int_t cclosest;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xyc, 0, DT_INT, _state, ae_true);

    npoints = nc*100;
    restarts = 5;
    passcount = 10;
    ae_vector_set_length(&tmp, nvars-1+1, _state);
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Fill
         */
        ae_matrix_set_length(&xy, npoints-1+1, nvars-1+1, _state);
        majoraxis = ae_randominteger(nvars, _state);
        for(i=0; i<=npoints-1; i++)
        {
            testkmeansunit_rsphere(&xy, nvars, i, _state);
            xy.ptr.pp_double[i][majoraxis] = nc*xy.ptr.pp_double[i][majoraxis];
        }
        
        /*
         * Test
         */
        kmeansgenerate(&xy, npoints, nvars, nc, restarts, &info, &c, &xyc, _state);
        if( info<0 )
        {
            *converrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Test that XYC is correct mapping to cluster centers
         */
        for(i=0; i<=npoints-1; i++)
        {
            cclosest = -1;
            dclosest = ae_maxrealnumber;
            for(j=0; j<=nc-1; j++)
            {
                ae_v_move(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                ae_v_sub(&tmp.ptr.p_double[0], 1, &c.ptr.pp_double[0][j], c.stride, ae_v_len(0,nvars-1));
                v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                if( ae_fp_less(v,dclosest) )
                {
                    cclosest = j;
                    dclosest = v;
                }
            }
            if( cclosest!=xyc.ptr.p_int[i] )
            {
                *othererrors = ae_true;
                ae_frame_leave(_state);
                return;
            }
        }
        
        /*
         * Use first NC rows of XY as random centers
         * (XY is totally random, so it is as good as any other choice).
         *
         * Compare potential functions.
         */
        ekmeans = 0;
        for(i=0; i<=npoints-1; i++)
        {
            ae_v_move(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            ae_v_sub(&tmp.ptr.p_double[0], 1, &c.ptr.pp_double[0][xyc.ptr.p_int[i]], c.stride, ae_v_len(0,nvars-1));
            v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
            ekmeans = ekmeans+v;
        }
        erandom = 0;
        for(i=0; i<=npoints-1; i++)
        {
            dclosest = ae_maxrealnumber;
            for(j=0; j<=nc-1; j++)
            {
                ae_v_move(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                ae_v_sub(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
                v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                if( ae_fp_less(v,dclosest) )
                {
                    dclosest = v;
                }
            }
            erandom = erandom+v;
        }
        if( ae_fp_less(erandom,ekmeans) )
        {
            *simpleerrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This non-deterministic test checks that Restarts>1 significantly  improves
quality of results.

Subroutine generates random task 3 unit balls in 2D, each with 20  points,
separated by 5 units wide gaps, and solves it  with  Restarts=1  and  with
Restarts=5. Potential functions are compared,  outcome  of  the  trial  is
either 0 or 1 (depending on what is better).

Sequence of 1000 such tasks is  solved.  If  Restarts>1  actually  improve
quality of solution, sum of outcome will be non-binomial.  If  it  doesn't
matter, it will be binomially distributed.

P.S. This test was added after report from Gianluca  Borello  who  noticed
error in the handling of multiple restarts.
*************************************************************************/
static void testkmeansunit_restartstest(ae_bool* converrors,
     ae_bool* restartserrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t npoints;
    ae_int_t nvars;
    ae_int_t nclusters;
    ae_int_t clustersize;
    ae_int_t restarts;
    ae_int_t passcount;
    double sigmathreshold;
    double p;
    double s;
    ae_matrix xy;
    ae_matrix ca;
    ae_matrix cb;
    ae_vector xyca;
    ae_vector xycb;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    ae_int_t info;
    ae_int_t pass;
    double ea;
    double eb;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ca, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cb, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xyca, 0, DT_INT, _state, ae_true);
    ae_vector_init(&xycb, 0, DT_INT, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    restarts = 5;
    passcount = 1000;
    clustersize = 20;
    nclusters = 3;
    nvars = 2;
    npoints = nclusters*clustersize;
    sigmathreshold = 5;
    ae_matrix_set_length(&xy, npoints, nvars, _state);
    ae_vector_set_length(&tmp, nvars, _state);
    p = 0;
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Fill
         */
        for(i=0; i<=npoints-1; i++)
        {
            testkmeansunit_rsphere(&xy, nvars, i, _state);
            for(j=0; j<=nvars-1; j++)
            {
                xy.ptr.pp_double[i][j] = xy.ptr.pp_double[i][j]+(double)i/(double)clustersize*5;
            }
        }
        
        /*
         * Test: Restarts=1
         */
        kmeansgenerate(&xy, npoints, nvars, nclusters, 1, &info, &ca, &xyca, _state);
        if( info<0 )
        {
            *converrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        ea = 0;
        for(i=0; i<=npoints-1; i++)
        {
            ae_v_move(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            ae_v_sub(&tmp.ptr.p_double[0], 1, &ca.ptr.pp_double[0][xyca.ptr.p_int[i]], ca.stride, ae_v_len(0,nvars-1));
            v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
            ea = ea+v;
        }
        
        /*
         * Test: Restarts>1
         */
        kmeansgenerate(&xy, npoints, nvars, nclusters, restarts, &info, &cb, &xycb, _state);
        if( info<0 )
        {
            *converrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        eb = 0;
        for(i=0; i<=npoints-1; i++)
        {
            ae_v_move(&tmp.ptr.p_double[0], 1, &xy.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            ae_v_sub(&tmp.ptr.p_double[0], 1, &cb.ptr.pp_double[0][xycb.ptr.p_int[i]], cb.stride, ae_v_len(0,nvars-1));
            v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
            eb = eb+v;
        }
        
        /*
         * Calculate statistic.
         */
        if( ae_fp_less(ea,eb) )
        {
            p = p+1;
        }
        if( ae_fp_eq(ea,eb) )
        {
            p = p+0.5;
        }
    }
    
    /*
     * If Restarts doesn't influence quality of centers found, P must be
     * binomially distributed random value with mean 0.5*PassCount and
     * standard deviation Sqrt(PassCount/4).
     *
     * If Restarts do influence quality of solution, P must be significantly
     * lower than 0.5*PassCount.
     */
    s = (p-0.5*passcount)/ae_sqrt((double)passcount/(double)4, _state);
    *restartserrors = *restartserrors||ae_fp_greater(s,-sigmathreshold);
    ae_frame_leave(_state);
}


/*************************************************************************
Random normal number
*************************************************************************/
static double testkmeansunit_rnormal(ae_state *_state)
{
    double u;
    double v;
    double s;
    double x1;
    double x2;
    double result;


    for(;;)
    {
        u = 2*ae_randomreal(_state)-1;
        v = 2*ae_randomreal(_state)-1;
        s = ae_sqr(u, _state)+ae_sqr(v, _state);
        if( ae_fp_greater(s,0)&&ae_fp_less(s,1) )
        {
            s = ae_sqrt(-2*ae_log(s, _state)/s, _state);
            x1 = u*s;
            x2 = v*s;
            break;
        }
    }
    result = x1;
    return result;
}


/*************************************************************************
Random point from sphere
*************************************************************************/
static void testkmeansunit_rsphere(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t i,
     ae_state *_state)
{
    ae_int_t j;
    double v;


    for(j=0; j<=n-1; j++)
    {
        xy->ptr.pp_double[i][j] = testkmeansunit_rnormal(_state);
    }
    v = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
    v = ae_randomreal(_state)/ae_sqrt(v, _state);
    ae_v_muld(&xy->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
}


/*$ End $*/
