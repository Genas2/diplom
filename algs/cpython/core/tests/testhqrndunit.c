

#include <stdafx.h>
#include <stdio.h>
#include "testhqrndunit.h"


/*$ Declarations $*/
static void testhqrndunit_calculatemv(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* mean,
     double* means,
     double* stddev,
     double* stddevs,
     ae_state *_state);
static void testhqrndunit_unsetstate(hqrndstate* state, ae_state *_state);


/*$ Body $*/


ae_bool testhqrnd(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_int_t samplesize;
    double sigmathreshold;
    ae_int_t passcount;
    ae_int_t n;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t s1;
    ae_int_t s2;
    ae_int_t i1;
    ae_int_t i2;
    double r1;
    double r2;
    ae_vector x;
    double mean;
    double means;
    double stddev;
    double stddevs;
    double lambdav;
    ae_bool seederrors;
    ae_bool urerrors;
    double ursigmaerr;
    ae_bool uierrors;
    double uisigmaerr;
    ae_bool normerrors;
    double normsigmaerr;
    ae_bool experrors;
    double expsigmaerr;
    hqrndstate state;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    _hqrndstate_init(&state, _state, ae_true);

    waserrors = ae_false;
    sigmathreshold = 7;
    samplesize = 100000;
    passcount = 50;
    seederrors = ae_false;
    urerrors = ae_false;
    uierrors = ae_false;
    normerrors = ae_false;
    experrors = ae_false;
    ae_vector_set_length(&x, samplesize-1+1, _state);
    
    /*
     * Test seed errors
     */
    for(pass=1; pass<=passcount; pass++)
    {
        s1 = 1+ae_randominteger(32000, _state);
        s2 = 1+ae_randominteger(32000, _state);
        testhqrndunit_unsetstate(&state, _state);
        hqrndseed(s1, s2, &state, _state);
        i1 = hqrnduniformi(&state, 100, _state);
        testhqrndunit_unsetstate(&state, _state);
        hqrndseed(s1, s2, &state, _state);
        i2 = hqrnduniformi(&state, 100, _state);
        seederrors = seederrors||i1!=i2;
        testhqrndunit_unsetstate(&state, _state);
        hqrndseed(s1, s2, &state, _state);
        r1 = hqrnduniformr(&state, _state);
        testhqrndunit_unsetstate(&state, _state);
        hqrndseed(s1, s2, &state, _state);
        r2 = hqrnduniformr(&state, _state);
        seederrors = seederrors||ae_fp_neq(r1,r2);
    }
    
    /*
     * Test HQRNDRandomize() and real uniform generator
     */
    testhqrndunit_unsetstate(&state, _state);
    hqrndrandomize(&state, _state);
    ursigmaerr = 0;
    for(i=0; i<=samplesize-1; i++)
    {
        x.ptr.p_double[i] = hqrnduniformr(&state, _state);
    }
    for(i=0; i<=samplesize-1; i++)
    {
        urerrors = (urerrors||ae_fp_less_eq(x.ptr.p_double[i],0))||ae_fp_greater_eq(x.ptr.p_double[i],1);
    }
    testhqrndunit_calculatemv(&x, samplesize, &mean, &means, &stddev, &stddevs, _state);
    if( ae_fp_neq(means,0) )
    {
        ursigmaerr = ae_maxreal(ursigmaerr, ae_fabs((mean-0.5)/means, _state), _state);
    }
    else
    {
        urerrors = ae_true;
    }
    if( ae_fp_neq(stddevs,0) )
    {
        ursigmaerr = ae_maxreal(ursigmaerr, ae_fabs((stddev-ae_sqrt((double)1/(double)12, _state))/stddevs, _state), _state);
    }
    else
    {
        urerrors = ae_true;
    }
    urerrors = urerrors||ae_fp_greater(ursigmaerr,sigmathreshold);
    
    /*
     * Test HQRNDRandomize() and integer uniform
     */
    testhqrndunit_unsetstate(&state, _state);
    hqrndrandomize(&state, _state);
    uisigmaerr = 0;
    for(n=2; n<=10; n++)
    {
        for(i=0; i<=samplesize-1; i++)
        {
            x.ptr.p_double[i] = hqrnduniformi(&state, n, _state);
        }
        for(i=0; i<=samplesize-1; i++)
        {
            uierrors = (uierrors||ae_fp_less(x.ptr.p_double[i],0))||ae_fp_greater_eq(x.ptr.p_double[i],n);
        }
        testhqrndunit_calculatemv(&x, samplesize, &mean, &means, &stddev, &stddevs, _state);
        if( ae_fp_neq(means,0) )
        {
            uisigmaerr = ae_maxreal(uisigmaerr, ae_fabs((mean-0.5*(n-1))/means, _state), _state);
        }
        else
        {
            uierrors = ae_true;
        }
        if( ae_fp_neq(stddevs,0) )
        {
            uisigmaerr = ae_maxreal(uisigmaerr, ae_fabs((stddev-ae_sqrt((ae_sqr(n, _state)-1)/12, _state))/stddevs, _state), _state);
        }
        else
        {
            uierrors = ae_true;
        }
    }
    uierrors = uierrors||ae_fp_greater(uisigmaerr,sigmathreshold);
    
    /*
     * Special 'close-to-limit' test on uniformity of integers
     * (straightforward implementation like 'RND mod N' will return
     *  non-uniform numbers for N=2/3*LIMIT)
     */
    testhqrndunit_unsetstate(&state, _state);
    hqrndrandomize(&state, _state);
    uisigmaerr = 0;
    n = 1431655708;
    for(i=0; i<=samplesize-1; i++)
    {
        x.ptr.p_double[i] = hqrnduniformi(&state, n, _state);
    }
    for(i=0; i<=samplesize-1; i++)
    {
        uierrors = (uierrors||ae_fp_less(x.ptr.p_double[i],0))||ae_fp_greater_eq(x.ptr.p_double[i],n);
    }
    testhqrndunit_calculatemv(&x, samplesize, &mean, &means, &stddev, &stddevs, _state);
    if( ae_fp_neq(means,0) )
    {
        uisigmaerr = ae_maxreal(uisigmaerr, ae_fabs((mean-0.5*(n-1))/means, _state), _state);
    }
    else
    {
        uierrors = ae_true;
    }
    if( ae_fp_neq(stddevs,0) )
    {
        uisigmaerr = ae_maxreal(uisigmaerr, ae_fabs((stddev-ae_sqrt((ae_sqr(n, _state)-1)/12, _state))/stddevs, _state), _state);
    }
    else
    {
        uierrors = ae_true;
    }
    uierrors = uierrors||ae_fp_greater(uisigmaerr,sigmathreshold);
    
    /*
     * Test normal
     */
    testhqrndunit_unsetstate(&state, _state);
    hqrndrandomize(&state, _state);
    normsigmaerr = 0;
    i = 0;
    while(i<samplesize)
    {
        hqrndnormal2(&state, &r1, &r2, _state);
        x.ptr.p_double[i] = r1;
        if( i+1<samplesize )
        {
            x.ptr.p_double[i+1] = r2;
        }
        i = i+2;
    }
    testhqrndunit_calculatemv(&x, samplesize, &mean, &means, &stddev, &stddevs, _state);
    if( ae_fp_neq(means,0) )
    {
        normsigmaerr = ae_maxreal(normsigmaerr, ae_fabs((mean-0)/means, _state), _state);
    }
    else
    {
        normerrors = ae_true;
    }
    if( ae_fp_neq(stddevs,0) )
    {
        normsigmaerr = ae_maxreal(normsigmaerr, ae_fabs((stddev-1)/stddevs, _state), _state);
    }
    else
    {
        normerrors = ae_true;
    }
    normerrors = normerrors||ae_fp_greater(normsigmaerr,sigmathreshold);
    
    /*
     * Test exponential
     */
    testhqrndunit_unsetstate(&state, _state);
    hqrndrandomize(&state, _state);
    expsigmaerr = 0;
    lambdav = 2+5*ae_randomreal(_state);
    for(i=0; i<=samplesize-1; i++)
    {
        x.ptr.p_double[i] = hqrndexponential(&state, lambdav, _state);
    }
    for(i=0; i<=samplesize-1; i++)
    {
        uierrors = uierrors||ae_fp_less(x.ptr.p_double[i],0);
    }
    testhqrndunit_calculatemv(&x, samplesize, &mean, &means, &stddev, &stddevs, _state);
    if( ae_fp_neq(means,0) )
    {
        expsigmaerr = ae_maxreal(expsigmaerr, ae_fabs((mean-1.0/lambdav)/means, _state), _state);
    }
    else
    {
        experrors = ae_true;
    }
    if( ae_fp_neq(stddevs,0) )
    {
        expsigmaerr = ae_maxreal(expsigmaerr, ae_fabs((stddev-1.0/lambdav)/stddevs, _state), _state);
    }
    else
    {
        experrors = ae_true;
    }
    experrors = experrors||ae_fp_greater(expsigmaerr,sigmathreshold);
    
    /*
     * Final report
     */
    waserrors = (((seederrors||urerrors)||uierrors)||normerrors)||experrors;
    if( !silent )
    {
        printf("RNG TEST\n");
        printf("SEED TEST:                               ");
        if( !seederrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("UNIFORM CONTINUOUS:                      ");
        if( !urerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("UNIFORM INTEGER:                         ");
        if( !uierrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("NORMAL:                                  ");
        if( !normerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("EXPONENTIAL:                             ");
        if( !experrors )
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
    ae_frame_leave(_state);
    return result;
}


static void testhqrndunit_calculatemv(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* mean,
     double* means,
     double* stddev,
     double* stddevs,
     ae_state *_state)
{
    ae_int_t i;
    double v1;
    double v2;
    double variance;

    *mean = 0;
    *means = 0;
    *stddev = 0;
    *stddevs = 0;

    *mean = 0;
    *means = 1;
    *stddev = 0;
    *stddevs = 1;
    variance = 0;
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Mean
     */
    for(i=0; i<=n-1; i++)
    {
        *mean = *mean+x->ptr.p_double[i];
    }
    *mean = *mean/n;
    
    /*
     * Variance (using corrected two-pass algorithm)
     */
    if( n!=1 )
    {
        v1 = 0;
        for(i=0; i<=n-1; i++)
        {
            v1 = v1+ae_sqr(x->ptr.p_double[i]-(*mean), _state);
        }
        v2 = 0;
        for(i=0; i<=n-1; i++)
        {
            v2 = v2+(x->ptr.p_double[i]-(*mean));
        }
        v2 = ae_sqr(v2, _state)/n;
        variance = (v1-v2)/(n-1);
        if( ae_fp_less(variance,0) )
        {
            variance = 0;
        }
        *stddev = ae_sqrt(variance, _state);
    }
    
    /*
     * Errors
     */
    *means = *stddev/ae_sqrt(n, _state);
    *stddevs = *stddev*ae_sqrt(2, _state)/ae_sqrt(n-1, _state);
}


/*************************************************************************
Unsets HQRNDState structure
*************************************************************************/
static void testhqrndunit_unsetstate(hqrndstate* state, ae_state *_state)
{


    state->s1 = 0;
    state->s2 = 0;
    state->v = 0;
    state->magicv = 0;
}


/*$ End $*/
