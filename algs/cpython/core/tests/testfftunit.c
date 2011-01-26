

#include <stdafx.h>
#include <stdio.h>
#include "testfftunit.h"


/*$ Declarations $*/
static void testfftunit_reffftc1d(/* Complex */ ae_vector* a,
     ae_int_t n,
     ae_state *_state);
static void testfftunit_reffftc1dinv(/* Complex */ ae_vector* a,
     ae_int_t n,
     ae_state *_state);
static void testfftunit_refinternalcfft(/* Real    */ ae_vector* a,
     ae_int_t nn,
     ae_bool inversefft,
     ae_state *_state);
static void testfftunit_refinternalrfft(/* Real    */ ae_vector* a,
     ae_int_t nn,
     /* Complex */ ae_vector* f,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Test
*************************************************************************/
ae_bool testfft(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t n;
    ae_int_t i;
    ae_int_t k;
    ae_vector a1;
    ae_vector a2;
    ae_vector a3;
    ae_vector r1;
    ae_vector r2;
    ae_vector buf;
    ftplan plan;
    ae_int_t maxn;
    double bidierr;
    double bidirerr;
    double referr;
    double refrerr;
    double reinterr;
    double errtol;
    ae_bool referrors;
    ae_bool bidierrors;
    ae_bool refrerrors;
    ae_bool bidirerrors;
    ae_bool reinterrors;
    ae_bool waserrors;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&a1, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&a2, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&a3, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&r1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&r2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);
    _ftplan_init(&plan, _state, ae_true);

    maxn = 128;
    errtol = 100000*ae_pow(maxn, (double)3/(double)2, _state)*ae_machineepsilon;
    bidierrors = ae_false;
    referrors = ae_false;
    bidirerrors = ae_false;
    refrerrors = ae_false;
    reinterrors = ae_false;
    waserrors = ae_false;
    
    /*
     * Test bi-directional error: norm(x-invFFT(FFT(x)))
     */
    bidierr = 0;
    bidirerr = 0;
    for(n=1; n<=maxn; n++)
    {
        
        /*
         * Complex FFT/invFFT
         */
        ae_vector_set_length(&a1, n, _state);
        ae_vector_set_length(&a2, n, _state);
        ae_vector_set_length(&a3, n, _state);
        for(i=0; i<=n-1; i++)
        {
            a1.ptr.p_complex[i].x = 2*ae_randomreal(_state)-1;
            a1.ptr.p_complex[i].y = 2*ae_randomreal(_state)-1;
            a2.ptr.p_complex[i] = a1.ptr.p_complex[i];
            a3.ptr.p_complex[i] = a1.ptr.p_complex[i];
        }
        fftc1d(&a2, n, _state);
        fftc1dinv(&a2, n, _state);
        fftc1dinv(&a3, n, _state);
        fftc1d(&a3, n, _state);
        for(i=0; i<=n-1; i++)
        {
            bidierr = ae_maxreal(bidierr, ae_c_abs(ae_c_sub(a1.ptr.p_complex[i],a2.ptr.p_complex[i]), _state), _state);
            bidierr = ae_maxreal(bidierr, ae_c_abs(ae_c_sub(a1.ptr.p_complex[i],a3.ptr.p_complex[i]), _state), _state);
        }
        
        /*
         * Real
         */
        ae_vector_set_length(&r1, n, _state);
        ae_vector_set_length(&r2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            r1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            r2.ptr.p_double[i] = r1.ptr.p_double[i];
        }
        fftr1d(&r2, n, &a1, _state);
        ae_v_muld(&r2.ptr.p_double[0], 1, ae_v_len(0,n-1), 0);
        fftr1dinv(&a1, n, &r2, _state);
        for(i=0; i<=n-1; i++)
        {
            bidirerr = ae_maxreal(bidirerr, ae_c_abs(ae_complex_from_d(r1.ptr.p_double[i]-r2.ptr.p_double[i]), _state), _state);
        }
    }
    bidierrors = bidierrors||ae_fp_greater(bidierr,errtol);
    bidirerrors = bidirerrors||ae_fp_greater(bidirerr,errtol);
    
    /*
     * Test against reference O(N^2) implementation
     */
    referr = 0;
    refrerr = 0;
    for(n=1; n<=maxn; n++)
    {
        
        /*
         * Complex FFT
         */
        ae_vector_set_length(&a1, n, _state);
        ae_vector_set_length(&a2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            a1.ptr.p_complex[i].x = 2*ae_randomreal(_state)-1;
            a1.ptr.p_complex[i].y = 2*ae_randomreal(_state)-1;
            a2.ptr.p_complex[i] = a1.ptr.p_complex[i];
        }
        fftc1d(&a1, n, _state);
        testfftunit_reffftc1d(&a2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            referr = ae_maxreal(referr, ae_c_abs(ae_c_sub(a1.ptr.p_complex[i],a2.ptr.p_complex[i]), _state), _state);
        }
        
        /*
         * Complex inverse FFT
         */
        ae_vector_set_length(&a1, n, _state);
        ae_vector_set_length(&a2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            a1.ptr.p_complex[i].x = 2*ae_randomreal(_state)-1;
            a1.ptr.p_complex[i].y = 2*ae_randomreal(_state)-1;
            a2.ptr.p_complex[i] = a1.ptr.p_complex[i];
        }
        fftc1dinv(&a1, n, _state);
        testfftunit_reffftc1dinv(&a2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            referr = ae_maxreal(referr, ae_c_abs(ae_c_sub(a1.ptr.p_complex[i],a2.ptr.p_complex[i]), _state), _state);
        }
        
        /*
         * Real forward/inverse FFT:
         * * calculate and check forward FFT
         * * use precalculated FFT to check backward FFT
         *   fill unused parts of frequencies array with random numbers
         *   to ensure that they are not really used
         */
        ae_vector_set_length(&r1, n, _state);
        ae_vector_set_length(&r2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            r1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            r2.ptr.p_double[i] = r1.ptr.p_double[i];
        }
        fftr1d(&r1, n, &a1, _state);
        testfftunit_refinternalrfft(&r2, n, &a2, _state);
        for(i=0; i<=n-1; i++)
        {
            refrerr = ae_maxreal(refrerr, ae_c_abs(ae_c_sub(a1.ptr.p_complex[i],a2.ptr.p_complex[i]), _state), _state);
        }
        ae_vector_set_length(&a3, ae_ifloor((double)n/(double)2, _state)+1, _state);
        for(i=0; i<=ae_ifloor((double)n/(double)2, _state); i++)
        {
            a3.ptr.p_complex[i] = a2.ptr.p_complex[i];
        }
        a3.ptr.p_complex[0].y = 2*ae_randomreal(_state)-1;
        if( n%2==0 )
        {
            a3.ptr.p_complex[ae_ifloor((double)n/(double)2, _state)].y = 2*ae_randomreal(_state)-1;
        }
        for(i=0; i<=n-1; i++)
        {
            r1.ptr.p_double[i] = 0;
        }
        fftr1dinv(&a3, n, &r1, _state);
        for(i=0; i<=n-1; i++)
        {
            refrerr = ae_maxreal(refrerr, ae_fabs(r2.ptr.p_double[i]-r1.ptr.p_double[i], _state), _state);
        }
    }
    referrors = referrors||ae_fp_greater(referr,errtol);
    refrerrors = refrerrors||ae_fp_greater(refrerr,errtol);
    
    /*
     * test internal real even FFT
     */
    reinterr = 0;
    for(k=1; k<=maxn/2; k++)
    {
        n = 2*k;
        
        /*
         * Real forward FFT
         */
        ae_vector_set_length(&r1, n, _state);
        ae_vector_set_length(&r2, n, _state);
        for(i=0; i<=n-1; i++)
        {
            r1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            r2.ptr.p_double[i] = r1.ptr.p_double[i];
        }
        ftbasegeneratecomplexfftplan(n/2, &plan, _state);
        ae_vector_set_length(&buf, n, _state);
        fftr1dinternaleven(&r1, n, &buf, &plan, _state);
        testfftunit_refinternalrfft(&r2, n, &a2, _state);
        reinterr = ae_maxreal(reinterr, ae_fabs(r1.ptr.p_double[0]-a2.ptr.p_complex[0].x, _state), _state);
        reinterr = ae_maxreal(reinterr, ae_fabs(r1.ptr.p_double[1]-a2.ptr.p_complex[n/2].x, _state), _state);
        for(i=1; i<=n/2-1; i++)
        {
            reinterr = ae_maxreal(reinterr, ae_fabs(r1.ptr.p_double[2*i+0]-a2.ptr.p_complex[i].x, _state), _state);
            reinterr = ae_maxreal(reinterr, ae_fabs(r1.ptr.p_double[2*i+1]-a2.ptr.p_complex[i].y, _state), _state);
        }
        
        /*
         * Real backward FFT
         */
        ae_vector_set_length(&r1, n, _state);
        for(i=0; i<=n-1; i++)
        {
            r1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        ae_vector_set_length(&a2, ae_ifloor((double)n/(double)2, _state)+1, _state);
        a2.ptr.p_complex[0] = ae_complex_from_d(r1.ptr.p_double[0]);
        for(i=1; i<=ae_ifloor((double)n/(double)2, _state)-1; i++)
        {
            a2.ptr.p_complex[i].x = r1.ptr.p_double[2*i+0];
            a2.ptr.p_complex[i].y = r1.ptr.p_double[2*i+1];
        }
        a2.ptr.p_complex[ae_ifloor((double)n/(double)2, _state)] = ae_complex_from_d(r1.ptr.p_double[1]);
        ftbasegeneratecomplexfftplan(n/2, &plan, _state);
        ae_vector_set_length(&buf, n, _state);
        fftr1dinvinternaleven(&r1, n, &buf, &plan, _state);
        fftr1dinv(&a2, n, &r2, _state);
        for(i=0; i<=n-1; i++)
        {
            reinterr = ae_maxreal(reinterr, ae_fabs(r1.ptr.p_double[i]-r2.ptr.p_double[i], _state), _state);
        }
    }
    reinterrors = reinterrors||ae_fp_greater(reinterr,errtol);
    
    /*
     * end
     */
    waserrors = (((bidierrors||bidirerrors)||referrors)||refrerrors)||reinterrors;
    if( !silent )
    {
        printf("TESTING FFT\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* BI-DIRECTIONAL COMPLEX TEST:            ");
        if( bidierrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE COMPLEX FFT:          ");
        if( referrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* BI-DIRECTIONAL REAL TEST:               ");
        if( bidirerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* AGAINST REFERENCE REAL FFT:             ");
        if( refrerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* INTERNAL EVEN FFT:                      ");
        if( reinterrors )
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
    }
    result = !waserrors;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Reference FFT
*************************************************************************/
static void testfftunit_reffftc1d(/* Complex */ ae_vector* a,
     ae_int_t n,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);

    ae_assert(n>0, "FFTC1D: incorrect N!", _state);
    ae_vector_set_length(&buf, 2*n, _state);
    for(i=0; i<=n-1; i++)
    {
        buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
        buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
    }
    testfftunit_refinternalcfft(&buf, n, ae_false, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_complex[i].x = buf.ptr.p_double[2*i+0];
        a->ptr.p_complex[i].y = buf.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Reference inverse FFT
*************************************************************************/
static void testfftunit_reffftc1dinv(/* Complex */ ae_vector* a,
     ae_int_t n,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);

    ae_assert(n>0, "FFTC1DInv: incorrect N!", _state);
    ae_vector_set_length(&buf, 2*n, _state);
    for(i=0; i<=n-1; i++)
    {
        buf.ptr.p_double[2*i+0] = a->ptr.p_complex[i].x;
        buf.ptr.p_double[2*i+1] = a->ptr.p_complex[i].y;
    }
    testfftunit_refinternalcfft(&buf, n, ae_true, _state);
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_complex[i].x = buf.ptr.p_double[2*i+0];
        a->ptr.p_complex[i].y = buf.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal complex FFT stub.
Uses straightforward formula with O(N^2) complexity.
*************************************************************************/
static void testfftunit_refinternalcfft(/* Real    */ ae_vector* a,
     ae_int_t nn,
     ae_bool inversefft,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double hre;
    double him;
    double c;
    double s;
    double re;
    double im;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&tmp, 2*nn-1+1, _state);
    if( !inversefft )
    {
        for(i=0; i<=nn-1; i++)
        {
            hre = 0;
            him = 0;
            for(k=0; k<=nn-1; k++)
            {
                re = a->ptr.p_double[2*k];
                im = a->ptr.p_double[2*k+1];
                c = ae_cos(-2*ae_pi*k*i/nn, _state);
                s = ae_sin(-2*ae_pi*k*i/nn, _state);
                hre = hre+c*re-s*im;
                him = him+c*im+s*re;
            }
            tmp.ptr.p_double[2*i] = hre;
            tmp.ptr.p_double[2*i+1] = him;
        }
        for(i=0; i<=2*nn-1; i++)
        {
            a->ptr.p_double[i] = tmp.ptr.p_double[i];
        }
    }
    else
    {
        for(k=0; k<=nn-1; k++)
        {
            hre = 0;
            him = 0;
            for(i=0; i<=nn-1; i++)
            {
                re = a->ptr.p_double[2*i];
                im = a->ptr.p_double[2*i+1];
                c = ae_cos(2*ae_pi*k*i/nn, _state);
                s = ae_sin(2*ae_pi*k*i/nn, _state);
                hre = hre+c*re-s*im;
                him = him+c*im+s*re;
            }
            tmp.ptr.p_double[2*k] = hre/nn;
            tmp.ptr.p_double[2*k+1] = him/nn;
        }
        for(i=0; i<=2*nn-1; i++)
        {
            a->ptr.p_double[i] = tmp.ptr.p_double[i];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal real FFT stub.
Uses straightforward formula with O(N^2) complexity.
*************************************************************************/
static void testfftunit_refinternalrfft(/* Real    */ ae_vector* a,
     ae_int_t nn,
     /* Complex */ ae_vector* f,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(f);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&tmp, 2*nn-1+1, _state);
    for(i=0; i<=nn-1; i++)
    {
        tmp.ptr.p_double[2*i] = a->ptr.p_double[i];
        tmp.ptr.p_double[2*i+1] = 0;
    }
    testfftunit_refinternalcfft(&tmp, nn, ae_false, _state);
    ae_vector_set_length(f, nn, _state);
    for(i=0; i<=nn-1; i++)
    {
        f->ptr.p_complex[i].x = tmp.ptr.p_double[2*i+0];
        f->ptr.p_complex[i].y = tmp.ptr.p_double[2*i+1];
    }
    ae_frame_leave(_state);
}


/*$ End $*/
