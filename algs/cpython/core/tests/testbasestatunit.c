

#include <stdafx.h>
#include <stdio.h>
#include "testbasestatunit.h"


/*$ Declarations $*/


/*$ Body $*/


ae_bool testbasestat(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_bool s1errors;
    ae_bool covcorrerrors;
    double threshold;
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    ae_int_t kx;
    ae_int_t ky;
    ae_int_t ctype;
    ae_int_t cidxx;
    ae_int_t cidxy;
    ae_vector x;
    ae_vector y;
    ae_matrix mx;
    ae_matrix my;
    ae_matrix cc;
    ae_matrix cp;
    ae_matrix cs;
    double mean;
    double variance;
    double skewness;
    double kurtosis;
    double adev;
    double median;
    double pv;
    double v;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&mx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&my, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cc, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cp, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cs, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * Primary settings
     */
    waserrors = ae_false;
    s1errors = ae_false;
    covcorrerrors = ae_false;
    threshold = 1000*ae_machineepsilon;
    
    /*
     * * prepare X and Y - two test samples
     * * test 1-sample coefficients
     */
    n = 10;
    ae_vector_set_length(&x, n, _state);
    for(i=0; i<=n-1; i++)
    {
        x.ptr.p_double[i] = ae_sqr(i, _state);
    }
    samplemoments(&x, n, &mean, &variance, &skewness, &kurtosis, _state);
    s1errors = s1errors||ae_fp_greater(ae_fabs(mean-28.5, _state),0.001);
    s1errors = s1errors||ae_fp_greater(ae_fabs(variance-801.1667, _state),0.001);
    s1errors = s1errors||ae_fp_greater(ae_fabs(skewness-0.5751, _state),0.001);
    s1errors = s1errors||ae_fp_greater(ae_fabs(kurtosis+1.2666, _state),0.001);
    sampleadev(&x, n, &adev, _state);
    s1errors = s1errors||ae_fp_greater(ae_fabs(adev-23.2000, _state),0.001);
    samplemedian(&x, n, &median, _state);
    s1errors = s1errors||ae_fp_greater(ae_fabs(median-0.5*(16+25), _state),0.001);
    for(i=0; i<=n-1; i++)
    {
        samplepercentile(&x, n, (double)i/(double)(n-1), &pv, _state);
        s1errors = s1errors||ae_fp_greater(ae_fabs(pv-x.ptr.p_double[i], _state),0.001);
    }
    samplepercentile(&x, n, 0.5, &pv, _state);
    s1errors = s1errors||ae_fp_greater(ae_fabs(pv-0.5*(16+25), _state),0.001);
    
    /*
     * test covariance/correlation:
     * * 2-sample coefficients
     *
     * We generate random matrices MX and MY
     */
    n = 10;
    ae_vector_set_length(&x, n, _state);
    ae_vector_set_length(&y, n, _state);
    for(i=0; i<=n-1; i++)
    {
        x.ptr.p_double[i] = ae_sqr(i, _state);
        y.ptr.p_double[i] = i;
    }
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(pearsoncorr2(&x, &y, n, _state)-0.9627, _state),0.0001);
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(spearmancorr2(&x, &y, n, _state)-1.0000, _state),0.0001);
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(cov2(&x, &y, n, _state)-82.5000, _state),0.0001);
    for(i=0; i<=n-1; i++)
    {
        x.ptr.p_double[i] = ae_sqr(i-0.5*n, _state);
        y.ptr.p_double[i] = i;
    }
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(pearsoncorr2(&x, &y, n, _state)+0.3676, _state),0.0001);
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(spearmancorr2(&x, &y, n, _state)+0.2761, _state),0.0001);
    covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(cov2(&x, &y, n, _state)+9.1667, _state),0.0001);
    
    /*
     * test covariance/correlation:
     * * matrix covariance/correlation
     * * matrix cross-covariance/cross-correlation
     *
     * We generate random matrices MX and MY which contain KX (KY)
     * columns, all except one are random, one of them is constant.
     * We test that function (a) do not crash on constant column,
     * and (b) return variances/correlations that are exactly zero
     * for this column.
     *
     * CType control variable controls type of constant: 0 - no constant
     * column, 1 - zero column, 2 - nonzero column with value whose
     * binary representation contains many non-zero bits. Using such
     * type of constant column we are able to ensure than even in the
     * presense of roundoff error functions correctly detect constant
     * columns.
     */
    for(n=0; n<=10; n++)
    {
        if( n>0 )
        {
            ae_vector_set_length(&x, n, _state);
            ae_vector_set_length(&y, n, _state);
        }
        for(ctype=0; ctype<=2; ctype++)
        {
            for(kx=1; kx<=10; kx++)
            {
                for(ky=1; ky<=10; ky++)
                {
                    
                    /*
                     * Fill matrices, add constant column (when CType=1 or =2)
                     */
                    if( n>0 )
                    {
                        ae_matrix_set_length(&mx, n, kx, _state);
                        ae_matrix_set_length(&my, n, ky, _state);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=kx-1; j++)
                            {
                                mx.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                            }
                            for(j=0; j<=ky-1; j++)
                            {
                                my.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                            }
                        }
                        if( ctype==1 )
                        {
                            cidxx = ae_randominteger(kx, _state);
                            cidxy = ae_randominteger(ky, _state);
                            for(i=0; i<=n-1; i++)
                            {
                                mx.ptr.pp_double[i][cidxx] = 0.0;
                                my.ptr.pp_double[i][cidxy] = 0.0;
                            }
                        }
                        if( ctype==2 )
                        {
                            cidxx = ae_randominteger(kx, _state);
                            cidxy = ae_randominteger(ky, _state);
                            v = ae_sqrt((double)(ae_randominteger(kx, _state)+1)/(double)kx, _state);
                            for(i=0; i<=n-1; i++)
                            {
                                mx.ptr.pp_double[i][cidxx] = v;
                                my.ptr.pp_double[i][cidxy] = v;
                            }
                        }
                    }
                    
                    /*
                     * test covariance/correlation matrix using
                     * 2-sample functions as reference point.
                     *
                     * We also test that coefficients for constant variables
                     * are exactly zero.
                     */
                    covm(&mx, n, kx, &cc, _state);
                    pearsoncorrm(&mx, n, kx, &cp, _state);
                    spearmancorrm(&mx, n, kx, &cs, _state);
                    for(i=0; i<=kx-1; i++)
                    {
                        for(j=0; j<=kx-1; j++)
                        {
                            if( n>0 )
                            {
                                ae_v_move(&x.ptr.p_double[0], 1, &mx.ptr.pp_double[0][i], mx.stride, ae_v_len(0,n-1));
                                ae_v_move(&y.ptr.p_double[0], 1, &mx.ptr.pp_double[0][j], mx.stride, ae_v_len(0,n-1));
                            }
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(cov2(&x, &y, n, _state)-cc.ptr.pp_double[i][j], _state),threshold);
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(pearsoncorr2(&x, &y, n, _state)-cp.ptr.pp_double[i][j], _state),threshold);
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(spearmancorr2(&x, &y, n, _state)-cs.ptr.pp_double[i][j], _state),threshold);
                        }
                    }
                    if( ctype!=0&&n>0 )
                    {
                        for(i=0; i<=kx-1; i++)
                        {
                            covcorrerrors = covcorrerrors||ae_fp_neq(cc.ptr.pp_double[i][cidxx],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cc.ptr.pp_double[cidxx][i],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cp.ptr.pp_double[i][cidxx],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cp.ptr.pp_double[cidxx][i],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cs.ptr.pp_double[i][cidxx],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cs.ptr.pp_double[cidxx][i],0);
                        }
                    }
                    
                    /*
                     * test cross-covariance/cross-correlation matrix using
                     * 2-sample functions as reference point.
                     *
                     * We also test that coefficients for constant variables
                     * are exactly zero.
                     */
                    covm2(&mx, &my, n, kx, ky, &cc, _state);
                    pearsoncorrm2(&mx, &my, n, kx, ky, &cp, _state);
                    spearmancorrm2(&mx, &my, n, kx, ky, &cs, _state);
                    for(i=0; i<=kx-1; i++)
                    {
                        for(j=0; j<=ky-1; j++)
                        {
                            if( n>0 )
                            {
                                ae_v_move(&x.ptr.p_double[0], 1, &mx.ptr.pp_double[0][i], mx.stride, ae_v_len(0,n-1));
                                ae_v_move(&y.ptr.p_double[0], 1, &my.ptr.pp_double[0][j], my.stride, ae_v_len(0,n-1));
                            }
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(cov2(&x, &y, n, _state)-cc.ptr.pp_double[i][j], _state),threshold);
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(pearsoncorr2(&x, &y, n, _state)-cp.ptr.pp_double[i][j], _state),threshold);
                            covcorrerrors = covcorrerrors||ae_fp_greater(ae_fabs(spearmancorr2(&x, &y, n, _state)-cs.ptr.pp_double[i][j], _state),threshold);
                        }
                    }
                    if( ctype!=0&&n>0 )
                    {
                        for(i=0; i<=kx-1; i++)
                        {
                            covcorrerrors = covcorrerrors||ae_fp_neq(cc.ptr.pp_double[i][cidxy],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cp.ptr.pp_double[i][cidxy],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cs.ptr.pp_double[i][cidxy],0);
                        }
                        for(j=0; j<=ky-1; j++)
                        {
                            covcorrerrors = covcorrerrors||ae_fp_neq(cc.ptr.pp_double[cidxx][j],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cp.ptr.pp_double[cidxx][j],0);
                            covcorrerrors = covcorrerrors||ae_fp_neq(cs.ptr.pp_double[cidxx][j],0);
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Final report
     */
    waserrors = s1errors||covcorrerrors;
    if( !silent )
    {
        printf("DESC.STAT TEST\n");
        printf("TOTAL RESULTS:                           ");
        if( !waserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* 1-SAMPLE FUNCTIONALITY:                ");
        if( !s1errors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* CORRELATION/COVARIATION:               ");
        if( !covcorrerrors )
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


/*$ End $*/
