

#include <stdafx.h>
#include <stdio.h>
#include "testtrfacunit.h"


/*$ Declarations $*/
static void testtrfacunit_testcluproblem(/* Complex */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     double threshold,
     ae_bool* err,
     ae_bool* properr,
     ae_state *_state);
static void testtrfacunit_testrluproblem(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     double threshold,
     ae_bool* err,
     ae_bool* properr,
     ae_state *_state);


/*$ Body $*/


ae_bool testtrfac(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix ra;
    ae_matrix ral;
    ae_matrix rau;
    ae_matrix ca;
    ae_matrix cal;
    ae_matrix cau;
    ae_int_t m;
    ae_int_t n;
    ae_int_t mx;
    ae_int_t maxmn;
    ae_int_t i;
    ae_int_t j;
    ae_int_t minij;
    ae_int_t pass;
    ae_complex vc;
    double vr;
    ae_bool waserrors;
    ae_bool spderr;
    ae_bool hpderr;
    ae_bool rerr;
    ae_bool cerr;
    ae_bool properr;
    double threshold;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&ra, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ral, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rau, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ca, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&cal, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&cau, 0, 0, DT_COMPLEX, _state, ae_true);

    rerr = ae_false;
    spderr = ae_false;
    cerr = ae_false;
    hpderr = ae_false;
    properr = ae_false;
    waserrors = ae_false;
    maxmn = 4*ablasblocksize(&ra, _state)+1;
    threshold = 1000*ae_machineepsilon*maxmn;
    
    /*
     * test LU
     */
    for(mx=1; mx<=maxmn; mx++)
    {
        
        /*
         * Initialize N/M, both are <=MX,
         * at least one of them is exactly equal to MX
         */
        n = 1+ae_randominteger(mx, _state);
        m = 1+ae_randominteger(mx, _state);
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            n = mx;
        }
        else
        {
            m = mx;
        }
        
        /*
         * First, test on zero matrix
         */
        ae_matrix_set_length(&ra, m, n, _state);
        ae_matrix_set_length(&ca, m, n, _state);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ra.ptr.pp_double[i][j] = 0;
                ca.ptr.pp_complex[i][j] = ae_complex_from_d(0);
            }
        }
        testtrfacunit_testcluproblem(&ca, m, n, threshold, &cerr, &properr, _state);
        testtrfacunit_testrluproblem(&ra, m, n, threshold, &rerr, &properr, _state);
        
        /*
         * Second, random matrix with moderate condition number
         */
        ae_matrix_set_length(&ra, m, n, _state);
        ae_matrix_set_length(&ca, m, n, _state);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ra.ptr.pp_double[i][j] = 0;
                ca.ptr.pp_complex[i][j] = ae_complex_from_d(0);
            }
        }
        for(i=0; i<=ae_minint(m, n, _state)-1; i++)
        {
            ra.ptr.pp_double[i][i] = 1+10*ae_randomreal(_state);
            ca.ptr.pp_complex[i][i] = ae_complex_from_d(1+10*ae_randomreal(_state));
        }
        cmatrixrndorthogonalfromtheleft(&ca, m, n, _state);
        cmatrixrndorthogonalfromtheright(&ca, m, n, _state);
        rmatrixrndorthogonalfromtheleft(&ra, m, n, _state);
        rmatrixrndorthogonalfromtheright(&ra, m, n, _state);
        testtrfacunit_testcluproblem(&ca, m, n, threshold, &cerr, &properr, _state);
        testtrfacunit_testrluproblem(&ra, m, n, threshold, &rerr, &properr, _state);
    }
    
    /*
     * Test Cholesky
     */
    for(n=1; n<=maxmn; n++)
    {
        
        /*
         * Load CA (HPD matrix with low condition number),
         *      CAL and CAU - its lower and upper triangles
         */
        hpdmatrixrndcond(n, 1+50*ae_randomreal(_state), &ca, _state);
        ae_matrix_set_length(&cal, n, n, _state);
        ae_matrix_set_length(&cau, n, n, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                cal.ptr.pp_complex[i][j] = ae_complex_from_d(i);
                cau.ptr.pp_complex[i][j] = ae_complex_from_d(j);
            }
        }
        for(i=0; i<=n-1; i++)
        {
            ae_v_cmove(&cal.ptr.pp_complex[i][0], 1, &ca.ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i));
            ae_v_cmove(&cau.ptr.pp_complex[i][i], 1, &ca.ptr.pp_complex[i][i], 1, "N", ae_v_len(i,n-1));
        }
        
        /*
         * Test HPDMatrixCholesky:
         * 1. it must leave upper (lower) part unchanged
         * 2. max(A-L*L^H) must be small
         */
        if( hpdmatrixcholesky(&cal, n, ae_false, _state) )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( j>i )
                    {
                        hpderr = hpderr||ae_c_neq_d(cal.ptr.pp_complex[i][j],i);
                    }
                    else
                    {
                        vc = ae_v_cdotproduct(&cal.ptr.pp_complex[i][0], 1, "N", &cal.ptr.pp_complex[j][0], 1, "Conj", ae_v_len(0,j));
                        hpderr = hpderr||ae_fp_greater(ae_c_abs(ae_c_sub(ca.ptr.pp_complex[i][j],vc), _state),threshold);
                    }
                }
            }
        }
        else
        {
            hpderr = ae_true;
        }
        if( hpdmatrixcholesky(&cau, n, ae_true, _state) )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( j<i )
                    {
                        hpderr = hpderr||ae_c_neq_d(cau.ptr.pp_complex[i][j],j);
                    }
                    else
                    {
                        vc = ae_v_cdotproduct(&cau.ptr.pp_complex[0][i], cau.stride, "Conj", &cau.ptr.pp_complex[0][j], cau.stride, "N", ae_v_len(0,i));
                        hpderr = hpderr||ae_fp_greater(ae_c_abs(ae_c_sub(ca.ptr.pp_complex[i][j],vc), _state),threshold);
                    }
                }
            }
        }
        else
        {
            hpderr = ae_true;
        }
        
        /*
         * Load RA (SPD matrix with low condition number),
         *      RAL and RAU - its lower and upper triangles
         */
        spdmatrixrndcond(n, 1+50*ae_randomreal(_state), &ra, _state);
        ae_matrix_set_length(&ral, n, n, _state);
        ae_matrix_set_length(&rau, n, n, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ral.ptr.pp_double[i][j] = i;
                rau.ptr.pp_double[i][j] = j;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            ae_v_move(&ral.ptr.pp_double[i][0], 1, &ra.ptr.pp_double[i][0], 1, ae_v_len(0,i));
            ae_v_move(&rau.ptr.pp_double[i][i], 1, &ra.ptr.pp_double[i][i], 1, ae_v_len(i,n-1));
        }
        
        /*
         * Test SPDMatrixCholesky:
         * 1. it must leave upper (lower) part unchanged
         * 2. max(A-L*L^H) must be small
         */
        if( spdmatrixcholesky(&ral, n, ae_false, _state) )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( j>i )
                    {
                        spderr = spderr||ae_fp_neq(ral.ptr.pp_double[i][j],i);
                    }
                    else
                    {
                        vr = ae_v_dotproduct(&ral.ptr.pp_double[i][0], 1, &ral.ptr.pp_double[j][0], 1, ae_v_len(0,j));
                        spderr = spderr||ae_fp_greater(ae_fabs(ra.ptr.pp_double[i][j]-vr, _state),threshold);
                    }
                }
            }
        }
        else
        {
            spderr = ae_true;
        }
        if( spdmatrixcholesky(&rau, n, ae_true, _state) )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( j<i )
                    {
                        spderr = spderr||ae_fp_neq(rau.ptr.pp_double[i][j],j);
                    }
                    else
                    {
                        vr = ae_v_dotproduct(&rau.ptr.pp_double[0][i], rau.stride, &rau.ptr.pp_double[0][j], rau.stride, ae_v_len(0,i));
                        spderr = spderr||ae_fp_greater(ae_fabs(ra.ptr.pp_double[i][j]-vr, _state),threshold);
                    }
                }
            }
        }
        else
        {
            spderr = ae_true;
        }
    }
    
    /*
     * report
     */
    waserrors = (((rerr||spderr)||cerr)||hpderr)||properr;
    if( !silent )
    {
        printf("TESTING TRIANGULAR FACTORIZATIONS\n");
        printf("* REAL:                                  ");
        if( rerr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SPD:                                   ");
        if( spderr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COMPLEX:                               ");
        if( cerr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* HPD:                                   ");
        if( hpderr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* OTHER PROPERTIES:                      ");
        if( properr )
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


static void testtrfacunit_testcluproblem(/* Complex */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     double threshold,
     ae_bool* err,
     ae_bool* properr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix ca;
    ae_matrix cl;
    ae_matrix cu;
    ae_matrix ca2;
    ae_vector ct;
    ae_int_t i;
    ae_int_t j;
    ae_int_t minmn;
    ae_complex v;
    ae_vector p;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&ca, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&cl, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&cu, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&ca2, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&ct, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&p, 0, DT_INT, _state, ae_true);

    minmn = ae_minint(m, n, _state);
    
    /*
     * PLU test
     */
    ae_matrix_set_length(&ca, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        ae_v_cmove(&ca.ptr.pp_complex[i][0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,n-1));
    }
    cmatrixplu(&ca, m, n, &p, _state);
    for(i=0; i<=minmn-1; i++)
    {
        if( p.ptr.p_int[i]<i||p.ptr.p_int[i]>=m )
        {
            *properr = ae_false;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(&cl, m, minmn, _state);
    for(j=0; j<=minmn-1; j++)
    {
        for(i=0; i<=j-1; i++)
        {
            cl.ptr.pp_complex[i][j] = ae_complex_from_d(0.0);
        }
        cl.ptr.pp_complex[j][j] = ae_complex_from_d(1.0);
        for(i=j+1; i<=m-1; i++)
        {
            cl.ptr.pp_complex[i][j] = ca.ptr.pp_complex[i][j];
        }
    }
    ae_matrix_set_length(&cu, minmn, n, _state);
    for(i=0; i<=minmn-1; i++)
    {
        for(j=0; j<=i-1; j++)
        {
            cu.ptr.pp_complex[i][j] = ae_complex_from_d(0.0);
        }
        for(j=i; j<=n-1; j++)
        {
            cu.ptr.pp_complex[i][j] = ca.ptr.pp_complex[i][j];
        }
    }
    ae_matrix_set_length(&ca2, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            v = ae_v_cdotproduct(&cl.ptr.pp_complex[i][0], 1, "N", &cu.ptr.pp_complex[0][j], cu.stride, "N", ae_v_len(0,minmn-1));
            ca2.ptr.pp_complex[i][j] = v;
        }
    }
    ae_vector_set_length(&ct, n, _state);
    for(i=minmn-1; i>=0; i--)
    {
        if( i!=p.ptr.p_int[i] )
        {
            ae_v_cmove(&ct.ptr.p_complex[0], 1, &ca2.ptr.pp_complex[i][0], 1, "N", ae_v_len(0,n-1));
            ae_v_cmove(&ca2.ptr.pp_complex[i][0], 1, &ca2.ptr.pp_complex[p.ptr.p_int[i]][0], 1, "N", ae_v_len(0,n-1));
            ae_v_cmove(&ca2.ptr.pp_complex[p.ptr.p_int[i]][0], 1, &ct.ptr.p_complex[0], 1, "N", ae_v_len(0,n-1));
        }
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            *err = *err||ae_fp_greater(ae_c_abs(ae_c_sub(a->ptr.pp_complex[i][j],ca2.ptr.pp_complex[i][j]), _state),threshold);
        }
    }
    
    /*
     * LUP test
     */
    ae_matrix_set_length(&ca, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        ae_v_cmove(&ca.ptr.pp_complex[i][0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,n-1));
    }
    cmatrixlup(&ca, m, n, &p, _state);
    for(i=0; i<=minmn-1; i++)
    {
        if( p.ptr.p_int[i]<i||p.ptr.p_int[i]>=n )
        {
            *properr = ae_false;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(&cl, m, minmn, _state);
    for(j=0; j<=minmn-1; j++)
    {
        for(i=0; i<=j-1; i++)
        {
            cl.ptr.pp_complex[i][j] = ae_complex_from_d(0.0);
        }
        for(i=j; i<=m-1; i++)
        {
            cl.ptr.pp_complex[i][j] = ca.ptr.pp_complex[i][j];
        }
    }
    ae_matrix_set_length(&cu, minmn, n, _state);
    for(i=0; i<=minmn-1; i++)
    {
        for(j=0; j<=i-1; j++)
        {
            cu.ptr.pp_complex[i][j] = ae_complex_from_d(0.0);
        }
        cu.ptr.pp_complex[i][i] = ae_complex_from_d(1.0);
        for(j=i+1; j<=n-1; j++)
        {
            cu.ptr.pp_complex[i][j] = ca.ptr.pp_complex[i][j];
        }
    }
    ae_matrix_set_length(&ca2, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            v = ae_v_cdotproduct(&cl.ptr.pp_complex[i][0], 1, "N", &cu.ptr.pp_complex[0][j], cu.stride, "N", ae_v_len(0,minmn-1));
            ca2.ptr.pp_complex[i][j] = v;
        }
    }
    ae_vector_set_length(&ct, m, _state);
    for(i=minmn-1; i>=0; i--)
    {
        if( i!=p.ptr.p_int[i] )
        {
            ae_v_cmove(&ct.ptr.p_complex[0], 1, &ca2.ptr.pp_complex[0][i], ca2.stride, "N", ae_v_len(0,m-1));
            ae_v_cmove(&ca2.ptr.pp_complex[0][i], ca2.stride, &ca2.ptr.pp_complex[0][p.ptr.p_int[i]], ca2.stride, "N", ae_v_len(0,m-1));
            ae_v_cmove(&ca2.ptr.pp_complex[0][p.ptr.p_int[i]], ca2.stride, &ct.ptr.p_complex[0], 1, "N", ae_v_len(0,m-1));
        }
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            *err = *err||ae_fp_greater(ae_c_abs(ae_c_sub(a->ptr.pp_complex[i][j],ca2.ptr.pp_complex[i][j]), _state),threshold);
        }
    }
    ae_frame_leave(_state);
}


static void testtrfacunit_testrluproblem(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     double threshold,
     ae_bool* err,
     ae_bool* properr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix ca;
    ae_matrix cl;
    ae_matrix cu;
    ae_matrix ca2;
    ae_vector ct;
    ae_int_t i;
    ae_int_t j;
    ae_int_t minmn;
    double v;
    ae_vector p;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&ca, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cl, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cu, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ca2, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&p, 0, DT_INT, _state, ae_true);

    minmn = ae_minint(m, n, _state);
    
    /*
     * PLU test
     */
    ae_matrix_set_length(&ca, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        ae_v_move(&ca.ptr.pp_double[i][0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
    }
    rmatrixplu(&ca, m, n, &p, _state);
    for(i=0; i<=minmn-1; i++)
    {
        if( p.ptr.p_int[i]<i||p.ptr.p_int[i]>=m )
        {
            *properr = ae_false;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(&cl, m, minmn, _state);
    for(j=0; j<=minmn-1; j++)
    {
        for(i=0; i<=j-1; i++)
        {
            cl.ptr.pp_double[i][j] = 0.0;
        }
        cl.ptr.pp_double[j][j] = 1.0;
        for(i=j+1; i<=m-1; i++)
        {
            cl.ptr.pp_double[i][j] = ca.ptr.pp_double[i][j];
        }
    }
    ae_matrix_set_length(&cu, minmn, n, _state);
    for(i=0; i<=minmn-1; i++)
    {
        for(j=0; j<=i-1; j++)
        {
            cu.ptr.pp_double[i][j] = 0.0;
        }
        for(j=i; j<=n-1; j++)
        {
            cu.ptr.pp_double[i][j] = ca.ptr.pp_double[i][j];
        }
    }
    ae_matrix_set_length(&ca2, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            v = ae_v_dotproduct(&cl.ptr.pp_double[i][0], 1, &cu.ptr.pp_double[0][j], cu.stride, ae_v_len(0,minmn-1));
            ca2.ptr.pp_double[i][j] = v;
        }
    }
    ae_vector_set_length(&ct, n, _state);
    for(i=minmn-1; i>=0; i--)
    {
        if( i!=p.ptr.p_int[i] )
        {
            ae_v_move(&ct.ptr.p_double[0], 1, &ca2.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
            ae_v_move(&ca2.ptr.pp_double[i][0], 1, &ca2.ptr.pp_double[p.ptr.p_int[i]][0], 1, ae_v_len(0,n-1));
            ae_v_move(&ca2.ptr.pp_double[p.ptr.p_int[i]][0], 1, &ct.ptr.p_double[0], 1, ae_v_len(0,n-1));
        }
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            *err = *err||ae_fp_greater(ae_fabs(a->ptr.pp_double[i][j]-ca2.ptr.pp_double[i][j], _state),threshold);
        }
    }
    
    /*
     * LUP test
     */
    ae_matrix_set_length(&ca, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        ae_v_move(&ca.ptr.pp_double[i][0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
    }
    rmatrixlup(&ca, m, n, &p, _state);
    for(i=0; i<=minmn-1; i++)
    {
        if( p.ptr.p_int[i]<i||p.ptr.p_int[i]>=n )
        {
            *properr = ae_false;
            ae_frame_leave(_state);
            return;
        }
    }
    ae_matrix_set_length(&cl, m, minmn, _state);
    for(j=0; j<=minmn-1; j++)
    {
        for(i=0; i<=j-1; i++)
        {
            cl.ptr.pp_double[i][j] = 0.0;
        }
        for(i=j; i<=m-1; i++)
        {
            cl.ptr.pp_double[i][j] = ca.ptr.pp_double[i][j];
        }
    }
    ae_matrix_set_length(&cu, minmn, n, _state);
    for(i=0; i<=minmn-1; i++)
    {
        for(j=0; j<=i-1; j++)
        {
            cu.ptr.pp_double[i][j] = 0.0;
        }
        cu.ptr.pp_double[i][i] = 1.0;
        for(j=i+1; j<=n-1; j++)
        {
            cu.ptr.pp_double[i][j] = ca.ptr.pp_double[i][j];
        }
    }
    ae_matrix_set_length(&ca2, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            v = ae_v_dotproduct(&cl.ptr.pp_double[i][0], 1, &cu.ptr.pp_double[0][j], cu.stride, ae_v_len(0,minmn-1));
            ca2.ptr.pp_double[i][j] = v;
        }
    }
    ae_vector_set_length(&ct, m, _state);
    for(i=minmn-1; i>=0; i--)
    {
        if( i!=p.ptr.p_int[i] )
        {
            ae_v_move(&ct.ptr.p_double[0], 1, &ca2.ptr.pp_double[0][i], ca2.stride, ae_v_len(0,m-1));
            ae_v_move(&ca2.ptr.pp_double[0][i], ca2.stride, &ca2.ptr.pp_double[0][p.ptr.p_int[i]], ca2.stride, ae_v_len(0,m-1));
            ae_v_move(&ca2.ptr.pp_double[0][p.ptr.p_int[i]], ca2.stride, &ct.ptr.p_double[0], 1, ae_v_len(0,m-1));
        }
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            *err = *err||ae_fp_greater(ae_fabs(a->ptr.pp_double[i][j]-ca2.ptr.pp_double[i][j], _state),threshold);
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
