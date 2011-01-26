

#include <stdafx.h>
#include <stdio.h>
#include "testtsortunit.h"


/*$ Declarations $*/
static void testtsortunit_unset2d(/* Complex */ ae_matrix* a,
     ae_state *_state);
static void testtsortunit_unset1d(/* Real    */ ae_vector* a,
     ae_state *_state);
static void testtsortunit_unset1di(/* Integer */ ae_vector* a,
     ae_state *_state);
static void testtsortunit_testsortresults(/* Real    */ ae_vector* asorted,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     /* Real    */ ae_vector* aoriginal,
     ae_int_t n,
     ae_bool* waserrors,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Testing tag sort
*************************************************************************/
ae_bool testtsort(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_int_t n;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t passcount;
    ae_int_t maxn;
    ae_vector a;
    ae_vector a0;
    ae_vector a1;
    ae_vector a2;
    ae_vector a3;
    ae_vector ar;
    ae_vector ai;
    ae_vector p1;
    ae_vector p2;
    ae_vector bufr1;
    ae_vector bufr2;
    ae_vector bufi1;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&a, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&a0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&a1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&a2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&a3, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ar, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ai, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p1, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p2, 0, DT_INT, _state, ae_true);
    ae_vector_init(&bufr1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bufr2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bufi1, 0, DT_INT, _state, ae_true);

    waserrors = ae_false;
    maxn = 100;
    passcount = 10;
    
    /*
     * Test tagsort
     */
    for(n=1; n<=maxn; n++)
    {
        for(pass=1; pass<=passcount; pass++)
        {
            
            /*
             * (probably) distinct sort:
             * * first sort A0 using TagSort and test sort results
             * * now we can use A0 as reference point and test other functions
             */
            testtsortunit_unset1di(&p1, _state);
            testtsortunit_unset1di(&p2, _state);
            ae_vector_set_length(&a, n, _state);
            ae_vector_set_length(&a0, n, _state);
            ae_vector_set_length(&a1, n, _state);
            ae_vector_set_length(&a2, n, _state);
            ae_vector_set_length(&a3, n, _state);
            ae_vector_set_length(&ar, n, _state);
            ae_vector_set_length(&ai, n, _state);
            for(i=0; i<=n-1; i++)
            {
                a.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                a0.ptr.p_double[i] = a.ptr.p_double[i];
                a1.ptr.p_double[i] = a.ptr.p_double[i];
                a2.ptr.p_double[i] = a.ptr.p_double[i];
                a3.ptr.p_double[i] = a.ptr.p_double[i];
                ar.ptr.p_double[i] = i;
                ai.ptr.p_int[i] = i;
            }
            tagsort(&a0, n, &p1, &p2, _state);
            testtsortunit_testsortresults(&a0, &p1, &p2, &a, n, &waserrors, _state);
            tagsortfasti(&a1, &ai, &bufr1, &bufi1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a1.ptr.p_double[i],a0.ptr.p_double[i]))||ai.ptr.p_int[i]!=p1.ptr.p_int[i];
            }
            tagsortfastr(&a2, &ar, &bufr1, &bufr2, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a2.ptr.p_double[i],a0.ptr.p_double[i]))||ae_fp_neq(ar.ptr.p_double[i],p1.ptr.p_int[i]);
            }
            tagsortfast(&a3, &bufr1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = waserrors||ae_fp_neq(a3.ptr.p_double[i],a0.ptr.p_double[i]);
            }
            
            /*
             * non-distinct sort
             */
            testtsortunit_unset1di(&p1, _state);
            testtsortunit_unset1di(&p2, _state);
            ae_vector_set_length(&a, n, _state);
            ae_vector_set_length(&a0, n, _state);
            ae_vector_set_length(&a1, n, _state);
            ae_vector_set_length(&a2, n, _state);
            ae_vector_set_length(&a3, n, _state);
            ae_vector_set_length(&ar, n, _state);
            ae_vector_set_length(&ai, n, _state);
            for(i=0; i<=n-1; i++)
            {
                a.ptr.p_double[i] = i/2;
                a0.ptr.p_double[i] = a.ptr.p_double[i];
                a1.ptr.p_double[i] = a.ptr.p_double[i];
                a2.ptr.p_double[i] = a.ptr.p_double[i];
                a3.ptr.p_double[i] = a.ptr.p_double[i];
                ar.ptr.p_double[i] = i;
                ai.ptr.p_int[i] = i;
            }
            tagsort(&a0, n, &p1, &p2, _state);
            testtsortunit_testsortresults(&a0, &p1, &p2, &a, n, &waserrors, _state);
            tagsortfasti(&a1, &ai, &bufr1, &bufi1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a1.ptr.p_double[i],a0.ptr.p_double[i]))||ai.ptr.p_int[i]!=p1.ptr.p_int[i];
            }
            tagsortfastr(&a2, &ar, &bufr1, &bufr2, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a2.ptr.p_double[i],a0.ptr.p_double[i]))||ae_fp_neq(ar.ptr.p_double[i],p1.ptr.p_int[i]);
            }
            tagsortfast(&a3, &bufr1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = waserrors||ae_fp_neq(a3.ptr.p_double[i],a0.ptr.p_double[i]);
            }
            
            /*
             * 'All same' sort
             */
            testtsortunit_unset1di(&p1, _state);
            testtsortunit_unset1di(&p2, _state);
            ae_vector_set_length(&a, n, _state);
            ae_vector_set_length(&a0, n, _state);
            ae_vector_set_length(&a1, n, _state);
            ae_vector_set_length(&a2, n, _state);
            ae_vector_set_length(&a3, n, _state);
            ae_vector_set_length(&ar, n, _state);
            ae_vector_set_length(&ai, n, _state);
            for(i=0; i<=n-1; i++)
            {
                a.ptr.p_double[i] = 0;
                a0.ptr.p_double[i] = a.ptr.p_double[i];
                a1.ptr.p_double[i] = a.ptr.p_double[i];
                a2.ptr.p_double[i] = a.ptr.p_double[i];
                a3.ptr.p_double[i] = a.ptr.p_double[i];
                ar.ptr.p_double[i] = i;
                ai.ptr.p_int[i] = i;
            }
            tagsort(&a0, n, &p1, &p2, _state);
            testtsortunit_testsortresults(&a0, &p1, &p2, &a, n, &waserrors, _state);
            tagsortfasti(&a1, &ai, &bufr1, &bufi1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a1.ptr.p_double[i],a0.ptr.p_double[i]))||ai.ptr.p_int[i]!=p1.ptr.p_int[i];
            }
            tagsortfastr(&a2, &ar, &bufr1, &bufr2, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a2.ptr.p_double[i],a0.ptr.p_double[i]))||ae_fp_neq(ar.ptr.p_double[i],p1.ptr.p_int[i]);
            }
            tagsortfast(&a3, &bufr1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = waserrors||ae_fp_neq(a3.ptr.p_double[i],a0.ptr.p_double[i]);
            }
            
            /*
             * 0-1 sort
             */
            testtsortunit_unset1di(&p1, _state);
            testtsortunit_unset1di(&p2, _state);
            ae_vector_set_length(&a, n, _state);
            ae_vector_set_length(&a0, n, _state);
            ae_vector_set_length(&a1, n, _state);
            ae_vector_set_length(&a2, n, _state);
            ae_vector_set_length(&a3, n, _state);
            ae_vector_set_length(&ar, n, _state);
            ae_vector_set_length(&ai, n, _state);
            for(i=0; i<=n-1; i++)
            {
                a.ptr.p_double[i] = ae_randominteger(2, _state);
                a0.ptr.p_double[i] = a.ptr.p_double[i];
                a1.ptr.p_double[i] = a.ptr.p_double[i];
                a2.ptr.p_double[i] = a.ptr.p_double[i];
                a3.ptr.p_double[i] = a.ptr.p_double[i];
                ar.ptr.p_double[i] = i;
                ai.ptr.p_int[i] = i;
            }
            tagsort(&a0, n, &p1, &p2, _state);
            testtsortunit_testsortresults(&a0, &p1, &p2, &a, n, &waserrors, _state);
            tagsortfasti(&a1, &ai, &bufr1, &bufi1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a1.ptr.p_double[i],a0.ptr.p_double[i]))||ai.ptr.p_int[i]!=p1.ptr.p_int[i];
            }
            tagsortfastr(&a2, &ar, &bufr1, &bufr2, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = (waserrors||ae_fp_neq(a2.ptr.p_double[i],a0.ptr.p_double[i]))||ae_fp_neq(ar.ptr.p_double[i],p1.ptr.p_int[i]);
            }
            tagsortfast(&a3, &bufr1, n, _state);
            for(i=0; i<=n-1; i++)
            {
                waserrors = waserrors||ae_fp_neq(a3.ptr.p_double[i],a0.ptr.p_double[i]);
            }
        }
    }
    
    /*
     * report
     */
    if( !silent )
    {
        printf("TESTING TAGSORT\n");
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
Unsets 2D array.
*************************************************************************/
static void testtsortunit_unset2d(/* Complex */ ae_matrix* a,
     ae_state *_state)
{


    ae_matrix_set_length(a, 0+1, 0+1, _state);
    a->ptr.pp_complex[0][0] = ae_complex_from_d(2*ae_randomreal(_state)-1);
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void testtsortunit_unset1d(/* Real    */ ae_vector* a,
     ae_state *_state)
{


    ae_vector_set_length(a, 0+1, _state);
    a->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void testtsortunit_unset1di(/* Integer */ ae_vector* a,
     ae_state *_state)
{


    ae_vector_set_length(a, 0+1, _state);
    a->ptr.p_int[0] = ae_randominteger(3, _state)-1;
}


static void testtsortunit_testsortresults(/* Real    */ ae_vector* asorted,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     /* Real    */ ae_vector* aoriginal,
     ae_int_t n,
     ae_bool* waserrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector a2;
    double t;
    ae_vector f;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&a2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&f, 0, DT_INT, _state, ae_true);

    ae_vector_set_length(&a2, n-1+1, _state);
    ae_vector_set_length(&f, n-1+1, _state);
    
    /*
     * is set ordered?
     */
    for(i=0; i<=n-2; i++)
    {
        *waserrors = *waserrors||ae_fp_greater(asorted->ptr.p_double[i],asorted->ptr.p_double[i+1]);
    }
    
    /*
     * P1 correctness
     */
    for(i=0; i<=n-1; i++)
    {
        *waserrors = *waserrors||ae_fp_neq(asorted->ptr.p_double[i],aoriginal->ptr.p_double[p1->ptr.p_int[i]]);
    }
    for(i=0; i<=n-1; i++)
    {
        f.ptr.p_int[i] = 0;
    }
    for(i=0; i<=n-1; i++)
    {
        f.ptr.p_int[p1->ptr.p_int[i]] = f.ptr.p_int[p1->ptr.p_int[i]]+1;
    }
    for(i=0; i<=n-1; i++)
    {
        *waserrors = *waserrors||f.ptr.p_int[i]!=1;
    }
    
    /*
     * P2 correctness
     */
    for(i=0; i<=n-1; i++)
    {
        a2.ptr.p_double[i] = aoriginal->ptr.p_double[i];
    }
    for(i=0; i<=n-1; i++)
    {
        if( p2->ptr.p_int[i]!=i )
        {
            t = a2.ptr.p_double[i];
            a2.ptr.p_double[i] = a2.ptr.p_double[p2->ptr.p_int[i]];
            a2.ptr.p_double[p2->ptr.p_int[i]] = t;
        }
    }
    for(i=0; i<=n-1; i++)
    {
        *waserrors = *waserrors||ae_fp_neq(asorted->ptr.p_double[i],a2.ptr.p_double[i]);
    }
    ae_frame_leave(_state);
}


/*$ End $*/
