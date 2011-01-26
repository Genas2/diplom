

#include <stdafx.h>
#include <stdio.h>
#include "testalglibbasicsunit.h"


/*$ Declarations $*/
static ae_bool testalglibbasicsunit_testcomplexarithmetics(ae_bool silent,
     ae_state *_state);
static ae_bool testalglibbasicsunit_testieeespecial(ae_bool silent,
     ae_state *_state);
static ae_bool testalglibbasicsunit_testswapfunctions(ae_bool silent,
     ae_state *_state);


/*$ Body $*/


ae_bool testalglibbasics(ae_bool silent, ae_state *_state)
{
    ae_bool result;


    result = ae_true;
    result = result&&testalglibbasicsunit_testcomplexarithmetics(silent, _state);
    result = result&&testalglibbasicsunit_testieeespecial(silent, _state);
    result = result&&testalglibbasicsunit_testswapfunctions(silent, _state);
    if( !silent )
    {
        printf("\n\n");
    }
    return result;
}


/*************************************************************************
Complex arithmetics test
*************************************************************************/
static ae_bool testalglibbasicsunit_testcomplexarithmetics(ae_bool silent,
     ae_state *_state)
{
    ae_bool absc;
    ae_bool addcc;
    ae_bool addcr;
    ae_bool addrc;
    ae_bool subcc;
    ae_bool subcr;
    ae_bool subrc;
    ae_bool mulcc;
    ae_bool mulcr;
    ae_bool mulrc;
    ae_bool divcc;
    ae_bool divcr;
    ae_bool divrc;
    ae_complex ca;
    ae_complex cb;
    ae_complex res;
    double ra;
    double rb;
    double threshold;
    ae_int_t pass;
    ae_int_t passcount;
    ae_bool result;


    threshold = 100*ae_machineepsilon;
    passcount = 1000;
    result = ae_true;
    absc = ae_true;
    addcc = ae_true;
    addcr = ae_true;
    addrc = ae_true;
    subcc = ae_true;
    subcr = ae_true;
    subrc = ae_true;
    mulcc = ae_true;
    mulcr = ae_true;
    mulrc = ae_true;
    divcc = ae_true;
    divcr = ae_true;
    divrc = ae_true;
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Test AbsC
         */
        ca.x = 2*ae_randomreal(_state)-1;
        ca.y = 2*ae_randomreal(_state)-1;
        ra = ae_c_abs(ca, _state);
        absc = absc&&ae_fp_less(ae_fabs(ra-ae_sqrt(ae_sqr(ca.x, _state)+ae_sqr(ca.y, _state), _state), _state),threshold);
        
        /*
         * test Add
         */
        ca.x = 2*ae_randomreal(_state)-1;
        ca.y = 2*ae_randomreal(_state)-1;
        cb.x = 2*ae_randomreal(_state)-1;
        cb.y = 2*ae_randomreal(_state)-1;
        ra = 2*ae_randomreal(_state)-1;
        rb = 2*ae_randomreal(_state)-1;
        res = ae_c_add(ca,cb);
        addcc = (addcc&&ae_fp_less(ae_fabs(res.x-ca.x-cb.x, _state),threshold))&&ae_fp_less(ae_fabs(res.y-ca.y-cb.y, _state),threshold);
        res = ae_c_add_d(ca,rb);
        addcr = (addcr&&ae_fp_less(ae_fabs(res.x-ca.x-rb, _state),threshold))&&ae_fp_less(ae_fabs(res.y-ca.y, _state),threshold);
        res = ae_c_add_d(cb,ra);
        addrc = (addrc&&ae_fp_less(ae_fabs(res.x-ra-cb.x, _state),threshold))&&ae_fp_less(ae_fabs(res.y-cb.y, _state),threshold);
        
        /*
         * test Sub
         */
        ca.x = 2*ae_randomreal(_state)-1;
        ca.y = 2*ae_randomreal(_state)-1;
        cb.x = 2*ae_randomreal(_state)-1;
        cb.y = 2*ae_randomreal(_state)-1;
        ra = 2*ae_randomreal(_state)-1;
        rb = 2*ae_randomreal(_state)-1;
        res = ae_c_sub(ca,cb);
        subcc = (subcc&&ae_fp_less(ae_fabs(res.x-(ca.x-cb.x), _state),threshold))&&ae_fp_less(ae_fabs(res.y-(ca.y-cb.y), _state),threshold);
        res = ae_c_sub_d(ca,rb);
        subcr = (subcr&&ae_fp_less(ae_fabs(res.x-(ca.x-rb), _state),threshold))&&ae_fp_less(ae_fabs(res.y-ca.y, _state),threshold);
        res = ae_c_d_sub(ra,cb);
        subrc = (subrc&&ae_fp_less(ae_fabs(res.x-(ra-cb.x), _state),threshold))&&ae_fp_less(ae_fabs(res.y+cb.y, _state),threshold);
        
        /*
         * test Mul
         */
        ca.x = 2*ae_randomreal(_state)-1;
        ca.y = 2*ae_randomreal(_state)-1;
        cb.x = 2*ae_randomreal(_state)-1;
        cb.y = 2*ae_randomreal(_state)-1;
        ra = 2*ae_randomreal(_state)-1;
        rb = 2*ae_randomreal(_state)-1;
        res = ae_c_mul(ca,cb);
        mulcc = (mulcc&&ae_fp_less(ae_fabs(res.x-(ca.x*cb.x-ca.y*cb.y), _state),threshold))&&ae_fp_less(ae_fabs(res.y-(ca.x*cb.y+ca.y*cb.x), _state),threshold);
        res = ae_c_mul_d(ca,rb);
        mulcr = (mulcr&&ae_fp_less(ae_fabs(res.x-ca.x*rb, _state),threshold))&&ae_fp_less(ae_fabs(res.y-ca.y*rb, _state),threshold);
        res = ae_c_mul_d(cb,ra);
        mulrc = (mulrc&&ae_fp_less(ae_fabs(res.x-ra*cb.x, _state),threshold))&&ae_fp_less(ae_fabs(res.y-ra*cb.y, _state),threshold);
        
        /*
         * test Div
         */
        ca.x = 2*ae_randomreal(_state)-1;
        ca.y = 2*ae_randomreal(_state)-1;
        do
        {
            cb.x = 2*ae_randomreal(_state)-1;
            cb.y = 2*ae_randomreal(_state)-1;
        }
        while(ae_fp_less_eq(ae_c_abs(cb, _state),0.5));
        ra = 2*ae_randomreal(_state)-1;
        do
        {
            rb = 2*ae_randomreal(_state)-1;
        }
        while(ae_fp_less_eq(ae_fabs(rb, _state),0.5));
        res = ae_c_div(ca,cb);
        divcc = (divcc&&ae_fp_less(ae_fabs(ae_c_mul(res,cb).x-ca.x, _state),threshold))&&ae_fp_less(ae_fabs(ae_c_mul(res,cb).y-ca.y, _state),threshold);
        res = ae_c_div_d(ca,rb);
        divcr = (divcr&&ae_fp_less(ae_fabs(res.x-ca.x/rb, _state),threshold))&&ae_fp_less(ae_fabs(res.y-ca.y/rb, _state),threshold);
        res = ae_c_d_div(ra,cb);
        divrc = (divrc&&ae_fp_less(ae_fabs(ae_c_mul(res,cb).x-ra, _state),threshold))&&ae_fp_less(ae_fabs(ae_c_mul(res,cb).y, _state),threshold);
    }
    
    /*
     * summary
     */
    result = result&&absc;
    result = result&&addcc;
    result = result&&addcr;
    result = result&&addrc;
    result = result&&subcc;
    result = result&&subcr;
    result = result&&subrc;
    result = result&&mulcc;
    result = result&&mulcr;
    result = result&&mulrc;
    result = result&&divcc;
    result = result&&divcr;
    result = result&&divrc;
    if( !silent )
    {
        if( result )
        {
            printf("COMPLEX ARITHMETICS:                     OK\n");
        }
        else
        {
            printf("COMPLEX ARITHMETICS:                     FAILED\n");
            printf("* AddCC - - - - - - - - - - - - - - - -  ");
            if( addcc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* AddCR - - - - - - - - - - - - - - - -  ");
            if( addcr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* AddRC - - - - - - - - - - - - - - - -  ");
            if( addrc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* SubCC - - - - - - - - - - - - - - - -  ");
            if( subcc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* SubCR - - - - - - - - - - - - - - - -  ");
            if( subcr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* SubRC - - - - - - - - - - - - - - - -  ");
            if( subrc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* MulCC - - - - - - - - - - - - - - - -  ");
            if( mulcc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* MulCR - - - - - - - - - - - - - - - -  ");
            if( mulcr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* MulRC - - - - - - - - - - - - - - - -  ");
            if( mulrc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* DivCC - - - - - - - - - - - - - - - -  ");
            if( divcc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* DivCR - - - - - - - - - - - - - - - -  ");
            if( divcr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* DivRC - - - - - - - - - - - - - - - -  ");
            if( divrc )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
        }
    }
    return result;
}


/*************************************************************************
Tests for IEEE special quantities
*************************************************************************/
static ae_bool testalglibbasicsunit_testieeespecial(ae_bool silent,
     ae_state *_state)
{
    ae_bool oknan;
    ae_bool okinf;
    ae_bool okother;
    double v1;
    double v2;
    double v3;
    ae_bool result;


    result = ae_true;
    oknan = ae_true;
    okinf = ae_true;
    okother = ae_true;
    
    /*
     * Test classification functions
     */
    okother = okother&&!ae_isinf(_state->v_nan, _state);
    okother = okother&&ae_isinf(_state->v_posinf, _state);
    okother = okother&&!ae_isinf(ae_maxrealnumber, _state);
    okother = okother&&!ae_isinf(1.0, _state);
    okother = okother&&!ae_isinf(ae_minrealnumber, _state);
    okother = okother&&!ae_isinf(0.0, _state);
    okother = okother&&!ae_isinf(-ae_minrealnumber, _state);
    okother = okother&&!ae_isinf(-1.0, _state);
    okother = okother&&!ae_isinf(-ae_maxrealnumber, _state);
    okother = okother&&ae_isinf(_state->v_neginf, _state);
    okother = okother&&!ae_isposinf(_state->v_nan, _state);
    okother = okother&&ae_isposinf(_state->v_posinf, _state);
    okother = okother&&!ae_isposinf(ae_maxrealnumber, _state);
    okother = okother&&!ae_isposinf(1.0, _state);
    okother = okother&&!ae_isposinf(ae_minrealnumber, _state);
    okother = okother&&!ae_isposinf(0.0, _state);
    okother = okother&&!ae_isposinf(-ae_minrealnumber, _state);
    okother = okother&&!ae_isposinf(-1.0, _state);
    okother = okother&&!ae_isposinf(-ae_maxrealnumber, _state);
    okother = okother&&!ae_isposinf(_state->v_neginf, _state);
    okother = okother&&!ae_isneginf(_state->v_nan, _state);
    okother = okother&&!ae_isneginf(_state->v_posinf, _state);
    okother = okother&&!ae_isneginf(ae_maxrealnumber, _state);
    okother = okother&&!ae_isneginf(1.0, _state);
    okother = okother&&!ae_isneginf(ae_minrealnumber, _state);
    okother = okother&&!ae_isneginf(0.0, _state);
    okother = okother&&!ae_isneginf(-ae_minrealnumber, _state);
    okother = okother&&!ae_isneginf(-1.0, _state);
    okother = okother&&!ae_isneginf(-ae_maxrealnumber, _state);
    okother = okother&&ae_isneginf(_state->v_neginf, _state);
    okother = okother&&ae_isnan(_state->v_nan, _state);
    okother = okother&&!ae_isnan(_state->v_posinf, _state);
    okother = okother&&!ae_isnan(ae_maxrealnumber, _state);
    okother = okother&&!ae_isnan(1.0, _state);
    okother = okother&&!ae_isnan(ae_minrealnumber, _state);
    okother = okother&&!ae_isnan(0.0, _state);
    okother = okother&&!ae_isnan(-ae_minrealnumber, _state);
    okother = okother&&!ae_isnan(-1.0, _state);
    okother = okother&&!ae_isnan(-ae_maxrealnumber, _state);
    okother = okother&&!ae_isnan(_state->v_neginf, _state);
    okother = okother&&!ae_isfinite(_state->v_nan, _state);
    okother = okother&&!ae_isfinite(_state->v_posinf, _state);
    okother = okother&&ae_isfinite(ae_maxrealnumber, _state);
    okother = okother&&ae_isfinite(1.0, _state);
    okother = okother&&ae_isfinite(ae_minrealnumber, _state);
    okother = okother&&ae_isfinite(0.0, _state);
    okother = okother&&ae_isfinite(-ae_minrealnumber, _state);
    okother = okother&&ae_isfinite(-1.0, _state);
    okother = okother&&ae_isfinite(-ae_maxrealnumber, _state);
    okother = okother&&!ae_isfinite(_state->v_neginf, _state);
    
    /*
     * Test NAN
     */
    v1 = _state->v_nan;
    v2 = _state->v_nan;
    oknan = oknan&&ae_isnan(v1, _state);
    oknan = oknan&&ae_fp_neq(v1,v2);
    oknan = oknan&&!ae_fp_eq(v1,v2);
    
    /*
     * Test INF:
     * * basic properties
     * * comparisons involving PosINF on one of the sides
     * * comparisons involving NegINF on one of the sides
     */
    v1 = _state->v_posinf;
    v2 = _state->v_neginf;
    okinf = okinf&&ae_isinf(_state->v_posinf, _state);
    okinf = okinf&&ae_isinf(v1, _state);
    okinf = okinf&&ae_isinf(_state->v_neginf, _state);
    okinf = okinf&&ae_isinf(v2, _state);
    okinf = okinf&&ae_isposinf(_state->v_posinf, _state);
    okinf = okinf&&ae_isposinf(v1, _state);
    okinf = okinf&&!ae_isposinf(_state->v_neginf, _state);
    okinf = okinf&&!ae_isposinf(v2, _state);
    okinf = okinf&&!ae_isneginf(_state->v_posinf, _state);
    okinf = okinf&&!ae_isneginf(v1, _state);
    okinf = okinf&&ae_isneginf(_state->v_neginf, _state);
    okinf = okinf&&ae_isneginf(v2, _state);
    okinf = okinf&&ae_fp_eq(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&ae_fp_eq(_state->v_posinf,v1);
    okinf = okinf&&!ae_fp_eq(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&!ae_fp_eq(_state->v_posinf,v2);
    okinf = okinf&&!ae_fp_eq(_state->v_posinf,0);
    okinf = okinf&&!ae_fp_eq(_state->v_posinf,1.2);
    okinf = okinf&&!ae_fp_eq(_state->v_posinf,-1.2);
    okinf = okinf&&ae_fp_eq(v1,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(v2,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(0,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(-1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_neq(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&!ae_fp_neq(_state->v_posinf,v1);
    okinf = okinf&&ae_fp_neq(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&ae_fp_neq(_state->v_posinf,v2);
    okinf = okinf&&ae_fp_neq(_state->v_posinf,0);
    okinf = okinf&&ae_fp_neq(_state->v_posinf,1.2);
    okinf = okinf&&ae_fp_neq(_state->v_posinf,-1.2);
    okinf = okinf&&!ae_fp_neq(v1,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(v2,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(0,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(1.2,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(-1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,v1);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,v2);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,0);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,1.2);
    okinf = okinf&&!ae_fp_less(_state->v_posinf,-1.2);
    okinf = okinf&&!ae_fp_less(v1,_state->v_posinf);
    okinf = okinf&&ae_fp_less(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_less(v2,_state->v_posinf);
    okinf = okinf&&ae_fp_less(0,_state->v_posinf);
    okinf = okinf&&ae_fp_less(1.2,_state->v_posinf);
    okinf = okinf&&ae_fp_less(-1.2,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(_state->v_posinf,v1);
    okinf = okinf&&!ae_fp_less_eq(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&!ae_fp_less_eq(_state->v_posinf,v2);
    okinf = okinf&&!ae_fp_less_eq(_state->v_posinf,0);
    okinf = okinf&&!ae_fp_less_eq(_state->v_posinf,1.2);
    okinf = okinf&&!ae_fp_less_eq(_state->v_posinf,-1.2);
    okinf = okinf&&ae_fp_less_eq(v1,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(v2,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(0,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(1.2,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(-1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(_state->v_posinf,v1);
    okinf = okinf&&ae_fp_greater(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&ae_fp_greater(_state->v_posinf,v2);
    okinf = okinf&&ae_fp_greater(_state->v_posinf,0);
    okinf = okinf&&ae_fp_greater(_state->v_posinf,1.2);
    okinf = okinf&&ae_fp_greater(_state->v_posinf,-1.2);
    okinf = okinf&&!ae_fp_greater(v1,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(v2,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(0,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(-1.2,_state->v_posinf);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,_state->v_posinf);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,v1);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,v2);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,0);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,1.2);
    okinf = okinf&&ae_fp_greater_eq(_state->v_posinf,-1.2);
    okinf = okinf&&ae_fp_greater_eq(v1,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(v2,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(0,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(-1.2,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,v1);
    okinf = okinf&&ae_fp_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_eq(_state->v_neginf,v2);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,0);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,1.2);
    okinf = okinf&&!ae_fp_eq(_state->v_neginf,-1.2);
    okinf = okinf&&!ae_fp_eq(v1,_state->v_neginf);
    okinf = okinf&&ae_fp_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_eq(v2,_state->v_neginf);
    okinf = okinf&&!ae_fp_eq(0,_state->v_neginf);
    okinf = okinf&&!ae_fp_eq(1.2,_state->v_neginf);
    okinf = okinf&&!ae_fp_eq(-1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,v1);
    okinf = okinf&&!ae_fp_neq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_neq(_state->v_neginf,v2);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,0);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,1.2);
    okinf = okinf&&ae_fp_neq(_state->v_neginf,-1.2);
    okinf = okinf&&ae_fp_neq(v1,_state->v_neginf);
    okinf = okinf&&!ae_fp_neq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_neq(v2,_state->v_neginf);
    okinf = okinf&&ae_fp_neq(0,_state->v_neginf);
    okinf = okinf&&ae_fp_neq(1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_neq(-1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_less(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_less(_state->v_neginf,v1);
    okinf = okinf&&!ae_fp_less(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(_state->v_neginf,v2);
    okinf = okinf&&ae_fp_less(_state->v_neginf,0);
    okinf = okinf&&ae_fp_less(_state->v_neginf,1.2);
    okinf = okinf&&ae_fp_less(_state->v_neginf,-1.2);
    okinf = okinf&&!ae_fp_less(v1,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(v2,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(0,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(1.2,_state->v_neginf);
    okinf = okinf&&!ae_fp_less(-1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,v1);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,v2);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,0);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,1.2);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,-1.2);
    okinf = okinf&&!ae_fp_less_eq(v1,_state->v_neginf);
    okinf = okinf&&ae_fp_less_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_less_eq(v2,_state->v_neginf);
    okinf = okinf&&!ae_fp_less_eq(0,_state->v_neginf);
    okinf = okinf&&!ae_fp_less_eq(1.2,_state->v_neginf);
    okinf = okinf&&!ae_fp_less_eq(-1.2,_state->v_neginf);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,v1);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,v2);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,0);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,1.2);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,-1.2);
    okinf = okinf&&ae_fp_greater(v1,_state->v_neginf);
    okinf = okinf&&!ae_fp_greater(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&!ae_fp_greater(v2,_state->v_neginf);
    okinf = okinf&&ae_fp_greater(0,_state->v_neginf);
    okinf = okinf&&ae_fp_greater(1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_greater(-1.2,_state->v_neginf);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,_state->v_posinf);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,v1);
    okinf = okinf&&ae_fp_greater_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(_state->v_neginf,v2);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,0);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,1.2);
    okinf = okinf&&!ae_fp_greater_eq(_state->v_neginf,-1.2);
    okinf = okinf&&ae_fp_greater_eq(v1,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(_state->v_neginf,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(v2,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(0,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(1.2,_state->v_neginf);
    okinf = okinf&&ae_fp_greater_eq(-1.2,_state->v_neginf);
    
    /*
     * summary
     */
    result = result&&oknan;
    result = result&&okinf;
    result = result&&okother;
    if( !silent )
    {
        if( result )
        {
            printf("IEEE SPECIAL VALUES:                     OK\n");
        }
        else
        {
            printf("IEEE SPECIAL VALUES:                     FAILED\n");
            printf("* NAN - - - - - - - - - - - - - - - - -  ");
            if( oknan )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* INF - - - - - - - - - - - - - - - - -  ");
            if( okinf )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("* FUNCTIONS - - - - - - - - - - - - - -  ");
            if( okother )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
        }
    }
    return result;
}


/*************************************************************************
Tests for swapping functions
*************************************************************************/
static ae_bool testalglibbasicsunit_testswapfunctions(ae_bool silent,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool okb1;
    ae_bool okb2;
    ae_bool oki1;
    ae_bool oki2;
    ae_bool okr1;
    ae_bool okr2;
    ae_bool okc1;
    ae_bool okc2;
    ae_vector b11;
    ae_vector b12;
    ae_vector i11;
    ae_vector i12;
    ae_vector r11;
    ae_vector r12;
    ae_vector c11;
    ae_vector c12;
    ae_matrix b21;
    ae_matrix b22;
    ae_matrix i21;
    ae_matrix i22;
    ae_matrix r21;
    ae_matrix r22;
    ae_matrix c21;
    ae_matrix c22;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&b11, 0, DT_BOOL, _state, ae_true);
    ae_vector_init(&b12, 0, DT_BOOL, _state, ae_true);
    ae_vector_init(&i11, 0, DT_INT, _state, ae_true);
    ae_vector_init(&i12, 0, DT_INT, _state, ae_true);
    ae_vector_init(&r11, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&r12, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c11, 0, DT_COMPLEX, _state, ae_true);
    ae_vector_init(&c12, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&b21, 0, 0, DT_BOOL, _state, ae_true);
    ae_matrix_init(&b22, 0, 0, DT_BOOL, _state, ae_true);
    ae_matrix_init(&i21, 0, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&i22, 0, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&r21, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&r22, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&c21, 0, 0, DT_COMPLEX, _state, ae_true);
    ae_matrix_init(&c22, 0, 0, DT_COMPLEX, _state, ae_true);

    result = ae_true;
    okb1 = ae_true;
    okb2 = ae_true;
    oki1 = ae_true;
    oki2 = ae_true;
    okr1 = ae_true;
    okr2 = ae_true;
    okc1 = ae_true;
    okc2 = ae_true;
    
    /*
     * Test B1 swaps
     */
    ae_vector_set_length(&b11, 1, _state);
    ae_vector_set_length(&b12, 2, _state);
    b11.ptr.p_bool[0] = ae_true;
    b12.ptr.p_bool[0] = ae_false;
    b12.ptr.p_bool[1] = ae_true;
    ae_swap_vectors(&b11, &b12);
    if( b11.cnt==2&&b12.cnt==1 )
    {
        okb1 = okb1&&!b11.ptr.p_bool[0];
        okb1 = okb1&&b11.ptr.p_bool[1];
        okb1 = okb1&&b12.ptr.p_bool[0];
    }
    else
    {
        okb1 = ae_false;
    }
    
    /*
     * Test I1 swaps
     */
    ae_vector_set_length(&i11, 1, _state);
    ae_vector_set_length(&i12, 2, _state);
    i11.ptr.p_int[0] = 1;
    i12.ptr.p_int[0] = 2;
    i12.ptr.p_int[1] = 3;
    ae_swap_vectors(&i11, &i12);
    if( i11.cnt==2&&i12.cnt==1 )
    {
        oki1 = oki1&&i11.ptr.p_int[0]==2;
        oki1 = oki1&&i11.ptr.p_int[1]==3;
        oki1 = oki1&&i12.ptr.p_int[0]==1;
    }
    else
    {
        oki1 = ae_false;
    }
    
    /*
     * Test R1 swaps
     */
    ae_vector_set_length(&r11, 1, _state);
    ae_vector_set_length(&r12, 2, _state);
    r11.ptr.p_double[0] = 1.5;
    r12.ptr.p_double[0] = 2.5;
    r12.ptr.p_double[1] = 3.5;
    ae_swap_vectors(&r11, &r12);
    if( r11.cnt==2&&r12.cnt==1 )
    {
        okr1 = okr1&&ae_fp_eq(r11.ptr.p_double[0],2.5);
        okr1 = okr1&&ae_fp_eq(r11.ptr.p_double[1],3.5);
        okr1 = okr1&&ae_fp_eq(r12.ptr.p_double[0],1.5);
    }
    else
    {
        okr1 = ae_false;
    }
    
    /*
     * Test C1 swaps
     */
    ae_vector_set_length(&c11, 1, _state);
    ae_vector_set_length(&c12, 2, _state);
    c11.ptr.p_complex[0] = ae_complex_from_d(1);
    c12.ptr.p_complex[0] = ae_complex_from_d(2);
    c12.ptr.p_complex[1] = ae_complex_from_d(3);
    ae_swap_vectors(&c11, &c12);
    if( c11.cnt==2&&c12.cnt==1 )
    {
        okc1 = okc1&&ae_c_eq_d(c11.ptr.p_complex[0],2);
        okc1 = okc1&&ae_c_eq_d(c11.ptr.p_complex[1],3);
        okc1 = okc1&&ae_c_eq_d(c12.ptr.p_complex[0],1);
    }
    else
    {
        okc1 = ae_false;
    }
    
    /*
     * Test B2 swaps
     */
    ae_matrix_set_length(&b21, 1, 2, _state);
    ae_matrix_set_length(&b22, 2, 1, _state);
    b21.ptr.pp_bool[0][0] = ae_true;
    b21.ptr.pp_bool[0][1] = ae_false;
    b22.ptr.pp_bool[0][0] = ae_false;
    b22.ptr.pp_bool[1][0] = ae_true;
    ae_swap_matrices(&b21, &b22);
    if( ((b21.rows==2&&b21.cols==1)&&b22.rows==1)&&b22.cols==2 )
    {
        okb2 = okb2&&!b21.ptr.pp_bool[0][0];
        okb2 = okb2&&b21.ptr.pp_bool[1][0];
        okb2 = okb2&&b22.ptr.pp_bool[0][0];
        okb2 = okb2&&!b22.ptr.pp_bool[0][1];
    }
    else
    {
        okb2 = ae_false;
    }
    
    /*
     * Test I2 swaps
     */
    ae_matrix_set_length(&i21, 1, 2, _state);
    ae_matrix_set_length(&i22, 2, 1, _state);
    i21.ptr.pp_int[0][0] = 1;
    i21.ptr.pp_int[0][1] = 2;
    i22.ptr.pp_int[0][0] = 3;
    i22.ptr.pp_int[1][0] = 4;
    ae_swap_matrices(&i21, &i22);
    if( ((i21.rows==2&&i21.cols==1)&&i22.rows==1)&&i22.cols==2 )
    {
        oki2 = oki2&&i21.ptr.pp_int[0][0]==3;
        oki2 = oki2&&i21.ptr.pp_int[1][0]==4;
        oki2 = oki2&&i22.ptr.pp_int[0][0]==1;
        oki2 = oki2&&i22.ptr.pp_int[0][1]==2;
    }
    else
    {
        oki2 = ae_false;
    }
    
    /*
     * Test R2 swaps
     */
    ae_matrix_set_length(&r21, 1, 2, _state);
    ae_matrix_set_length(&r22, 2, 1, _state);
    r21.ptr.pp_double[0][0] = 1;
    r21.ptr.pp_double[0][1] = 2;
    r22.ptr.pp_double[0][0] = 3;
    r22.ptr.pp_double[1][0] = 4;
    ae_swap_matrices(&r21, &r22);
    if( ((r21.rows==2&&r21.cols==1)&&r22.rows==1)&&r22.cols==2 )
    {
        okr2 = okr2&&ae_fp_eq(r21.ptr.pp_double[0][0],3);
        okr2 = okr2&&ae_fp_eq(r21.ptr.pp_double[1][0],4);
        okr2 = okr2&&ae_fp_eq(r22.ptr.pp_double[0][0],1);
        okr2 = okr2&&ae_fp_eq(r22.ptr.pp_double[0][1],2);
    }
    else
    {
        okr2 = ae_false;
    }
    
    /*
     * Test C2 swaps
     */
    ae_matrix_set_length(&c21, 1, 2, _state);
    ae_matrix_set_length(&c22, 2, 1, _state);
    c21.ptr.pp_complex[0][0] = ae_complex_from_d(1);
    c21.ptr.pp_complex[0][1] = ae_complex_from_d(2);
    c22.ptr.pp_complex[0][0] = ae_complex_from_d(3);
    c22.ptr.pp_complex[1][0] = ae_complex_from_d(4);
    ae_swap_matrices(&c21, &c22);
    if( ((c21.rows==2&&c21.cols==1)&&c22.rows==1)&&c22.cols==2 )
    {
        okc2 = okc2&&ae_c_eq_d(c21.ptr.pp_complex[0][0],3);
        okc2 = okc2&&ae_c_eq_d(c21.ptr.pp_complex[1][0],4);
        okc2 = okc2&&ae_c_eq_d(c22.ptr.pp_complex[0][0],1);
        okc2 = okc2&&ae_c_eq_d(c22.ptr.pp_complex[0][1],2);
    }
    else
    {
        okc2 = ae_false;
    }
    
    /*
     * summary
     */
    result = result&&okb1;
    result = result&&okb2;
    result = result&&oki1;
    result = result&&oki2;
    result = result&&okr1;
    result = result&&okr2;
    result = result&&okc1;
    result = result&&okc2;
    if( !silent )
    {
        if( result )
        {
            printf("SWAPPING FUNCTIONS:                      OK\n");
        }
        else
        {
            printf("SWAPPING FUNCTIONS:                      FAILED\n");
        }
    }
    ae_frame_leave(_state);
    return result;
}


ae_bool _rec1_init(rec1* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->b1field, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->r1field, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->i1field, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->c1field, 0, DT_COMPLEX, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->b2field, 0, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->r2field, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->i2field, 0, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->c2field, 0, 0, DT_COMPLEX, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _rec1_init_copy(rec1* dst, rec1* src, ae_state *_state, ae_bool make_automatic)
{
    dst->bfield = src->bfield;
    dst->rfield = src->rfield;
    dst->ifield = src->ifield;
    dst->cfield = src->cfield;
    if( !ae_vector_init_copy(&dst->b1field, &src->b1field, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->r1field, &src->r1field, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->i1field, &src->i1field, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->c1field, &src->c1field, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->b2field, &src->b2field, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->r2field, &src->r2field, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->i2field, &src->i2field, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->c2field, &src->c2field, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _rec1_clear(rec1* p)
{
    ae_vector_clear(&p->b1field);
    ae_vector_clear(&p->r1field);
    ae_vector_clear(&p->i1field);
    ae_vector_clear(&p->c1field);
    ae_matrix_clear(&p->b2field);
    ae_matrix_clear(&p->r2field);
    ae_matrix_clear(&p->i2field);
    ae_matrix_clear(&p->c2field);
}


/*$ End $*/
