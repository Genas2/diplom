

#include <stdafx.h>
#include <stdio.h>
#include "testmlptrainunit.h"


/*$ Declarations $*/
static void testmlptrainunit_createnetwork(multilayerperceptron* network,
     ae_int_t nkind,
     double a1,
     double a2,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_state *_state);
static void testmlptrainunit_unsetnetwork(multilayerperceptron* network,
     ae_state *_state);
static void testmlptrainunit_testinformational(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);
static void testmlptrainunit_testprocessing(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);
static void testmlptrainunit_testgradient(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);
static void testmlptrainunit_testhessian(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);


/*$ Body $*/


ae_bool testmlptrain(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_int_t passcount;
    ae_int_t maxn;
    ae_int_t maxhid;
    ae_int_t info;
    ae_int_t nf;
    ae_int_t nhid;
    ae_int_t nl;
    ae_int_t nhid1;
    ae_int_t nhid2;
    ae_int_t nkind;
    ae_int_t i;
    ae_int_t j;
    multilayerperceptron network;
    multilayerperceptron network2;
    mlpreport rep;
    mlpcvreport cvrep;
    ae_int_t ncount;
    ae_matrix xy;
    ae_matrix valxy;
    ae_int_t ssize;
    ae_int_t valsize;
    ae_bool allsame;
    ae_bool inferrors;
    ae_bool procerrors;
    ae_bool graderrors;
    ae_bool hesserrors;
    ae_bool trnerrors;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_init(&network, _state, ae_true);
    _multilayerperceptron_init(&network2, _state, ae_true);
    _mlpreport_init(&rep, _state, ae_true);
    _mlpcvreport_init(&cvrep, _state, ae_true);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&valxy, 0, 0, DT_REAL, _state, ae_true);

    waserrors = ae_false;
    inferrors = ae_false;
    procerrors = ae_false;
    graderrors = ae_false;
    hesserrors = ae_false;
    trnerrors = ae_false;
    passcount = 10;
    maxn = 4;
    maxhid = 4;
    
    /*
     * General multilayer network tests
     */
    for(nf=1; nf<=maxn; nf++)
    {
        for(nl=1; nl<=maxn; nl++)
        {
            for(nhid1=0; nhid1<=maxhid; nhid1++)
            {
                for(nhid2=0; nhid2<=0; nhid2++)
                {
                    for(nkind=0; nkind<=3; nkind++)
                    {
                        
                        /*
                         *  Skip meaningless parameters combinations
                         */
                        if( nkind==1&&nl<2 )
                        {
                            continue;
                        }
                        if( nhid1==0&&nhid2!=0 )
                        {
                            continue;
                        }
                        
                        /*
                         * Tests
                         */
                        testmlptrainunit_testinformational(nkind, nf, nhid1, nhid2, nl, passcount, &inferrors, _state);
                        testmlptrainunit_testprocessing(nkind, nf, nhid1, nhid2, nl, passcount, &procerrors, _state);
                        testmlptrainunit_testgradient(nkind, nf, nhid1, nhid2, nl, passcount, &graderrors, _state);
                        testmlptrainunit_testhessian(nkind, nf, nhid1, nhid2, nl, passcount, &hesserrors, _state);
                    }
                }
            }
        }
    }
    
    /*
     * Test network training on simple XOR problem
     */
    ae_matrix_set_length(&xy, 3+1, 2+1, _state);
    xy.ptr.pp_double[0][0] = -1;
    xy.ptr.pp_double[0][1] = -1;
    xy.ptr.pp_double[0][2] = -1;
    xy.ptr.pp_double[1][0] = 1;
    xy.ptr.pp_double[1][1] = -1;
    xy.ptr.pp_double[1][2] = 1;
    xy.ptr.pp_double[2][0] = -1;
    xy.ptr.pp_double[2][1] = 1;
    xy.ptr.pp_double[2][2] = 1;
    xy.ptr.pp_double[3][0] = 1;
    xy.ptr.pp_double[3][1] = 1;
    xy.ptr.pp_double[3][2] = -1;
    mlpcreate1(2, 2, 1, &network, _state);
    mlptrainlm(&network, &xy, 4, 0.001, 10, &info, &rep, _state);
    trnerrors = trnerrors||ae_fp_greater(mlprmserror(&network, &xy, 4, _state),0.1);
    
    /*
     * Test CV on random noisy problem
     */
    ncount = 100;
    ae_matrix_set_length(&xy, ncount-1+1, 1+1, _state);
    for(i=0; i<=ncount-1; i++)
    {
        xy.ptr.pp_double[i][0] = 2*ae_randomreal(_state)-1;
        xy.ptr.pp_double[i][1] = ae_randominteger(4, _state);
    }
    mlpcreatec0(1, 4, &network, _state);
    mlpkfoldcvlm(&network, &xy, ncount, 0.001, 5, 10, &info, &rep, &cvrep, _state);
    
    /*
     * Final report
     */
    waserrors = (((inferrors||procerrors)||graderrors)||hesserrors)||trnerrors;
    if( !silent )
    {
        printf("MLP TEST\n");
        printf("INFORMATIONAL FUNCTIONS:                 ");
        if( !inferrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("BASIC PROCESSING:                        ");
        if( !procerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("GRADIENT CALCULATION:                    ");
        if( !graderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HESSIAN CALCULATION:                     ");
        if( !hesserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("TRAINING:                                ");
        if( !trnerrors )
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


/*************************************************************************
Network creation
*************************************************************************/
static void testmlptrainunit_createnetwork(multilayerperceptron* network,
     ae_int_t nkind,
     double a1,
     double a2,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_state *_state)
{


    ae_assert(((nin>0&&nhid1>=0)&&nhid2>=0)&&nout>0, "CreateNetwork error", _state);
    ae_assert(nhid1!=0||nhid2==0, "CreateNetwork error", _state);
    ae_assert(nkind!=1||nout>=2, "CreateNetwork error", _state);
    if( nhid1==0 )
    {
        
        /*
         * No hidden layers
         */
        if( nkind==0 )
        {
            mlpcreate0(nin, nout, network, _state);
        }
        else
        {
            if( nkind==1 )
            {
                mlpcreatec0(nin, nout, network, _state);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpcreateb0(nin, nout, a1, a2, network, _state);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpcreater0(nin, nout, a1, a2, network, _state);
                    }
                }
            }
        }
        return;
    }
    if( nhid2==0 )
    {
        
        /*
         * One hidden layer
         */
        if( nkind==0 )
        {
            mlpcreate1(nin, nhid1, nout, network, _state);
        }
        else
        {
            if( nkind==1 )
            {
                mlpcreatec1(nin, nhid1, nout, network, _state);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpcreateb1(nin, nhid1, nout, a1, a2, network, _state);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpcreater1(nin, nhid1, nout, a1, a2, network, _state);
                    }
                }
            }
        }
        return;
    }
    
    /*
     * Two hidden layers
     */
    if( nkind==0 )
    {
        mlpcreate2(nin, nhid1, nhid2, nout, network, _state);
    }
    else
    {
        if( nkind==1 )
        {
            mlpcreatec2(nin, nhid1, nhid2, nout, network, _state);
        }
        else
        {
            if( nkind==2 )
            {
                mlpcreateb2(nin, nhid1, nhid2, nout, a1, a2, network, _state);
            }
            else
            {
                if( nkind==3 )
                {
                    mlpcreater2(nin, nhid1, nhid2, nout, a1, a2, network, _state);
                }
            }
        }
    }
}


/*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************/
static void testmlptrainunit_unsetnetwork(multilayerperceptron* network,
     ae_state *_state)
{


    mlpcreate0(1, 1, network, _state);
}


/*************************************************************************
Iformational functions test
*************************************************************************/
static void testmlptrainunit_testinformational(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron network;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t wcount;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_init(&network, _state, ae_true);

    testmlptrainunit_createnetwork(&network, nkind, 0.0, 0.0, nin, nhid1, nhid2, nout, _state);
    mlpproperties(&network, &n1, &n2, &wcount, _state);
    *err = ((*err||n1!=nin)||n2!=nout)||wcount<=0;
    ae_frame_leave(_state);
}


/*************************************************************************
Processing functions test
*************************************************************************/
static void testmlptrainunit_testprocessing(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron network;
    multilayerperceptron network2;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t wcount;
    ae_bool zeronet;
    double a1;
    double a2;
    ae_int_t pass;
    ae_int_t i;
    ae_bool allsame;
    ae_int_t rlen;
    ae_vector x1;
    ae_vector x2;
    ae_vector y1;
    ae_vector y2;
    ae_vector ra;
    ae_vector ra2;
    double v;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_init(&network, _state, ae_true);
    _multilayerperceptron_init(&network2, _state, ae_true);
    ae_vector_init(&x1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ra, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ra2, 0, DT_REAL, _state, ae_true);

    ae_assert(passcount>=2, "PassCount<2!", _state);
    
    /*
     * Prepare network
     */
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = 2*ae_randomreal(_state)-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = a1+(2*ae_randominteger(2, _state)-1)*(0.1+0.9*ae_randomreal(_state));
    }
    testmlptrainunit_createnetwork(&network, nkind, a1, a2, nin, nhid1, nhid2, nout, _state);
    mlpproperties(&network, &n1, &n2, &wcount, _state);
    
    /*
     * Initialize arrays
     */
    ae_vector_set_length(&x1, nin-1+1, _state);
    ae_vector_set_length(&x2, nin-1+1, _state);
    ae_vector_set_length(&y1, nout-1+1, _state);
    ae_vector_set_length(&y2, nout-1+1, _state);
    
    /*
     * Main cycle
     */
    for(pass=1; pass<=passcount; pass++)
    {
        
        /*
         * Last run is made on zero network
         */
        mlprandomizefull(&network, _state);
        zeronet = ae_false;
        if( pass==passcount )
        {
            ae_v_muld(&network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), 0);
            zeronet = ae_true;
        }
        
        /*
         * Same inputs leads to same outputs
         */
        for(i=0; i<=nin-1; i++)
        {
            x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            x2.ptr.p_double[i] = x1.ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            y2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        mlpprocess(&network, &x1, &y1, _state);
        mlpprocess(&network, &x2, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||!allsame;
        
        /*
         * Same inputs on original network leads to same outputs
         * on copy created using MLPCopy
         */
        testmlptrainunit_unsetnetwork(&network2, _state);
        mlpcopy(&network, &network2, _state);
        for(i=0; i<=nin-1; i++)
        {
            x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            x2.ptr.p_double[i] = x1.ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            y2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        mlpprocess(&network, &x1, &y1, _state);
        mlpprocess(&network2, &x2, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||!allsame;
        
        /*
         * Same inputs on original network leads to same outputs
         * on copy created using MLPSerialize
         */
        testmlptrainunit_unsetnetwork(&network2, _state);
        mlpserialize(&network, &ra, &rlen, _state);
        ae_vector_set_length(&ra2, rlen-1+1, _state);
        for(i=0; i<=rlen-1; i++)
        {
            ra2.ptr.p_double[i] = ra.ptr.p_double[i];
        }
        mlpunserialize(&ra2, &network2, _state);
        for(i=0; i<=nin-1; i++)
        {
            x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            x2.ptr.p_double[i] = x1.ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            y2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        }
        mlpprocess(&network, &x1, &y1, _state);
        mlpprocess(&network2, &x2, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||!allsame;
        
        /*
         * Different inputs leads to different outputs (non-zero network)
         */
        if( !zeronet )
        {
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                x2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            for(i=0; i<=nout-1; i++)
            {
                y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                y2.ptr.p_double[i] = y1.ptr.p_double[i];
            }
            mlpprocess(&network, &x1, &y1, _state);
            mlpprocess(&network, &x2, &y2, _state);
            allsame = ae_true;
            for(i=0; i<=nout-1; i++)
            {
                allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
            }
            *err = *err||allsame;
        }
        
        /*
         * Randomization changes outputs (when inputs are unchanged, non-zero network)
         */
        if( !zeronet )
        {
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                x2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            for(i=0; i<=nout-1; i++)
            {
                y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                y2.ptr.p_double[i] = y1.ptr.p_double[i];
            }
            mlpcopy(&network, &network2, _state);
            mlprandomize(&network2, _state);
            mlpprocess(&network, &x1, &y1, _state);
            mlpprocess(&network2, &x1, &y2, _state);
            allsame = ae_true;
            for(i=0; i<=nout-1; i++)
            {
                allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
            }
            *err = *err||allsame;
        }
        
        /*
         * Full randomization changes outputs (when inputs are unchanged, non-zero network)
         */
        if( !zeronet )
        {
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                x2.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            for(i=0; i<=nout-1; i++)
            {
                y1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
                y2.ptr.p_double[i] = y1.ptr.p_double[i];
            }
            mlpcopy(&network, &network2, _state);
            mlprandomizefull(&network2, _state);
            mlpprocess(&network, &x1, &y1, _state);
            mlpprocess(&network2, &x1, &y2, _state);
            allsame = ae_true;
            for(i=0; i<=nout-1; i++)
            {
                allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
            }
            *err = *err||allsame;
        }
        
        /*
         * Normalization properties
         */
        if( nkind==1 )
        {
            
            /*
             * Classifier network outputs are normalized
             */
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            mlpprocess(&network, &x1, &y1, _state);
            v = 0;
            for(i=0; i<=nout-1; i++)
            {
                v = v+y1.ptr.p_double[i];
                *err = *err||ae_fp_less(y1.ptr.p_double[i],0);
            }
            *err = *err||ae_fp_greater(ae_fabs(v-1, _state),1000*ae_machineepsilon);
        }
        if( nkind==2 )
        {
            
            /*
             * B-type network outputs are bounded from above/below
             */
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            mlpprocess(&network, &x1, &y1, _state);
            for(i=0; i<=nout-1; i++)
            {
                if( ae_fp_greater_eq(a2,0) )
                {
                    *err = *err||ae_fp_less(y1.ptr.p_double[i],a1);
                }
                else
                {
                    *err = *err||ae_fp_greater(y1.ptr.p_double[i],a1);
                }
            }
        }
        if( nkind==3 )
        {
            
            /*
             * R-type network outputs are within [A1,A2] (or [A2,A1])
             */
            for(i=0; i<=nin-1; i++)
            {
                x1.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            mlpprocess(&network, &x1, &y1, _state);
            for(i=0; i<=nout-1; i++)
            {
                *err = (*err||ae_fp_less(y1.ptr.p_double[i],ae_minreal(a1, a2, _state)))||ae_fp_greater(y1.ptr.p_double[i],ae_maxreal(a1, a2, _state));
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Gradient functions test
*************************************************************************/
static void testmlptrainunit_testgradient(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron network;
    multilayerperceptron network2;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t wcount;
    ae_bool zeronet;
    double h;
    double etol;
    double a1;
    double a2;
    ae_int_t pass;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_bool allsame;
    ae_int_t ilen;
    ae_int_t rlen;
    ae_int_t ssize;
    ae_matrix xy;
    ae_vector grad1;
    ae_vector grad2;
    ae_vector x;
    ae_vector y;
    ae_vector x1;
    ae_vector x2;
    ae_vector y1;
    ae_vector y2;
    ae_vector ia;
    ae_vector ra;
    double v;
    double e;
    double e1;
    double e2;
    double v1;
    double v2;
    double v3;
    double v4;
    double wprev;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_init(&network, _state, ae_true);
    _multilayerperceptron_init(&network2, _state, ae_true);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&grad1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&grad2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ia, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ra, 0, DT_REAL, _state, ae_true);

    ae_assert(passcount>=2, "PassCount<2!", _state);
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = 2*ae_randomreal(_state)-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = a1+(2*ae_randominteger(2, _state)-1)*(0.1+0.9*ae_randomreal(_state));
    }
    testmlptrainunit_createnetwork(&network, nkind, a1, a2, nin, nhid1, nhid2, nout, _state);
    mlpproperties(&network, &n1, &n2, &wcount, _state);
    h = 0.0001;
    etol = 0.01;
    
    /*
     * Initialize
     */
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&x1, nin-1+1, _state);
    ae_vector_set_length(&x2, nin-1+1, _state);
    ae_vector_set_length(&y, nout-1+1, _state);
    ae_vector_set_length(&y1, nout-1+1, _state);
    ae_vector_set_length(&y2, nout-1+1, _state);
    ae_vector_set_length(&grad1, wcount-1+1, _state);
    ae_vector_set_length(&grad2, wcount-1+1, _state);
    
    /*
     * Process
     */
    for(pass=1; pass<=passcount; pass++)
    {
        mlprandomizefull(&network, _state);
        
        /*
         * Test error/gradient calculation (least squares)
         */
        ae_matrix_set_length(&xy, 0+1, nin+nout-1+1, _state);
        for(i=0; i<=nin-1; i++)
        {
            x.ptr.p_double[i] = 4*ae_randomreal(_state)-2;
        }
        ae_v_move(&xy.ptr.pp_double[0][0], 1, &x.ptr.p_double[0], 1, ae_v_len(0,nin-1));
        if( mlpissoftmax(&network, _state) )
        {
            for(i=0; i<=nout-1; i++)
            {
                y.ptr.p_double[i] = 0;
            }
            xy.ptr.pp_double[0][nin] = ae_randominteger(nout, _state);
            y.ptr.p_double[ae_round(xy.ptr.pp_double[0][nin], _state)] = 1;
        }
        else
        {
            for(i=0; i<=nout-1; i++)
            {
                y.ptr.p_double[i] = 4*ae_randomreal(_state)-2;
            }
            ae_v_move(&xy.ptr.pp_double[0][nin], 1, &y.ptr.p_double[0], 1, ae_v_len(nin,nin+nout-1));
        }
        mlpgrad(&network, &x, &y, &e, &grad2, _state);
        mlpprocess(&network, &x, &y2, _state);
        ae_v_sub(&y2.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
        v = ae_v_dotproduct(&y2.ptr.p_double[0], 1, &y2.ptr.p_double[0], 1, ae_v_len(0,nout-1));
        v = v/2;
        *err = *err||ae_fp_greater(ae_fabs((v-e)/v, _state),etol);
        *err = *err||ae_fp_greater(ae_fabs((mlperror(&network, &xy, 1, _state)-v)/v, _state),etol);
        for(i=0; i<=wcount-1; i++)
        {
            wprev = network.weights.ptr.p_double[i];
            network.weights.ptr.p_double[i] = wprev-2*h;
            mlpprocess(&network, &x, &y1, _state);
            ae_v_sub(&y1.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v1 = ae_v_dotproduct(&y1.ptr.p_double[0], 1, &y1.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v1 = v1/2;
            network.weights.ptr.p_double[i] = wprev-h;
            mlpprocess(&network, &x, &y1, _state);
            ae_v_sub(&y1.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v2 = ae_v_dotproduct(&y1.ptr.p_double[0], 1, &y1.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v2 = v2/2;
            network.weights.ptr.p_double[i] = wprev+h;
            mlpprocess(&network, &x, &y1, _state);
            ae_v_sub(&y1.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v3 = ae_v_dotproduct(&y1.ptr.p_double[0], 1, &y1.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v3 = v3/2;
            network.weights.ptr.p_double[i] = wprev+2*h;
            mlpprocess(&network, &x, &y1, _state);
            ae_v_sub(&y1.ptr.p_double[0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v4 = ae_v_dotproduct(&y1.ptr.p_double[0], 1, &y1.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            v4 = v4/2;
            network.weights.ptr.p_double[i] = wprev;
            grad1.ptr.p_double[i] = (v1-8*v2+8*v3-v4)/(12*h);
            if( ae_fp_greater(ae_fabs(grad1.ptr.p_double[i], _state),1.0E-3) )
            {
                *err = *err||ae_fp_greater(ae_fabs((grad2.ptr.p_double[i]-grad1.ptr.p_double[i])/grad1.ptr.p_double[i], _state),etol);
            }
            else
            {
                *err = *err||ae_fp_greater(ae_fabs(grad2.ptr.p_double[i]-grad1.ptr.p_double[i], _state),etol);
            }
        }
        
        /*
         * Test error/gradient calculation (natural).
         * Testing on non-random structure networks
         * (because NKind is representative only in that case).
         */
        ae_matrix_set_length(&xy, 0+1, nin+nout-1+1, _state);
        for(i=0; i<=nin-1; i++)
        {
            x.ptr.p_double[i] = 4*ae_randomreal(_state)-2;
        }
        ae_v_move(&xy.ptr.pp_double[0][0], 1, &x.ptr.p_double[0], 1, ae_v_len(0,nin-1));
        if( mlpissoftmax(&network, _state) )
        {
            for(i=0; i<=nout-1; i++)
            {
                y.ptr.p_double[i] = 0;
            }
            xy.ptr.pp_double[0][nin] = ae_randominteger(nout, _state);
            y.ptr.p_double[ae_round(xy.ptr.pp_double[0][nin], _state)] = 1;
        }
        else
        {
            for(i=0; i<=nout-1; i++)
            {
                y.ptr.p_double[i] = 4*ae_randomreal(_state)-2;
            }
            ae_v_move(&xy.ptr.pp_double[0][nin], 1, &y.ptr.p_double[0], 1, ae_v_len(nin,nin+nout-1));
        }
        mlpgradn(&network, &x, &y, &e, &grad2, _state);
        mlpprocess(&network, &x, &y2, _state);
        v = 0;
        if( nkind!=1 )
        {
            for(i=0; i<=nout-1; i++)
            {
                v = v+0.5*ae_sqr(y2.ptr.p_double[i]-y.ptr.p_double[i], _state);
            }
        }
        else
        {
            for(i=0; i<=nout-1; i++)
            {
                if( ae_fp_neq(y.ptr.p_double[i],0) )
                {
                    if( ae_fp_eq(y2.ptr.p_double[i],0) )
                    {
                        v = v+y.ptr.p_double[i]*ae_log(ae_maxrealnumber, _state);
                    }
                    else
                    {
                        v = v+y.ptr.p_double[i]*ae_log(y.ptr.p_double[i]/y2.ptr.p_double[i], _state);
                    }
                }
            }
        }
        *err = *err||ae_fp_greater(ae_fabs((v-e)/v, _state),etol);
        *err = *err||ae_fp_greater(ae_fabs((mlperrorn(&network, &xy, 1, _state)-v)/v, _state),etol);
        for(i=0; i<=wcount-1; i++)
        {
            wprev = network.weights.ptr.p_double[i];
            network.weights.ptr.p_double[i] = wprev+h;
            mlpprocess(&network, &x, &y2, _state);
            network.weights.ptr.p_double[i] = wprev-h;
            mlpprocess(&network, &x, &y1, _state);
            network.weights.ptr.p_double[i] = wprev;
            v = 0;
            if( nkind!=1 )
            {
                for(j=0; j<=nout-1; j++)
                {
                    v = v+0.5*(ae_sqr(y2.ptr.p_double[j]-y.ptr.p_double[j], _state)-ae_sqr(y1.ptr.p_double[j]-y.ptr.p_double[j], _state))/(2*h);
                }
            }
            else
            {
                for(j=0; j<=nout-1; j++)
                {
                    if( ae_fp_neq(y.ptr.p_double[j],0) )
                    {
                        if( ae_fp_eq(y2.ptr.p_double[j],0) )
                        {
                            v = v+y.ptr.p_double[j]*ae_log(ae_maxrealnumber, _state);
                        }
                        else
                        {
                            v = v+y.ptr.p_double[j]*ae_log(y.ptr.p_double[j]/y2.ptr.p_double[j], _state);
                        }
                        if( ae_fp_eq(y1.ptr.p_double[j],0) )
                        {
                            v = v-y.ptr.p_double[j]*ae_log(ae_maxrealnumber, _state);
                        }
                        else
                        {
                            v = v-y.ptr.p_double[j]*ae_log(y.ptr.p_double[j]/y1.ptr.p_double[j], _state);
                        }
                    }
                }
                v = v/(2*h);
            }
            grad1.ptr.p_double[i] = v;
            if( ae_fp_greater(ae_fabs(grad1.ptr.p_double[i], _state),1.0E-3) )
            {
                *err = *err||ae_fp_greater(ae_fabs((grad2.ptr.p_double[i]-grad1.ptr.p_double[i])/grad1.ptr.p_double[i], _state),etol);
            }
            else
            {
                *err = *err||ae_fp_greater(ae_fabs(grad2.ptr.p_double[i]-grad1.ptr.p_double[i], _state),etol);
            }
        }
        
        /*
         * Test gradient calculation: batch (least squares)
         */
        ssize = 1+ae_randominteger(10, _state);
        ae_matrix_set_length(&xy, ssize-1+1, nin+nout-1+1, _state);
        for(i=0; i<=wcount-1; i++)
        {
            grad1.ptr.p_double[i] = 0;
        }
        e1 = 0;
        for(i=0; i<=ssize-1; i++)
        {
            for(j=0; j<=nin-1; j++)
            {
                x1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
            }
            ae_v_move(&xy.ptr.pp_double[i][0], 1, &x1.ptr.p_double[0], 1, ae_v_len(0,nin-1));
            if( mlpissoftmax(&network, _state) )
            {
                for(j=0; j<=nout-1; j++)
                {
                    y1.ptr.p_double[j] = 0;
                }
                xy.ptr.pp_double[i][nin] = ae_randominteger(nout, _state);
                y1.ptr.p_double[ae_round(xy.ptr.pp_double[i][nin], _state)] = 1;
            }
            else
            {
                for(j=0; j<=nout-1; j++)
                {
                    y1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
                }
                ae_v_move(&xy.ptr.pp_double[i][nin], 1, &y1.ptr.p_double[0], 1, ae_v_len(nin,nin+nout-1));
            }
            mlpgrad(&network, &x1, &y1, &v, &grad2, _state);
            e1 = e1+v;
            ae_v_add(&grad1.ptr.p_double[0], 1, &grad2.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
        mlpgradbatch(&network, &xy, ssize, &e2, &grad2, _state);
        *err = *err||ae_fp_greater(ae_fabs(e1-e2, _state)/e1,0.01);
        for(i=0; i<=wcount-1; i++)
        {
            if( ae_fp_neq(grad1.ptr.p_double[i],0) )
            {
                *err = *err||ae_fp_greater(ae_fabs((grad2.ptr.p_double[i]-grad1.ptr.p_double[i])/grad1.ptr.p_double[i], _state),etol);
            }
            else
            {
                *err = *err||ae_fp_neq(grad2.ptr.p_double[i],grad1.ptr.p_double[i]);
            }
        }
        
        /*
         * Test gradient calculation: batch (natural error func)
         */
        ssize = 1+ae_randominteger(10, _state);
        ae_matrix_set_length(&xy, ssize-1+1, nin+nout-1+1, _state);
        for(i=0; i<=wcount-1; i++)
        {
            grad1.ptr.p_double[i] = 0;
        }
        e1 = 0;
        for(i=0; i<=ssize-1; i++)
        {
            for(j=0; j<=nin-1; j++)
            {
                x1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
            }
            ae_v_move(&xy.ptr.pp_double[i][0], 1, &x1.ptr.p_double[0], 1, ae_v_len(0,nin-1));
            if( mlpissoftmax(&network, _state) )
            {
                for(j=0; j<=nout-1; j++)
                {
                    y1.ptr.p_double[j] = 0;
                }
                xy.ptr.pp_double[i][nin] = ae_randominteger(nout, _state);
                y1.ptr.p_double[ae_round(xy.ptr.pp_double[i][nin], _state)] = 1;
            }
            else
            {
                for(j=0; j<=nout-1; j++)
                {
                    y1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
                }
                ae_v_move(&xy.ptr.pp_double[i][nin], 1, &y1.ptr.p_double[0], 1, ae_v_len(nin,nin+nout-1));
            }
            mlpgradn(&network, &x1, &y1, &v, &grad2, _state);
            e1 = e1+v;
            ae_v_add(&grad1.ptr.p_double[0], 1, &grad2.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
        mlpgradnbatch(&network, &xy, ssize, &e2, &grad2, _state);
        *err = *err||ae_fp_greater(ae_fabs(e1-e2, _state)/e1,etol);
        for(i=0; i<=wcount-1; i++)
        {
            if( ae_fp_neq(grad1.ptr.p_double[i],0) )
            {
                *err = *err||ae_fp_greater(ae_fabs((grad2.ptr.p_double[i]-grad1.ptr.p_double[i])/grad1.ptr.p_double[i], _state),etol);
            }
            else
            {
                *err = *err||ae_fp_neq(grad2.ptr.p_double[i],grad1.ptr.p_double[i]);
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Hessian functions test
*************************************************************************/
static void testmlptrainunit_testhessian(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron network;
    multilayerperceptron network2;
    ae_int_t hkind;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t wcount;
    ae_bool zeronet;
    double h;
    double etol;
    ae_int_t pass;
    ae_int_t i;
    ae_int_t j;
    ae_bool allsame;
    ae_int_t ilen;
    ae_int_t rlen;
    ae_int_t ssize;
    double a1;
    double a2;
    ae_matrix xy;
    ae_matrix h1;
    ae_matrix h2;
    ae_vector grad1;
    ae_vector grad2;
    ae_vector grad3;
    ae_vector x;
    ae_vector y;
    ae_vector x1;
    ae_vector x2;
    ae_vector y1;
    ae_vector y2;
    ae_vector ia;
    ae_vector ra;
    double v;
    double e;
    double e1;
    double e2;
    double v1;
    double v2;
    double v3;
    double v4;
    double wprev;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_init(&network, _state, ae_true);
    _multilayerperceptron_init(&network2, _state, ae_true);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&h1, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&h2, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&grad1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&grad2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&grad3, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ia, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ra, 0, DT_REAL, _state, ae_true);

    ae_assert(passcount>=2, "PassCount<2!", _state);
    a1 = 0;
    a2 = 0;
    if( nkind==2 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = 2*ae_randomreal(_state)-1;
    }
    if( nkind==3 )
    {
        a1 = 1000*ae_randomreal(_state)-500;
        a2 = a1+(2*ae_randominteger(2, _state)-1)*(0.1+0.9*ae_randomreal(_state));
    }
    testmlptrainunit_createnetwork(&network, nkind, a1, a2, nin, nhid1, nhid2, nout, _state);
    mlpproperties(&network, &n1, &n2, &wcount, _state);
    h = 0.0001;
    etol = 0.05;
    
    /*
     * Initialize
     */
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&x1, nin-1+1, _state);
    ae_vector_set_length(&x2, nin-1+1, _state);
    ae_vector_set_length(&y, nout-1+1, _state);
    ae_vector_set_length(&y1, nout-1+1, _state);
    ae_vector_set_length(&y2, nout-1+1, _state);
    ae_vector_set_length(&grad1, wcount-1+1, _state);
    ae_vector_set_length(&grad2, wcount-1+1, _state);
    ae_vector_set_length(&grad3, wcount-1+1, _state);
    ae_matrix_set_length(&h1, wcount-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&h2, wcount-1+1, wcount-1+1, _state);
    
    /*
     * Process
     */
    for(pass=1; pass<=passcount; pass++)
    {
        mlprandomizefull(&network, _state);
        
        /*
         * Test hessian calculation .
         * E1 contains total error (calculated using MLPGrad/MLPGradN)
         * Grad1 contains total gradient (calculated using MLPGrad/MLPGradN)
         * H1 contains Hessian calculated using differences of gradients
         *
         * E2, Grad2 and H2 contains corresponing values calculated using MLPHessianBatch/MLPHessianNBatch
         */
        for(hkind=0; hkind<=1; hkind++)
        {
            ssize = 1+ae_randominteger(10, _state);
            ae_matrix_set_length(&xy, ssize-1+1, nin+nout-1+1, _state);
            for(i=0; i<=wcount-1; i++)
            {
                grad1.ptr.p_double[i] = 0;
            }
            for(i=0; i<=wcount-1; i++)
            {
                for(j=0; j<=wcount-1; j++)
                {
                    h1.ptr.pp_double[i][j] = 0;
                }
            }
            e1 = 0;
            for(i=0; i<=ssize-1; i++)
            {
                
                /*
                 * X, Y
                 */
                for(j=0; j<=nin-1; j++)
                {
                    x1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
                }
                ae_v_move(&xy.ptr.pp_double[i][0], 1, &x1.ptr.p_double[0], 1, ae_v_len(0,nin-1));
                if( mlpissoftmax(&network, _state) )
                {
                    for(j=0; j<=nout-1; j++)
                    {
                        y1.ptr.p_double[j] = 0;
                    }
                    xy.ptr.pp_double[i][nin] = ae_randominteger(nout, _state);
                    y1.ptr.p_double[ae_round(xy.ptr.pp_double[i][nin], _state)] = 1;
                }
                else
                {
                    for(j=0; j<=nout-1; j++)
                    {
                        y1.ptr.p_double[j] = 4*ae_randomreal(_state)-2;
                    }
                    ae_v_move(&xy.ptr.pp_double[i][nin], 1, &y1.ptr.p_double[0], 1, ae_v_len(nin,nin+nout-1));
                }
                
                /*
                 * E1, Grad1
                 */
                if( hkind==0 )
                {
                    mlpgrad(&network, &x1, &y1, &v, &grad2, _state);
                }
                else
                {
                    mlpgradn(&network, &x1, &y1, &v, &grad2, _state);
                }
                e1 = e1+v;
                ae_v_add(&grad1.ptr.p_double[0], 1, &grad2.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                
                /*
                 * H1
                 */
                for(j=0; j<=wcount-1; j++)
                {
                    wprev = network.weights.ptr.p_double[j];
                    network.weights.ptr.p_double[j] = wprev-2*h;
                    if( hkind==0 )
                    {
                        mlpgrad(&network, &x1, &y1, &v, &grad2, _state);
                    }
                    else
                    {
                        mlpgradn(&network, &x1, &y1, &v, &grad2, _state);
                    }
                    network.weights.ptr.p_double[j] = wprev-h;
                    if( hkind==0 )
                    {
                        mlpgrad(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    else
                    {
                        mlpgradn(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    ae_v_subd(&grad2.ptr.p_double[0], 1, &grad3.ptr.p_double[0], 1, ae_v_len(0,wcount-1), 8);
                    network.weights.ptr.p_double[j] = wprev+h;
                    if( hkind==0 )
                    {
                        mlpgrad(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    else
                    {
                        mlpgradn(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    ae_v_addd(&grad2.ptr.p_double[0], 1, &grad3.ptr.p_double[0], 1, ae_v_len(0,wcount-1), 8);
                    network.weights.ptr.p_double[j] = wprev+2*h;
                    if( hkind==0 )
                    {
                        mlpgrad(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    else
                    {
                        mlpgradn(&network, &x1, &y1, &v, &grad3, _state);
                    }
                    ae_v_sub(&grad2.ptr.p_double[0], 1, &grad3.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                    v = 1/(12*h);
                    ae_v_addd(&h1.ptr.pp_double[j][0], 1, &grad2.ptr.p_double[0], 1, ae_v_len(0,wcount-1), v);
                    network.weights.ptr.p_double[j] = wprev;
                }
            }
            if( hkind==0 )
            {
                mlphessianbatch(&network, &xy, ssize, &e2, &grad2, &h2, _state);
            }
            else
            {
                mlphessiannbatch(&network, &xy, ssize, &e2, &grad2, &h2, _state);
            }
            *err = *err||ae_fp_greater(ae_fabs(e1-e2, _state)/e1,etol);
            for(i=0; i<=wcount-1; i++)
            {
                if( ae_fp_greater(ae_fabs(grad1.ptr.p_double[i], _state),1.0E-2) )
                {
                    *err = *err||ae_fp_greater(ae_fabs((grad2.ptr.p_double[i]-grad1.ptr.p_double[i])/grad1.ptr.p_double[i], _state),etol);
                }
                else
                {
                    *err = *err||ae_fp_greater(ae_fabs(grad2.ptr.p_double[i]-grad1.ptr.p_double[i], _state),etol);
                }
            }
            for(i=0; i<=wcount-1; i++)
            {
                for(j=0; j<=wcount-1; j++)
                {
                    if( ae_fp_greater(ae_fabs(h1.ptr.pp_double[i][j], _state),5.0E-2) )
                    {
                        *err = *err||ae_fp_greater(ae_fabs((h1.ptr.pp_double[i][j]-h2.ptr.pp_double[i][j])/h1.ptr.pp_double[i][j], _state),etol);
                    }
                    else
                    {
                        *err = *err||ae_fp_greater(ae_fabs(h2.ptr.pp_double[i][j]-h1.ptr.pp_double[i][j], _state),etol);
                    }
                }
            }
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
