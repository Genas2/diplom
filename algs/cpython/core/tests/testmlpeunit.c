

#include <stdafx.h>
#include <stdio.h>
#include "testmlpeunit.h"


/*$ Declarations $*/
static void testmlpeunit_createensemble(mlpensemble* ensemble,
     ae_int_t nkind,
     double a1,
     double a2,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
     ae_state *_state);
static void testmlpeunit_unsetensemble(mlpensemble* ensemble,
     ae_state *_state);
static void testmlpeunit_testinformational(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);
static void testmlpeunit_testprocessing(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state);


/*$ Body $*/


ae_bool testmlpe(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool waserrors;
    ae_int_t passcount;
    ae_int_t maxn;
    ae_int_t maxhid;
    ae_int_t nf;
    ae_int_t nhid;
    ae_int_t nl;
    ae_int_t nhid1;
    ae_int_t nhid2;
    ae_int_t ec;
    ae_int_t nkind;
    ae_int_t algtype;
    ae_int_t tasktype;
    ae_int_t pass;
    mlpensemble ensemble;
    mlpreport rep;
    mlpcvreport oobrep;
    ae_matrix xy;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t npoints;
    double e;
    ae_int_t info;
    ae_int_t nless;
    ae_int_t nall;
    ae_int_t nclasses;
    ae_bool allsame;
    ae_bool inferrors;
    ae_bool procerrors;
    ae_bool trnerrors;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_init(&ensemble, _state, ae_true);
    _mlpreport_init(&rep, _state, ae_true);
    _mlpcvreport_init(&oobrep, _state, ae_true);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);

    waserrors = ae_false;
    inferrors = ae_false;
    procerrors = ae_false;
    trnerrors = ae_false;
    passcount = 10;
    maxn = 4;
    maxhid = 4;
    
    /*
     * General MLP ensembles tests
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
                        for(ec=1; ec<=3; ec++)
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
                            testmlpeunit_testinformational(nkind, nf, nhid1, nhid2, nl, ec, passcount, &inferrors, _state);
                            testmlpeunit_testprocessing(nkind, nf, nhid1, nhid2, nl, ec, passcount, &procerrors, _state);
                        }
                    }
                }
            }
        }
    }
    
    /*
     * network training must reduce error
     * test on random regression task
     */
    nin = 3;
    nout = 2;
    nhid = 5;
    npoints = 100;
    nless = 0;
    nall = 0;
    for(pass=1; pass<=10; pass++)
    {
        for(algtype=0; algtype<=1; algtype++)
        {
            for(tasktype=0; tasktype<=1; tasktype++)
            {
                if( tasktype==0 )
                {
                    ae_matrix_set_length(&xy, npoints-1+1, nin+nout-1+1, _state);
                    for(i=0; i<=npoints-1; i++)
                    {
                        for(j=0; j<=nin+nout-1; j++)
                        {
                            xy.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                        }
                    }
                    mlpecreate1(nin, nhid, nout, 1+ae_randominteger(3, _state), &ensemble, _state);
                }
                else
                {
                    ae_matrix_set_length(&xy, npoints-1+1, nin+1, _state);
                    nclasses = 2+ae_randominteger(2, _state);
                    for(i=0; i<=npoints-1; i++)
                    {
                        for(j=0; j<=nin-1; j++)
                        {
                            xy.ptr.pp_double[i][j] = 2*ae_randomreal(_state)-1;
                        }
                        xy.ptr.pp_double[i][nin] = ae_randominteger(nclasses, _state);
                    }
                    mlpecreatec1(nin, nhid, nclasses, 1+ae_randominteger(3, _state), &ensemble, _state);
                }
                e = mlpermserror(&ensemble, &xy, npoints, _state);
                if( algtype==0 )
                {
                    mlpebagginglm(&ensemble, &xy, npoints, 0.001, 1, &info, &rep, &oobrep, _state);
                }
                else
                {
                    mlpebagginglbfgs(&ensemble, &xy, npoints, 0.001, 1, 0.01, 0, &info, &rep, &oobrep, _state);
                }
                if( info<0 )
                {
                    trnerrors = ae_true;
                }
                else
                {
                    if( ae_fp_less(mlpermserror(&ensemble, &xy, npoints, _state),e) )
                    {
                        nless = nless+1;
                    }
                }
                nall = nall+1;
            }
        }
    }
    trnerrors = trnerrors||ae_fp_greater(nall-nless,0.3*nall);
    
    /*
     * Final report
     */
    waserrors = (inferrors||procerrors)||trnerrors;
    if( !silent )
    {
        printf("MLP ENSEMBLE TEST\n");
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
static void testmlpeunit_createensemble(mlpensemble* ensemble,
     ae_int_t nkind,
     double a1,
     double a2,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
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
            mlpecreate0(nin, nout, ec, ensemble, _state);
        }
        else
        {
            if( nkind==1 )
            {
                mlpecreatec0(nin, nout, ec, ensemble, _state);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpecreateb0(nin, nout, a1, a2, ec, ensemble, _state);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpecreater0(nin, nout, a1, a2, ec, ensemble, _state);
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
            mlpecreate1(nin, nhid1, nout, ec, ensemble, _state);
        }
        else
        {
            if( nkind==1 )
            {
                mlpecreatec1(nin, nhid1, nout, ec, ensemble, _state);
            }
            else
            {
                if( nkind==2 )
                {
                    mlpecreateb1(nin, nhid1, nout, a1, a2, ec, ensemble, _state);
                }
                else
                {
                    if( nkind==3 )
                    {
                        mlpecreater1(nin, nhid1, nout, a1, a2, ec, ensemble, _state);
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
        mlpecreate2(nin, nhid1, nhid2, nout, ec, ensemble, _state);
    }
    else
    {
        if( nkind==1 )
        {
            mlpecreatec2(nin, nhid1, nhid2, nout, ec, ensemble, _state);
        }
        else
        {
            if( nkind==2 )
            {
                mlpecreateb2(nin, nhid1, nhid2, nout, a1, a2, ec, ensemble, _state);
            }
            else
            {
                if( nkind==3 )
                {
                    mlpecreater2(nin, nhid1, nhid2, nout, a1, a2, ec, ensemble, _state);
                }
            }
        }
    }
}


/*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************/
static void testmlpeunit_unsetensemble(mlpensemble* ensemble,
     ae_state *_state)
{


    mlpecreate0(1, 1, 1, ensemble, _state);
}


/*************************************************************************
Iformational functions test
*************************************************************************/
static void testmlpeunit_testinformational(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    mlpensemble ensemble;
    ae_int_t n1;
    ae_int_t n2;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_init(&ensemble, _state, ae_true);

    testmlpeunit_createensemble(&ensemble, nkind, -1.0, 1.0, nin, nhid1, nhid2, nout, ec, _state);
    mlpeproperties(&ensemble, &n1, &n2, _state);
    *err = (*err||n1!=nin)||n2!=nout;
    ae_frame_leave(_state);
}


/*************************************************************************
Processing functions test
*************************************************************************/
static void testmlpeunit_testprocessing(ae_int_t nkind,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ec,
     ae_int_t passcount,
     ae_bool* err,
     ae_state *_state)
{
    ae_frame _frame_block;
    mlpensemble ensemble;
    mlpensemble ensemble2;
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
    _mlpensemble_init(&ensemble, _state, ae_true);
    _mlpensemble_init(&ensemble2, _state, ae_true);
    ae_vector_init(&x1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ra, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ra2, 0, DT_REAL, _state, ae_true);

    
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
        testmlpeunit_createensemble(&ensemble, nkind, a1, a2, nin, nhid1, nhid2, nout, ec, _state);
        
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
        mlpeprocess(&ensemble, &x1, &y1, _state);
        mlpeprocess(&ensemble, &x2, &y2, _state);
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
        testmlpeunit_unsetensemble(&ensemble2, _state);
        mlpecopy(&ensemble, &ensemble2, _state);
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
        mlpeprocess(&ensemble, &x1, &y1, _state);
        mlpeprocess(&ensemble2, &x2, &y2, _state);
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
        testmlpeunit_unsetensemble(&ensemble2, _state);
        mlpeserialize(&ensemble, &ra, &rlen, _state);
        ae_vector_set_length(&ra2, rlen-1+1, _state);
        for(i=0; i<=rlen-1; i++)
        {
            ra2.ptr.p_double[i] = ra.ptr.p_double[i];
        }
        mlpeunserialize(&ra2, &ensemble2, _state);
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
        mlpeprocess(&ensemble, &x1, &y1, _state);
        mlpeprocess(&ensemble2, &x2, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||!allsame;
        
        /*
         * Different inputs leads to different outputs (non-zero network)
         */
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
        mlpeprocess(&ensemble, &x1, &y1, _state);
        mlpeprocess(&ensemble, &x2, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||allsame;
        
        /*
         * Randomization changes outputs (when inputs are unchanged, non-zero network)
         */
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
        mlpecopy(&ensemble, &ensemble2, _state);
        mlperandomize(&ensemble2, _state);
        mlpeprocess(&ensemble, &x1, &y1, _state);
        mlpeprocess(&ensemble2, &x1, &y2, _state);
        allsame = ae_true;
        for(i=0; i<=nout-1; i++)
        {
            allsame = allsame&&ae_fp_eq(y1.ptr.p_double[i],y2.ptr.p_double[i]);
        }
        *err = *err||allsame;
        
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
            mlpeprocess(&ensemble, &x1, &y1, _state);
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
            mlpeprocess(&ensemble, &x1, &y1, _state);
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
            mlpeprocess(&ensemble, &x1, &y1, _state);
            for(i=0; i<=nout-1; i++)
            {
                *err = (*err||ae_fp_less(y1.ptr.p_double[i],ae_minreal(a1, a2, _state)))||ae_fp_greater(y1.ptr.p_double[i],ae_maxreal(a1, a2, _state));
            }
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
