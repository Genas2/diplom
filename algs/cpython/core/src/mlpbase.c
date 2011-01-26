/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/


#include <stdafx.h>
#include "mlpbase.h"


/*$ Declarations $*/
static ae_int_t mlpbase_mlpvnum = 7;
static ae_int_t mlpbase_nfieldwidth = 4;
static ae_int_t mlpbase_chunksize = 32;
static void mlpbase_addinputlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addbiasedsummatorlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addactivationlayer(ae_int_t functype,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addzerolayer(/* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_mlpcreate(ae_int_t nin,
     ae_int_t nout,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t layerscount,
     ae_bool isclsnet,
     multilayerperceptron* network,
     ae_state *_state);
static void mlpbase_mlpactivationfunction(double net,
     ae_int_t k,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state);
static void mlpbase_mlphessianbatchinternal(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_bool naturalerr,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
static void mlpbase_mlpinternalcalculategradient(multilayerperceptron* network,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* derror,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state);
static void mlpbase_mlpchunkedgradient(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state);
static double mlpbase_safecrossentropy(double t,
     double z,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+2;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+2;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+2;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values.

Activation function of the output layer takes values:

    (B, +INF), if D>=0

or

    (-INF, B), if D<0.


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC0: NOut<2!", _state);
    layerscount = 1+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC1: NOut<2!", _state);
    layerscount = 1+3+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC2: NOut<2!", _state);
    layerscount = 1+3+3+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Copying of neural network

INPUT PARAMETERS:
    Network1 -   original

OUTPUT PARAMETERS:
    Network2 -   copy

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcopy(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    _multilayerperceptron_clear(network2);

    
    /*
     * Unload info
     */
    ssize = network1->structinfo.ptr.p_int[0];
    nin = network1->structinfo.ptr.p_int[1];
    nout = network1->structinfo.ptr.p_int[2];
    ntotal = network1->structinfo.ptr.p_int[3];
    wcount = network1->structinfo.ptr.p_int[4];
    
    /*
     * Allocate space
     */
    ae_vector_set_length(&network2->structinfo, ssize-1+1, _state);
    ae_vector_set_length(&network2->weights, wcount-1+1, _state);
    if( mlpissoftmax(network1, _state) )
    {
        ae_vector_set_length(&network2->columnmeans, nin-1+1, _state);
        ae_vector_set_length(&network2->columnsigmas, nin-1+1, _state);
    }
    else
    {
        ae_vector_set_length(&network2->columnmeans, nin+nout-1+1, _state);
        ae_vector_set_length(&network2->columnsigmas, nin+nout-1+1, _state);
    }
    ae_vector_set_length(&network2->neurons, ntotal-1+1, _state);
    ae_matrix_set_length(&network2->chunks, 3*ntotal+1, mlpbase_chunksize-1+1, _state);
    ae_vector_set_length(&network2->nwbuf, ae_maxint(wcount, 2*nout, _state)-1+1, _state);
    ae_vector_set_length(&network2->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&network2->x, nin-1+1, _state);
    ae_vector_set_length(&network2->y, nout-1+1, _state);
    ae_vector_set_length(&network2->derror, ntotal-1+1, _state);
    
    /*
     * Copy
     */
    for(i=0; i<=ssize-1; i++)
    {
        network2->structinfo.ptr.p_int[i] = network1->structinfo.ptr.p_int[i];
    }
    ae_v_move(&network2->weights.ptr.p_double[0], 1, &network1->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    if( mlpissoftmax(network1, _state) )
    {
        ae_v_move(&network2->columnmeans.ptr.p_double[0], 1, &network1->columnmeans.ptr.p_double[0], 1, ae_v_len(0,nin-1));
        ae_v_move(&network2->columnsigmas.ptr.p_double[0], 1, &network1->columnsigmas.ptr.p_double[0], 1, ae_v_len(0,nin-1));
    }
    else
    {
        ae_v_move(&network2->columnmeans.ptr.p_double[0], 1, &network1->columnmeans.ptr.p_double[0], 1, ae_v_len(0,nin+nout-1));
        ae_v_move(&network2->columnsigmas.ptr.p_double[0], 1, &network1->columnsigmas.ptr.p_double[0], 1, ae_v_len(0,nin+nout-1));
    }
    ae_v_move(&network2->neurons.ptr.p_double[0], 1, &network1->neurons.ptr.p_double[0], 1, ae_v_len(0,ntotal-1));
    ae_v_move(&network2->dfdnet.ptr.p_double[0], 1, &network1->dfdnet.ptr.p_double[0], 1, ae_v_len(0,ntotal-1));
    ae_v_move(&network2->x.ptr.p_double[0], 1, &network1->x.ptr.p_double[0], 1, ae_v_len(0,nin-1));
    ae_v_move(&network2->y.ptr.p_double[0], 1, &network1->y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
    ae_v_move(&network2->derror.ptr.p_double[0], 1, &network1->derror.ptr.p_double[0], 1, ae_v_len(0,ntotal-1));
}


/*************************************************************************
Serialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    Network -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores network,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpserialize(multilayerperceptron* network,
     /* Real    */ ae_vector* ra,
     ae_int_t* rlen,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t sigmalen;
    ae_int_t offs;

    ae_vector_clear(ra);
    *rlen = 0;

    
    /*
     * Unload info
     */
    ssize = network->structinfo.ptr.p_int[0];
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    wcount = network->structinfo.ptr.p_int[4];
    if( mlpissoftmax(network, _state) )
    {
        sigmalen = nin;
    }
    else
    {
        sigmalen = nin+nout;
    }
    
    /*
     *  RA format:
     *      LEN         DESRC.
     *      1           RLen
     *      1           version (MLPVNum)
     *      1           StructInfo size
     *      SSize       StructInfo
     *      WCount      Weights
     *      SigmaLen    ColumnMeans
     *      SigmaLen    ColumnSigmas
     */
    *rlen = 3+ssize+wcount+2*sigmalen;
    ae_vector_set_length(ra, *rlen-1+1, _state);
    ra->ptr.p_double[0] = *rlen;
    ra->ptr.p_double[1] = mlpbase_mlpvnum;
    ra->ptr.p_double[2] = ssize;
    offs = 3;
    for(i=0; i<=ssize-1; i++)
    {
        ra->ptr.p_double[offs+i] = network->structinfo.ptr.p_int[i];
    }
    offs = offs+ssize;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(offs,offs+wcount-1));
    offs = offs+wcount;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->columnmeans.ptr.p_double[0], 1, ae_v_len(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->columnsigmas.ptr.p_double[0], 1, ae_v_len(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Unserialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    RA      -   real array which stores network

OUTPUT PARAMETERS:
    Network -   restored network

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpunserialize(/* Real    */ ae_vector* ra,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t sigmalen;
    ae_int_t offs;

    _multilayerperceptron_clear(network);

    ae_assert(ae_round(ra->ptr.p_double[1], _state)==mlpbase_mlpvnum, "MLPUnserialize: incorrect array!", _state);
    
    /*
     * Unload StructInfo from IA
     */
    offs = 3;
    ssize = ae_round(ra->ptr.p_double[2], _state);
    ae_vector_set_length(&network->structinfo, ssize-1+1, _state);
    for(i=0; i<=ssize-1; i++)
    {
        network->structinfo.ptr.p_int[i] = ae_round(ra->ptr.p_double[offs+i], _state);
    }
    offs = offs+ssize;
    
    /*
     * Unload info from StructInfo
     */
    ssize = network->structinfo.ptr.p_int[0];
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    wcount = network->structinfo.ptr.p_int[4];
    if( network->structinfo.ptr.p_int[6]==0 )
    {
        sigmalen = nin+nout;
    }
    else
    {
        sigmalen = nin;
    }
    
    /*
     * Allocate space for other fields
     */
    ae_vector_set_length(&network->weights, wcount-1+1, _state);
    ae_vector_set_length(&network->columnmeans, sigmalen-1+1, _state);
    ae_vector_set_length(&network->columnsigmas, sigmalen-1+1, _state);
    ae_vector_set_length(&network->neurons, ntotal-1+1, _state);
    ae_matrix_set_length(&network->chunks, 3*ntotal+1, mlpbase_chunksize-1+1, _state);
    ae_vector_set_length(&network->nwbuf, ae_maxint(wcount, 2*nout, _state)-1+1, _state);
    ae_vector_set_length(&network->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&network->x, nin-1+1, _state);
    ae_vector_set_length(&network->y, nout-1+1, _state);
    ae_vector_set_length(&network->derror, ntotal-1+1, _state);
    
    /*
     * Copy parameters from RA
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,wcount-1));
    offs = offs+wcount;
    ae_v_move(&network->columnmeans.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,sigmalen-1));
    offs = offs+sigmalen;
    ae_v_move(&network->columnsigmas.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(multilayerperceptron* network, ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
}


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(multilayerperceptron* network, ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t offs;
    ae_int_t ntype;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Process network
     */
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        network->columnsigmas.ptr.p_double[i] = 1.5*ae_randomreal(_state)+0.5;
    }
    if( !mlpissoftmax(network, _state) )
    {
        for(i=0; i<=nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*mlpbase_nfieldwidth;
            ntype = network->structinfo.ptr.p_int[offs+0];
            if( ntype==0 )
            {
                
                /*
                 * Shifts are changed only for linear outputs neurons
                 */
                network->columnmeans.ptr.p_double[nin+i] = 2*ae_randomreal(_state)-1;
            }
            if( ntype==0||ntype==3 )
            {
                
                /*
                 * Scales are changed only for linear or bounded outputs neurons.
                 * Note that scale randomization preserves sign.
                 */
                network->columnsigmas.ptr.p_double[nin+i] = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state)*(1.5*ae_randomreal(_state)+0.5);
            }
        }
    }
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpinitpreprocessor(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t jmax;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t offs;
    ae_int_t ntype;
    ae_vector means;
    ae_vector sigmas;
    double s;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&means, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sigmas, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Means/Sigmas
     */
    if( mlpissoftmax(network, _state) )
    {
        jmax = nin-1;
    }
    else
    {
        jmax = nin+nout-1;
    }
    ae_vector_set_length(&means, jmax+1, _state);
    ae_vector_set_length(&sigmas, jmax+1, _state);
    for(j=0; j<=jmax; j++)
    {
        means.ptr.p_double[j] = 0;
        for(i=0; i<=ssize-1; i++)
        {
            means.ptr.p_double[j] = means.ptr.p_double[j]+xy->ptr.pp_double[i][j];
        }
        means.ptr.p_double[j] = means.ptr.p_double[j]/ssize;
        sigmas.ptr.p_double[j] = 0;
        for(i=0; i<=ssize-1; i++)
        {
            sigmas.ptr.p_double[j] = sigmas.ptr.p_double[j]+ae_sqr(xy->ptr.pp_double[i][j]-means.ptr.p_double[j], _state);
        }
        sigmas.ptr.p_double[j] = ae_sqrt(sigmas.ptr.p_double[j]/ssize, _state);
    }
    
    /*
     * Inputs
     */
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = means.ptr.p_double[i];
        network->columnsigmas.ptr.p_double[i] = sigmas.ptr.p_double[i];
        if( ae_fp_eq(network->columnsigmas.ptr.p_double[i],0) )
        {
            network->columnsigmas.ptr.p_double[i] = 1;
        }
    }
    
    /*
     * Outputs
     */
    if( !mlpissoftmax(network, _state) )
    {
        for(i=0; i<=nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*mlpbase_nfieldwidth;
            ntype = network->structinfo.ptr.p_int[offs+0];
            
            /*
             * Linear outputs
             */
            if( ntype==0 )
            {
                network->columnmeans.ptr.p_double[nin+i] = means.ptr.p_double[nin+i];
                network->columnsigmas.ptr.p_double[nin+i] = sigmas.ptr.p_double[nin+i];
                if( ae_fp_eq(network->columnsigmas.ptr.p_double[nin+i],0) )
                {
                    network->columnsigmas.ptr.p_double[nin+i] = 1;
                }
            }
            
            /*
             * Bounded outputs (half-interval)
             */
            if( ntype==3 )
            {
                s = means.ptr.p_double[nin+i]-network->columnmeans.ptr.p_double[nin+i];
                if( ae_fp_eq(s,0) )
                {
                    s = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state);
                }
                if( ae_fp_eq(s,0) )
                {
                    s = 1.0;
                }
                network->columnsigmas.ptr.p_double[nin+i] = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state)*ae_fabs(s, _state);
                if( ae_fp_eq(network->columnsigmas.ptr.p_double[nin+i],0) )
                {
                    network->columnsigmas.ptr.p_double[nin+i] = 1;
                }
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(multilayerperceptron* network,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_int_t* wcount,
     ae_state *_state)
{

    *nin = 0;
    *nout = 0;
    *wcount = 0;

    *nin = network->structinfo.ptr.p_int[1];
    *nout = network->structinfo.ptr.p_int[2];
    *wcount = network->structinfo.ptr.p_int[4];
}


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_bool mlpissoftmax(multilayerperceptron* network, ae_state *_state)
{
    ae_bool result;


    result = network->structinfo.ptr.p_int[6]==1;
    return result;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Network -   neural network
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also MLPProcessI

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpprocess(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{


    if( y->cnt<network->structinfo.ptr.p_int[2] )
    {
        ae_vector_set_length(y, network->structinfo.ptr.p_int[2], _state);
    }
    mlpinternalprocessvector(&network->structinfo, &network->weights, &network->columnmeans, &network->columnsigmas, &network->neurons, &network->dfdnet, x, y, _state);
}


/*************************************************************************
'interactive'  variant  of  MLPProcess  for  languages  like  Python which
support constructs like "Y = MLPProcess(NN,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 21.09.2010 by Bochkanov Sergey
*************************************************************************/
void mlpprocessi(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    mlpprocess(network, x, y, _state);
}


/*************************************************************************
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double e;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels outputs
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            if( k>=0&&k<nout )
            {
                network->y.ptr.p_double[k] = network->y.ptr.p_double[k]-1;
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            ae_v_sub(&network->y.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1));
        }
        e = ae_v_dotproduct(&network->y.ptr.p_double[0], 1, &network->y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
        result = result+e/2;
    }
    return result;
}


/*************************************************************************
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double e;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        
        /*
         * Process vector
         */
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        
        /*
         * Update error function
         */
        if( network->structinfo.ptr.p_int[6]==0 )
        {
            
            /*
             * Least squares error function
             */
            ae_v_sub(&network->y.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1));
            e = ae_v_dotproduct(&network->y.ptr.p_double[0], 1, &network->y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            result = result+e/2;
        }
        else
        {
            
            /*
             * Cross-entropy error function
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            if( k>=0&&k<nout )
            {
                result = result+mlpbase_safecrossentropy(1, network->y.ptr.p_double[k], _state);
            }
        }
    }
    return result;
}


/*************************************************************************
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector workx;
    ae_vector worky;
    ae_int_t nn;
    ae_int_t ns;
    ae_int_t nmax;
    ae_int_t result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&worky, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ae_vector_set_length(&workx, nin-1+1, _state);
    ae_vector_set_length(&worky, nout-1+1, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        
        /*
         * Process
         */
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &workx, &worky, _state);
        
        /*
         * Network version of the answer
         */
        nmax = 0;
        for(j=0; j<=nout-1; j++)
        {
            if( ae_fp_greater(worky.ptr.p_double[j],worky.ptr.p_double[nmax]) )
            {
                nmax = j;
            }
        }
        nn = nmax;
        
        /*
         * Right answer
         */
        if( mlpissoftmax(network, _state) )
        {
            ns = ae_round(xy->ptr.pp_double[i][nin], _state);
        }
        else
        {
            nmax = 0;
            for(j=0; j<=nout-1; j++)
            {
                if( ae_fp_greater(xy->ptr.pp_double[i][nin+j],xy->ptr.pp_double[i][nin+nmax]) )
                {
                    nmax = j;
                }
            }
            ns = nmax;
        }
        
        /*
         * compare
         */
        if( nn!=ns )
        {
            result = result+1;
        }
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Network -   network
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases. Works both for
    classifier networks and general purpose networks used as
    classifiers.

  -- ALGLIB --
     Copyright 25.12.2008 by Bochkanov Sergey
*************************************************************************/
double mlprelclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double result;


    result = (double)mlpclserror(network, xy, npoints, _state)/(double)npoints;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if network solves regression task.

  -- ALGLIB --
     Copyright 08.01.2009 by Bochkanov Sergey
*************************************************************************/
double mlpavgce(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    if( mlpissoftmax(network, _state) )
    {
        mlpproperties(network, &nin, &nout, &wcount, _state);
        result = mlperrorn(network, xy, npoints, _state)/(npoints*ae_log(2, _state));
    }
    else
    {
        result = 0;
    }
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlprmserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = ae_sqrt(2*mlperror(network, xy, npoints, _state)/(npoints*nout), _state);
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            for(j=0; j<=nout-1; j++)
            {
                if( j==k )
                {
                    result = result+ae_fabs(1-network->y.ptr.p_double[j], _state);
                }
                else
                {
                    result = result+ae_fabs(network->y.ptr.p_double[j], _state);
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            for(j=0; j<=nout-1; j++)
            {
                result = result+ae_fabs(xy->ptr.pp_double[i][nin+j]-network->y.ptr.p_double[j], _state);
            }
        }
    }
    result = result/(npoints*nout);
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t lk;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    k = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels
             */
            lk = ae_round(xy->ptr.pp_double[i][nin], _state);
            for(j=0; j<=nout-1; j++)
            {
                if( j==lk )
                {
                    result = result+ae_fabs(1-network->y.ptr.p_double[j], _state);
                    k = k+1;
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            for(j=0; j<=nout-1; j++)
            {
                if( ae_fp_neq(xy->ptr.pp_double[i][nin+j],0) )
                {
                    result = result+ae_fabs(xy->ptr.pp_double[i][nin+j]-network->y.ptr.p_double[j], _state)/ae_fabs(xy->ptr.pp_double[i][nin+j], _state);
                    k = k+1;
                }
            }
        }
    }
    if( k!=0 )
    {
        result = result/k;
    }
    return result;
}


/*************************************************************************
Gradient calculation

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]
    
  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nout;
    ae_int_t ntotal;

    *e = 0;

    
    /*
     * Alloc
     */
    if( grad->cnt<network->structinfo.ptr.p_int[4] )
    {
        ae_vector_set_length(grad, network->structinfo.ptr.p_int[4], _state);
    }
    
    /*
     * Prepare dError/dOut, internal structures
     */
    mlpprocess(network, x, &network->y, _state);
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    *e = 0;
    for(i=0; i<=ntotal-1; i++)
    {
        network->derror.ptr.p_double[i] = 0;
    }
    for(i=0; i<=nout-1; i++)
    {
        network->derror.ptr.p_double[ntotal-nout+i] = network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
        *e = *e+ae_sqr(network->y.ptr.p_double[i]-desiredy->ptr.p_double[i], _state)/2;
    }
    
    /*
     * gradient
     */
    mlpbase_mlpinternalcalculategradient(network, &network->neurons, &network->weights, &network->derror, grad, ae_false, _state);
}


/*************************************************************************
Gradient calculation (natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    double s;
    ae_int_t i;
    ae_int_t nout;
    ae_int_t ntotal;

    *e = 0;

    
    /*
     * Alloc
     */
    if( grad->cnt<network->structinfo.ptr.p_int[4] )
    {
        ae_vector_set_length(grad, network->structinfo.ptr.p_int[4], _state);
    }
    
    /*
     * Prepare dError/dOut, internal structures
     */
    mlpprocess(network, x, &network->y, _state);
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    for(i=0; i<=ntotal-1; i++)
    {
        network->derror.ptr.p_double[i] = 0;
    }
    *e = 0;
    if( network->structinfo.ptr.p_int[6]==0 )
    {
        
        /*
         * Regression network, least squares
         */
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
            *e = *e+ae_sqr(network->y.ptr.p_double[i]-desiredy->ptr.p_double[i], _state)/2;
        }
    }
    else
    {
        
        /*
         * Classification network, cross-entropy
         */
        s = 0;
        for(i=0; i<=nout-1; i++)
        {
            s = s+desiredy->ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = s*network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
            *e = *e+mlpbase_safecrossentropy(desiredy->ptr.p_double[i], network->y.ptr.p_double[i], _state);
        }
    }
    
    /*
     * gradient
     */
    mlpbase_mlpinternalcalculategradient(network, &network->neurons, &network->weights, &network->derror, grad, ae_true, _state);
}


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    *e = 0;

    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        grad->ptr.p_double[i] = 0;
    }
    *e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpbase_mlpchunkedgradient(network, xy, i, ae_minint(ssize, i+mlpbase_chunksize, _state)-i, e, grad, ae_false, _state);
        i = i+mlpbase_chunksize;
    }
}


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs
(natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    *e = 0;

    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        grad->ptr.p_double[i] = 0;
    }
    *e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpbase_mlpchunkedgradient(network, xy, i, ae_minint(ssize, i+mlpbase_chunksize, _state)-i, e, grad, ae_true, _state);
        i = i+mlpbase_chunksize;
    }
}


/*************************************************************************
Batch Hessian calculation (natural error function) using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.
     
     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessiannbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{

    *e = 0;

    mlpbase_mlphessianbatchinternal(network, xy, ssize, ae_true, e, grad, h, _state);
}


/*************************************************************************
Batch Hessian calculation using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessianbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{

    *e = 0;

    mlpbase_mlphessianbatchinternal(network, xy, ssize, ae_false, e, grad, h, _state);
}


/*************************************************************************
Internal subroutine, shouldn't be called by user.
*************************************************************************/
void mlpinternalprocessvector(/* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* dfdnet,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t istart;
    ae_int_t offs;
    double net;
    double f;
    double df;
    double d2f;
    double mx;
    ae_bool perr;


    
    /*
     * Read network geometry
     */
    nin = structinfo->ptr.p_int[1];
    nout = structinfo->ptr.p_int[2];
    ntotal = structinfo->ptr.p_int[3];
    istart = structinfo->ptr.p_int[5];
    
    /*
     * Inputs standartisation and putting in the network
     */
    for(i=0; i<=nin-1; i++)
    {
        if( ae_fp_neq(columnsigmas->ptr.p_double[i],0) )
        {
            neurons->ptr.p_double[i] = (x->ptr.p_double[i]-columnmeans->ptr.p_double[i])/columnsigmas->ptr.p_double[i];
        }
        else
        {
            neurons->ptr.p_double[i] = x->ptr.p_double[i]-columnmeans->ptr.p_double[i];
        }
    }
    
    /*
     * Process network
     */
    for(i=0; i<=ntotal-1; i++)
    {
        offs = istart+i*mlpbase_nfieldwidth;
        if( structinfo->ptr.p_int[offs+0]>0 )
        {
            
            /*
             * Activation function
             */
            mlpbase_mlpactivationfunction(neurons->ptr.p_double[structinfo->ptr.p_int[offs+2]], structinfo->ptr.p_int[offs+0], &f, &df, &d2f, _state);
            neurons->ptr.p_double[i] = f;
            dfdnet->ptr.p_double[i] = df;
        }
        if( structinfo->ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator
             */
            n1 = structinfo->ptr.p_int[offs+2];
            n2 = n1+structinfo->ptr.p_int[offs+1]-1;
            w1 = structinfo->ptr.p_int[offs+3];
            w2 = w1+structinfo->ptr.p_int[offs+1]-1;
            net = ae_v_dotproduct(&weights->ptr.p_double[w1], 1, &neurons->ptr.p_double[n1], 1, ae_v_len(w1,w2));
            neurons->ptr.p_double[i] = net;
            dfdnet->ptr.p_double[i] = 1.0;
        }
        if( structinfo->ptr.p_int[offs+0]<0 )
        {
            perr = ae_true;
            if( structinfo->ptr.p_int[offs+0]==-2 )
            {
                
                /*
                 * input neuron, left unchanged
                 */
                perr = ae_false;
            }
            if( structinfo->ptr.p_int[offs+0]==-3 )
            {
                
                /*
                 * "-1" neuron
                 */
                neurons->ptr.p_double[i] = -1;
                perr = ae_false;
            }
            if( structinfo->ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * "0" neuron
                 */
                neurons->ptr.p_double[i] = 0;
                perr = ae_false;
            }
            ae_assert(!perr, "MLPInternalProcessVector: internal error - unknown neuron type!", _state);
        }
    }
    
    /*
     * Extract result
     */
    ae_v_move(&y->ptr.p_double[0], 1, &neurons->ptr.p_double[ntotal-nout], 1, ae_v_len(0,nout-1));
    
    /*
     * Softmax post-processing or standardisation if needed
     */
    ae_assert(structinfo->ptr.p_int[6]==0||structinfo->ptr.p_int[6]==1, "MLPInternalProcessVector: unknown normalization type!", _state);
    if( structinfo->ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax
         */
        mx = y->ptr.p_double[0];
        for(i=1; i<=nout-1; i++)
        {
            mx = ae_maxreal(mx, y->ptr.p_double[i], _state);
        }
        net = 0;
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = ae_exp(y->ptr.p_double[i]-mx, _state);
            net = net+y->ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = y->ptr.p_double[i]/net;
        }
    }
    else
    {
        
        /*
         * Standardisation
         */
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = y->ptr.p_double[i]*columnsigmas->ptr.p_double[nin+i]+columnmeans->ptr.p_double[nin+i];
        }
    }
}


/*************************************************************************
Internal subroutine: adding new input layer to network
*************************************************************************/
static void mlpbase_addinputlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state)
{


    lsizes->ptr.p_int[0] = ncount;
    ltypes->ptr.p_int[0] = -2;
    lconnfirst->ptr.p_int[0] = 0;
    lconnlast->ptr.p_int[0] = 0;
    *lastproc = 0;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void mlpbase_addbiasedsummatorlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state)
{


    lsizes->ptr.p_int[*lastproc+1] = 1;
    ltypes->ptr.p_int[*lastproc+1] = -3;
    lconnfirst->ptr.p_int[*lastproc+1] = 0;
    lconnlast->ptr.p_int[*lastproc+1] = 0;
    lsizes->ptr.p_int[*lastproc+2] = ncount;
    ltypes->ptr.p_int[*lastproc+2] = 0;
    lconnfirst->ptr.p_int[*lastproc+2] = *lastproc;
    lconnlast->ptr.p_int[*lastproc+2] = *lastproc+1;
    *lastproc = *lastproc+2;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void mlpbase_addactivationlayer(ae_int_t functype,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state)
{


    ae_assert(functype>0, "AddActivationLayer: incorrect function type", _state);
    lsizes->ptr.p_int[*lastproc+1] = lsizes->ptr.p_int[*lastproc];
    ltypes->ptr.p_int[*lastproc+1] = functype;
    lconnfirst->ptr.p_int[*lastproc+1] = *lastproc;
    lconnlast->ptr.p_int[*lastproc+1] = *lastproc;
    *lastproc = *lastproc+1;
}


/*************************************************************************
Internal subroutine: adding new zero layer to network
*************************************************************************/
static void mlpbase_addzerolayer(/* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state)
{


    lsizes->ptr.p_int[*lastproc+1] = 1;
    ltypes->ptr.p_int[*lastproc+1] = -4;
    lconnfirst->ptr.p_int[*lastproc+1] = 0;
    lconnlast->ptr.p_int[*lastproc+1] = 0;
    *lastproc = *lastproc+1;
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
static void mlpbase_mlpcreate(ae_int_t nin,
     ae_int_t nout,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t layerscount,
     ae_bool isclsnet,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t wcount;
    ae_int_t offs;
    ae_int_t nprocessed;
    ae_int_t wallocated;
    ae_vector localtemp;
    ae_vector lnfirst;
    ae_vector lnsyn;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&localtemp, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lnsyn, 0, DT_INT, _state, ae_true);

    
    /*
     * Check
     */
    ae_assert(layerscount>0, "MLPCreate: wrong parameters!", _state);
    ae_assert(ltypes->ptr.p_int[0]==-2, "MLPCreate: wrong LTypes[0] (must be -2)!", _state);
    for(i=0; i<=layerscount-1; i++)
    {
        ae_assert(lsizes->ptr.p_int[i]>0, "MLPCreate: wrong LSizes!", _state);
        ae_assert(lconnfirst->ptr.p_int[i]>=0&&(lconnfirst->ptr.p_int[i]<i||i==0), "MLPCreate: wrong LConnFirst!", _state);
        ae_assert(lconnlast->ptr.p_int[i]>=lconnfirst->ptr.p_int[i]&&(lconnlast->ptr.p_int[i]<i||i==0), "MLPCreate: wrong LConnLast!", _state);
    }
    
    /*
     * Build network geometry
     */
    ae_vector_set_length(&lnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lnsyn, layerscount-1+1, _state);
    ntotal = 0;
    wcount = 0;
    for(i=0; i<=layerscount-1; i++)
    {
        
        /*
         * Analyze connections.
         * This code must throw an assertion in case of unknown LTypes[I]
         */
        lnsyn.ptr.p_int[i] = -1;
        if( ltypes->ptr.p_int[i]>=0 )
        {
            lnsyn.ptr.p_int[i] = 0;
            for(j=lconnfirst->ptr.p_int[i]; j<=lconnlast->ptr.p_int[i]; j++)
            {
                lnsyn.ptr.p_int[i] = lnsyn.ptr.p_int[i]+lsizes->ptr.p_int[j];
            }
        }
        else
        {
            if( (ltypes->ptr.p_int[i]==-2||ltypes->ptr.p_int[i]==-3)||ltypes->ptr.p_int[i]==-4 )
            {
                lnsyn.ptr.p_int[i] = 0;
            }
        }
        ae_assert(lnsyn.ptr.p_int[i]>=0, "MLPCreate: internal error #0!", _state);
        
        /*
         * Other info
         */
        lnfirst.ptr.p_int[i] = ntotal;
        ntotal = ntotal+lsizes->ptr.p_int[i];
        if( ltypes->ptr.p_int[i]==0 )
        {
            wcount = wcount+lnsyn.ptr.p_int[i]*lsizes->ptr.p_int[i];
        }
    }
    ssize = 7+ntotal*mlpbase_nfieldwidth;
    
    /*
     * Allocate
     */
    ae_vector_set_length(&network->structinfo, ssize-1+1, _state);
    ae_vector_set_length(&network->weights, wcount-1+1, _state);
    if( isclsnet )
    {
        ae_vector_set_length(&network->columnmeans, nin-1+1, _state);
        ae_vector_set_length(&network->columnsigmas, nin-1+1, _state);
    }
    else
    {
        ae_vector_set_length(&network->columnmeans, nin+nout-1+1, _state);
        ae_vector_set_length(&network->columnsigmas, nin+nout-1+1, _state);
    }
    ae_vector_set_length(&network->neurons, ntotal-1+1, _state);
    ae_matrix_set_length(&network->chunks, 3*ntotal+1, mlpbase_chunksize-1+1, _state);
    ae_vector_set_length(&network->nwbuf, ae_maxint(wcount, 2*nout, _state)-1+1, _state);
    ae_vector_set_length(&network->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&network->x, nin-1+1, _state);
    ae_vector_set_length(&network->y, nout-1+1, _state);
    ae_vector_set_length(&network->derror, ntotal-1+1, _state);
    
    /*
     * Fill structure: global info
     */
    network->structinfo.ptr.p_int[0] = ssize;
    network->structinfo.ptr.p_int[1] = nin;
    network->structinfo.ptr.p_int[2] = nout;
    network->structinfo.ptr.p_int[3] = ntotal;
    network->structinfo.ptr.p_int[4] = wcount;
    network->structinfo.ptr.p_int[5] = 7;
    if( isclsnet )
    {
        network->structinfo.ptr.p_int[6] = 1;
    }
    else
    {
        network->structinfo.ptr.p_int[6] = 0;
    }
    
    /*
     * Fill structure: neuron connections
     */
    nprocessed = 0;
    wallocated = 0;
    for(i=0; i<=layerscount-1; i++)
    {
        for(j=0; j<=lsizes->ptr.p_int[i]-1; j++)
        {
            offs = network->structinfo.ptr.p_int[5]+nprocessed*mlpbase_nfieldwidth;
            network->structinfo.ptr.p_int[offs+0] = ltypes->ptr.p_int[i];
            if( ltypes->ptr.p_int[i]==0 )
            {
                
                /*
                 * Adaptive summator:
                 * * connections with weights to previous neurons
                 */
                network->structinfo.ptr.p_int[offs+1] = lnsyn.ptr.p_int[i];
                network->structinfo.ptr.p_int[offs+2] = lnfirst.ptr.p_int[lconnfirst->ptr.p_int[i]];
                network->structinfo.ptr.p_int[offs+3] = wallocated;
                wallocated = wallocated+lnsyn.ptr.p_int[i];
                nprocessed = nprocessed+1;
            }
            if( ltypes->ptr.p_int[i]>0 )
            {
                
                /*
                 * Activation layer:
                 * * each neuron connected to one (only one) of previous neurons.
                 * * no weights
                 */
                network->structinfo.ptr.p_int[offs+1] = 1;
                network->structinfo.ptr.p_int[offs+2] = lnfirst.ptr.p_int[lconnfirst->ptr.p_int[i]]+j;
                network->structinfo.ptr.p_int[offs+3] = -1;
                nprocessed = nprocessed+1;
            }
            if( (ltypes->ptr.p_int[i]==-2||ltypes->ptr.p_int[i]==-3)||ltypes->ptr.p_int[i]==-4 )
            {
                nprocessed = nprocessed+1;
            }
        }
    }
    ae_assert(wallocated==wcount, "MLPCreate: internal error #1!", _state);
    ae_assert(nprocessed==ntotal, "MLPCreate: internal error #2!", _state);
    
    /*
     * Fill weights by small random values
     * Initialize means and sigmas
     */
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0;
        network->columnsigmas.ptr.p_double[i] = 1;
    }
    if( !isclsnet )
    {
        for(i=0; i<=nout-1; i++)
        {
            network->columnmeans.ptr.p_double[nin+i] = 0;
            network->columnsigmas.ptr.p_double[nin+i] = 1;
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
static void mlpbase_mlpactivationfunction(double net,
     ae_int_t k,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state)
{
    double net2;
    double arg;
    double root;
    double r;

    *f = 0;
    *df = 0;
    *d2f = 0;

    *f = 0;
    *df = 0;
    if( k==1 )
    {
        
        /*
         * TanH activation function
         */
        if( ae_fp_less(ae_fabs(net, _state),100) )
        {
            *f = ae_tanh(net, _state);
        }
        else
        {
            *f = ae_sign(net, _state);
        }
        *df = 1-ae_sqr(*f, _state);
        *d2f = -2*(*f)*(*df);
        return;
    }
    if( k==3 )
    {
        
        /*
         * EX activation function
         */
        if( ae_fp_greater_eq(net,0) )
        {
            net2 = net*net;
            arg = net2+1;
            root = ae_sqrt(arg, _state);
            *f = net+root;
            r = net/root;
            *df = 1+r;
            *d2f = (root-net*r)/arg;
        }
        else
        {
            *f = ae_exp(net, _state);
            *df = *f;
            *d2f = *f;
        }
        return;
    }
    if( k==2 )
    {
        *f = ae_exp(-ae_sqr(net, _state), _state);
        *df = -2*net*(*f);
        *d2f = -2*(*f+*df*net);
        return;
    }
}


/*************************************************************************
Internal subroutine for Hessian calculation.

WARNING!!! Unspeakable math far beyong human capabilities :)
*************************************************************************/
static void mlpbase_mlphessianbatchinternal(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_bool naturalerr,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t kl;
    ae_int_t offs;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    double s;
    double t;
    double v;
    double et;
    ae_bool bflag;
    double f;
    double df;
    double d2f;
    double deidyj;
    double mx;
    double q;
    double z;
    double s2;
    double expi;
    double expj;
    ae_vector x;
    ae_vector desiredy;
    ae_vector gt;
    ae_vector zeros;
    ae_matrix rx;
    ae_matrix ry;
    ae_matrix rdx;
    ae_matrix rdy;

    ae_frame_make(_state, &_frame_block);
    *e = 0;
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&desiredy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&gt, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&zeros, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ry, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rdx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rdy, 0, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Prepare
     */
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&desiredy, nout-1+1, _state);
    ae_vector_set_length(&zeros, wcount-1+1, _state);
    ae_vector_set_length(&gt, wcount-1+1, _state);
    ae_matrix_set_length(&rx, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&ry, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&rdx, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&rdy, ntotal+nout-1+1, wcount-1+1, _state);
    *e = 0;
    for(i=0; i<=wcount-1; i++)
    {
        zeros.ptr.p_double[i] = 0;
    }
    ae_v_move(&grad->ptr.p_double[0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    for(i=0; i<=wcount-1; i++)
    {
        ae_v_move(&h->ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    }
    
    /*
     * Process
     */
    for(k=0; k<=ssize-1; k++)
    {
        
        /*
         * Process vector with MLPGradN.
         * Now Neurons, DFDNET and DError contains results of the last run.
         */
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nin-1));
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels outputs
             */
            kl = ae_round(xy->ptr.pp_double[k][nin], _state);
            for(i=0; i<=nout-1; i++)
            {
                if( i==kl )
                {
                    desiredy.ptr.p_double[i] = 1;
                }
                else
                {
                    desiredy.ptr.p_double[i] = 0;
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            ae_v_move(&desiredy.ptr.p_double[0], 1, &xy->ptr.pp_double[k][nin], 1, ae_v_len(0,nout-1));
        }
        if( naturalerr )
        {
            mlpgradn(network, &x, &desiredy, &et, &gt, _state);
        }
        else
        {
            mlpgrad(network, &x, &desiredy, &et, &gt, _state);
        }
        
        /*
         * grad, error
         */
        *e = *e+et;
        ae_v_add(&grad->ptr.p_double[0], 1, &gt.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Hessian.
         * Forward pass of the R-algorithm
         */
        for(i=0; i<=ntotal-1; i++)
        {
            offs = istart+i*mlpbase_nfieldwidth;
            ae_v_move(&rx.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            ae_v_move(&ry.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            if( network->structinfo.ptr.p_int[offs+0]>0 )
            {
                
                /*
                 * Activation function
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                ae_v_move(&rx.ptr.pp_double[i][0], 1, &ry.ptr.pp_double[n1][0], 1, ae_v_len(0,wcount-1));
                v = network->dfdnet.ptr.p_double[i];
                ae_v_moved(&ry.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
            }
            if( network->structinfo.ptr.p_int[offs+0]==0 )
            {
                
                /*
                 * Adaptive summator
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
                w1 = network->structinfo.ptr.p_int[offs+3];
                w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
                for(j=n1; j<=n2; j++)
                {
                    v = network->weights.ptr.p_double[w1+j-n1];
                    ae_v_addd(&rx.ptr.pp_double[i][0], 1, &ry.ptr.pp_double[j][0], 1, ae_v_len(0,wcount-1), v);
                    rx.ptr.pp_double[i][w1+j-n1] = rx.ptr.pp_double[i][w1+j-n1]+network->neurons.ptr.p_double[j];
                }
                ae_v_move(&ry.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
            }
            if( network->structinfo.ptr.p_int[offs+0]<0 )
            {
                bflag = ae_true;
                if( network->structinfo.ptr.p_int[offs+0]==-2 )
                {
                    
                    /*
                     * input neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                if( network->structinfo.ptr.p_int[offs+0]==-3 )
                {
                    
                    /*
                     * "-1" neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                if( network->structinfo.ptr.p_int[offs+0]==-4 )
                {
                    
                    /*
                     * "0" neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                ae_assert(!bflag, "MLPHessianNBatch: internal error - unknown neuron type!", _state);
            }
        }
        
        /*
         * Hessian. Backward pass of the R-algorithm.
         *
         * Stage 1. Initialize RDY
         */
        for(i=0; i<=ntotal+nout-1; i++)
        {
            ae_v_move(&rdy.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
        if( network->structinfo.ptr.p_int[6]==0 )
        {
            
            /*
             * Standardisation.
             *
             * In context of the Hessian calculation standardisation
             * is considered as additional layer with weightless
             * activation function:
             *
             * F(NET) := Sigma*NET
             *
             * So we add one more layer to forward pass, and
             * make forward/backward pass through this layer.
             */
            for(i=0; i<=nout-1; i++)
            {
                n1 = ntotal-nout+i;
                n2 = ntotal+i;
                
                /*
                 * Forward pass from N1 to N2
                 */
                ae_v_move(&rx.ptr.pp_double[n2][0], 1, &ry.ptr.pp_double[n1][0], 1, ae_v_len(0,wcount-1));
                v = network->columnsigmas.ptr.p_double[nin+i];
                ae_v_moved(&ry.ptr.pp_double[n2][0], 1, &rx.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1), v);
                
                /*
                 * Initialization of RDY
                 */
                ae_v_move(&rdy.ptr.pp_double[n2][0], 1, &ry.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1));
                
                /*
                 * Backward pass from N2 to N1:
                 * 1. Calculate R(dE/dX).
                 * 2. No R(dE/dWij) is needed since weight of activation neuron
                 *    is fixed to 1. So we can update R(dE/dY) for
                 *    the connected neuron (note that Vij=0, Wij=1)
                 */
                df = network->columnsigmas.ptr.p_double[nin+i];
                ae_v_moved(&rdx.ptr.pp_double[n2][0], 1, &rdy.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1), df);
                ae_v_add(&rdy.ptr.pp_double[n1][0], 1, &rdx.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1));
            }
        }
        else
        {
            
            /*
             * Softmax.
             *
             * Initialize RDY using generalized expression for ei'(yi)
             * (see expression (9) from p. 5 of "Fast Exact Multiplication by the Hessian").
             *
             * When we are working with softmax network, generalized
             * expression for ei'(yi) is used because softmax
             * normalization leads to ei, which depends on all y's
             */
            if( naturalerr )
            {
                
                /*
                 * softmax + cross-entropy.
                 * We have:
                 *
                 * S = sum(exp(yk)),
                 * ei = sum(trn)*exp(yi)/S-trn_i
                 *
                 * j=i:   d(ei)/d(yj) = T*exp(yi)*(S-exp(yi))/S^2
                 * j<>i:  d(ei)/d(yj) = -T*exp(yi)*exp(yj)/S^2
                 */
                t = 0;
                for(i=0; i<=nout-1; i++)
                {
                    t = t+desiredy.ptr.p_double[i];
                }
                mx = network->neurons.ptr.p_double[ntotal-nout];
                for(i=0; i<=nout-1; i++)
                {
                    mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
                }
                s = 0;
                for(i=0; i<=nout-1; i++)
                {
                    network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                    s = s+network->nwbuf.ptr.p_double[i];
                }
                for(i=0; i<=nout-1; i++)
                {
                    for(j=0; j<=nout-1; j++)
                    {
                        if( j==i )
                        {
                            deidyj = t*network->nwbuf.ptr.p_double[i]*(s-network->nwbuf.ptr.p_double[i])/ae_sqr(s, _state);
                            ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+i][0], 1, ae_v_len(0,wcount-1), deidyj);
                        }
                        else
                        {
                            deidyj = -t*network->nwbuf.ptr.p_double[i]*network->nwbuf.ptr.p_double[j]/ae_sqr(s, _state);
                            ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+j][0], 1, ae_v_len(0,wcount-1), deidyj);
                        }
                    }
                }
            }
            else
            {
                
                /*
                 * For a softmax + squared error we have expression
                 * far beyond human imagination so we dont even try
                 * to comment on it. Just enjoy the code...
                 *
                 * P.S. That's why "natural error" is called "natural" -
                 * compact beatiful expressions, fast code....
                 */
                mx = network->neurons.ptr.p_double[ntotal-nout];
                for(i=0; i<=nout-1; i++)
                {
                    mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
                }
                s = 0;
                s2 = 0;
                for(i=0; i<=nout-1; i++)
                {
                    network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                    s = s+network->nwbuf.ptr.p_double[i];
                    s2 = s2+ae_sqr(network->nwbuf.ptr.p_double[i], _state);
                }
                q = 0;
                for(i=0; i<=nout-1; i++)
                {
                    q = q+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])*network->nwbuf.ptr.p_double[i];
                }
                for(i=0; i<=nout-1; i++)
                {
                    z = -q+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])*s;
                    expi = network->nwbuf.ptr.p_double[i];
                    for(j=0; j<=nout-1; j++)
                    {
                        expj = network->nwbuf.ptr.p_double[j];
                        if( j==i )
                        {
                            deidyj = expi/ae_sqr(s, _state)*((z+expi)*(s-2*expi)/s+expi*s2/ae_sqr(s, _state));
                        }
                        else
                        {
                            deidyj = expi*expj/ae_sqr(s, _state)*(s2/ae_sqr(s, _state)-2*z/s-(expi+expj)/s+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])-(network->y.ptr.p_double[j]-desiredy.ptr.p_double[j]));
                        }
                        ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+j][0], 1, ae_v_len(0,wcount-1), deidyj);
                    }
                }
            }
        }
        
        /*
         * Hessian. Backward pass of the R-algorithm
         *
         * Stage 2. Process.
         */
        for(i=ntotal-1; i>=0; i--)
        {
            
            /*
             * Possible variants:
             * 1. Activation function
             * 2. Adaptive summator
             * 3. Special neuron
             */
            offs = istart+i*mlpbase_nfieldwidth;
            if( network->structinfo.ptr.p_int[offs+0]>0 )
            {
                n1 = network->structinfo.ptr.p_int[offs+2];
                
                /*
                 * First, calculate R(dE/dX).
                 */
                mlpbase_mlpactivationfunction(network->neurons.ptr.p_double[n1], network->structinfo.ptr.p_int[offs+0], &f, &df, &d2f, _state);
                v = d2f*network->derror.ptr.p_double[i];
                ae_v_moved(&rdx.ptr.pp_double[i][0], 1, &rdy.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), df);
                ae_v_addd(&rdx.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                
                /*
                 * No R(dE/dWij) is needed since weight of activation neuron
                 * is fixed to 1.
                 *
                 * So we can update R(dE/dY) for the connected neuron.
                 * (note that Vij=0, Wij=1)
                 */
                ae_v_add(&rdy.ptr.pp_double[n1][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
            }
            if( network->structinfo.ptr.p_int[offs+0]==0 )
            {
                
                /*
                 * Adaptive summator
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
                w1 = network->structinfo.ptr.p_int[offs+3];
                w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
                
                /*
                 * First, calculate R(dE/dX).
                 */
                ae_v_move(&rdx.ptr.pp_double[i][0], 1, &rdy.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                
                /*
                 * Then, calculate R(dE/dWij)
                 */
                for(j=w1; j<=w2; j++)
                {
                    v = network->neurons.ptr.p_double[n1+j-w1];
                    ae_v_addd(&h->ptr.pp_double[j][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                    v = network->derror.ptr.p_double[i];
                    ae_v_addd(&h->ptr.pp_double[j][0], 1, &ry.ptr.pp_double[n1+j-w1][0], 1, ae_v_len(0,wcount-1), v);
                }
                
                /*
                 * And finally, update R(dE/dY) for connected neurons.
                 */
                for(j=w1; j<=w2; j++)
                {
                    v = network->weights.ptr.p_double[j];
                    ae_v_addd(&rdy.ptr.pp_double[n1+j-w1][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                    rdy.ptr.pp_double[n1+j-w1][j] = rdy.ptr.pp_double[n1+j-w1][j]+network->derror.ptr.p_double[i];
                }
            }
            if( network->structinfo.ptr.p_int[offs+0]<0 )
            {
                bflag = ae_false;
                if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
                {
                    
                    /*
                     * Special neuron type, no back-propagation required
                     */
                    bflag = ae_true;
                }
                ae_assert(bflag, "MLPHessianNBatch: unknown neuron type!", _state);
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine

Network must be processed by MLPProcess on X
*************************************************************************/
static void mlpbase_mlpinternalcalculategradient(multilayerperceptron* network,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* derror,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t offs;
    double dedf;
    double dfdnet;
    double v;
    double fown;
    double deown;
    double net;
    double mx;
    ae_bool bflag;


    
    /*
     * Read network geometry
     */
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Pre-processing of dError/dOut:
     * from dError/dOut(normalized) to dError/dOut(non-normalized)
     */
    ae_assert(network->structinfo.ptr.p_int[6]==0||network->structinfo.ptr.p_int[6]==1, "MLPInternalCalculateGradient: unknown normalization type!", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax
         */
        if( !naturalerrorfunc )
        {
            mx = network->neurons.ptr.p_double[ntotal-nout];
            for(i=0; i<=nout-1; i++)
            {
                mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
            }
            net = 0;
            for(i=0; i<=nout-1; i++)
            {
                network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                net = net+network->nwbuf.ptr.p_double[i];
            }
            v = ae_v_dotproduct(&network->derror.ptr.p_double[ntotal-nout], 1, &network->nwbuf.ptr.p_double[0], 1, ae_v_len(ntotal-nout,ntotal-1));
            for(i=0; i<=nout-1; i++)
            {
                fown = network->nwbuf.ptr.p_double[i];
                deown = network->derror.ptr.p_double[ntotal-nout+i];
                network->nwbuf.ptr.p_double[nout+i] = (-v+deown*fown+deown*(net-fown))*fown/ae_sqr(net, _state);
            }
            for(i=0; i<=nout-1; i++)
            {
                network->derror.ptr.p_double[ntotal-nout+i] = network->nwbuf.ptr.p_double[nout+i];
            }
        }
    }
    else
    {
        
        /*
         * Un-standardisation
         */
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = network->derror.ptr.p_double[ntotal-nout+i]*network->columnsigmas.ptr.p_double[nin+i];
        }
    }
    
    /*
     * Backpropagation
     */
    for(i=ntotal-1; i>=0; i--)
    {
        
        /*
         * Extract info
         */
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0 )
        {
            
            /*
             * Activation function
             */
            dedf = network->derror.ptr.p_double[i];
            dfdnet = network->dfdnet.ptr.p_double[i];
            derror->ptr.p_double[network->structinfo.ptr.p_int[offs+2]] = derror->ptr.p_double[network->structinfo.ptr.p_int[offs+2]]+dedf*dfdnet;
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            dedf = network->derror.ptr.p_double[i];
            dfdnet = 1.0;
            v = dedf*dfdnet;
            ae_v_moved(&grad->ptr.p_double[w1], 1, &neurons->ptr.p_double[n1], 1, ae_v_len(w1,w2), v);
            ae_v_addd(&derror->ptr.p_double[n1], 1, &weights->ptr.p_double[w1], 1, ae_v_len(n1,n2), v);
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * Special neuron type, no back-propagation required
                 */
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPInternalCalculateGradient: unknown neuron type!", _state);
        }
    }
}


/*************************************************************************
Internal subroutine, chunked gradient
*************************************************************************/
static void mlpbase_mlpchunkedgradient(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t kl;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t c1;
    ae_int_t c2;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t offs;
    double f;
    double df;
    double d2f;
    double v;
    double s;
    double fown;
    double deown;
    double net;
    double lnnet;
    double mx;
    ae_bool bflag;
    ae_int_t istart;
    ae_int_t ineurons;
    ae_int_t idfdnet;
    ae_int_t iderror;
    ae_int_t izeros;


    
    /*
     * Read network geometry, prepare data
     */
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    c1 = cstart;
    c2 = cstart+csize-1;
    ineurons = 0;
    idfdnet = ntotal;
    iderror = 2*ntotal;
    izeros = 3*ntotal;
    for(j=0; j<=csize-1; j++)
    {
        network->chunks.ptr.pp_double[izeros][j] = 0;
    }
    
    /*
     * Forward pass:
     * 1. Load inputs from XY to Chunks[0:NIn-1,0:CSize-1]
     * 2. Forward pass
     */
    for(i=0; i<=nin-1; i++)
    {
        for(j=0; j<=csize-1; j++)
        {
            if( ae_fp_neq(network->columnsigmas.ptr.p_double[i],0) )
            {
                network->chunks.ptr.pp_double[i][j] = (xy->ptr.pp_double[c1+j][i]-network->columnmeans.ptr.p_double[i])/network->columnsigmas.ptr.p_double[i];
            }
            else
            {
                network->chunks.ptr.pp_double[i][j] = xy->ptr.pp_double[c1+j][i]-network->columnmeans.ptr.p_double[i];
            }
        }
    }
    for(i=0; i<=ntotal-1; i++)
    {
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0 )
        {
            
            /*
             * Activation function:
             * * calculate F vector, F(i) = F(NET(i))
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            ae_v_move(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[n1][0], 1, ae_v_len(0,csize-1));
            for(j=0; j<=csize-1; j++)
            {
                mlpbase_mlpactivationfunction(network->chunks.ptr.pp_double[i][j], network->structinfo.ptr.p_int[offs+0], &f, &df, &d2f, _state);
                network->chunks.ptr.pp_double[i][j] = f;
                network->chunks.ptr.pp_double[idfdnet+i][j] = df;
            }
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator:
             * * calculate NET vector, NET(i) = SUM(W(j,i)*Neurons(j),j=N1..N2)
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            ae_v_move(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[izeros][0], 1, ae_v_len(0,csize-1));
            for(j=n1; j<=n2; j++)
            {
                v = network->weights.ptr.p_double[w1+j-n1];
                ae_v_addd(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[j][0], 1, ae_v_len(0,csize-1), v);
            }
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( network->structinfo.ptr.p_int[offs+0]==-2 )
            {
                
                /*
                 * input neuron, left unchanged
                 */
                bflag = ae_true;
            }
            if( network->structinfo.ptr.p_int[offs+0]==-3 )
            {
                
                /*
                 * "-1" neuron
                 */
                for(k=0; k<=csize-1; k++)
                {
                    network->chunks.ptr.pp_double[i][k] = -1;
                }
                bflag = ae_true;
            }
            if( network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * "0" neuron
                 */
                for(k=0; k<=csize-1; k++)
                {
                    network->chunks.ptr.pp_double[i][k] = 0;
                }
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPChunkedGradient: internal error - unknown neuron type!", _state);
        }
    }
    
    /*
     * Post-processing, error, dError/dOut
     */
    for(i=0; i<=ntotal-1; i++)
    {
        ae_v_move(&network->chunks.ptr.pp_double[iderror+i][0], 1, &network->chunks.ptr.pp_double[izeros][0], 1, ae_v_len(0,csize-1));
    }
    ae_assert(network->structinfo.ptr.p_int[6]==0||network->structinfo.ptr.p_int[6]==1, "MLPChunkedGradient: unknown normalization type!", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax output, classification network.
         *
         * For each K = 0..CSize-1 do:
         * 1. place exp(outputs[k]) to NWBuf[0:NOut-1]
         * 2. place sum(exp(..)) to NET
         * 3. calculate dError/dOut and place it to the second block of Chunks
         */
        for(k=0; k<=csize-1; k++)
        {
            
            /*
             * Normalize
             */
            mx = network->chunks.ptr.pp_double[ntotal-nout][k];
            for(i=1; i<=nout-1; i++)
            {
                mx = ae_maxreal(mx, network->chunks.ptr.pp_double[ntotal-nout+i][k], _state);
            }
            net = 0;
            for(i=0; i<=nout-1; i++)
            {
                network->nwbuf.ptr.p_double[i] = ae_exp(network->chunks.ptr.pp_double[ntotal-nout+i][k]-mx, _state);
                net = net+network->nwbuf.ptr.p_double[i];
            }
            
            /*
             * Calculate error function and dError/dOut
             */
            if( naturalerrorfunc )
            {
                
                /*
                 * Natural error func.
                 *
                 */
                s = 1;
                lnnet = ae_log(net, _state);
                kl = ae_round(xy->ptr.pp_double[cstart+k][nin], _state);
                for(i=0; i<=nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = 1;
                    }
                    else
                    {
                        v = 0;
                    }
                    network->chunks.ptr.pp_double[iderror+ntotal-nout+i][k] = s*network->nwbuf.ptr.p_double[i]/net-v;
                    *e = *e+mlpbase_safecrossentropy(v, network->nwbuf.ptr.p_double[i]/net, _state);
                }
            }
            else
            {
                
                /*
                 * Least squares error func
                 * Error, dError/dOut(normalized)
                 */
                kl = ae_round(xy->ptr.pp_double[cstart+k][nin], _state);
                for(i=0; i<=nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = network->nwbuf.ptr.p_double[i]/net-1;
                    }
                    else
                    {
                        v = network->nwbuf.ptr.p_double[i]/net;
                    }
                    network->nwbuf.ptr.p_double[nout+i] = v;
                    *e = *e+ae_sqr(v, _state)/2;
                }
                
                /*
                 * From dError/dOut(normalized) to dError/dOut(non-normalized)
                 */
                v = ae_v_dotproduct(&network->nwbuf.ptr.p_double[nout], 1, &network->nwbuf.ptr.p_double[0], 1, ae_v_len(nout,2*nout-1));
                for(i=0; i<=nout-1; i++)
                {
                    fown = network->nwbuf.ptr.p_double[i];
                    deown = network->nwbuf.ptr.p_double[nout+i];
                    network->chunks.ptr.pp_double[iderror+ntotal-nout+i][k] = (-v+deown*fown+deown*(net-fown))*fown/ae_sqr(net, _state);
                }
            }
        }
    }
    else
    {
        
        /*
         * Normal output, regression network
         *
         * For each K = 0..CSize-1 do:
         * 1. calculate dError/dOut and place it to the second block of Chunks
         */
        for(i=0; i<=nout-1; i++)
        {
            for(j=0; j<=csize-1; j++)
            {
                v = network->chunks.ptr.pp_double[ntotal-nout+i][j]*network->columnsigmas.ptr.p_double[nin+i]+network->columnmeans.ptr.p_double[nin+i]-xy->ptr.pp_double[cstart+j][nin+i];
                network->chunks.ptr.pp_double[iderror+ntotal-nout+i][j] = v*network->columnsigmas.ptr.p_double[nin+i];
                *e = *e+ae_sqr(v, _state)/2;
            }
        }
    }
    
    /*
     * Backpropagation
     */
    for(i=ntotal-1; i>=0; i--)
    {
        
        /*
         * Extract info
         */
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0 )
        {
            
            /*
             * Activation function
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            for(k=0; k<=csize-1; k++)
            {
                network->chunks.ptr.pp_double[iderror+i][k] = network->chunks.ptr.pp_double[iderror+i][k]*network->chunks.ptr.pp_double[idfdnet+i][k];
            }
            ae_v_add(&network->chunks.ptr.pp_double[iderror+n1][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1));
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * "Normal" activation function
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            for(j=w1; j<=w2; j++)
            {
                v = ae_v_dotproduct(&network->chunks.ptr.pp_double[n1+j-w1][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1));
                grad->ptr.p_double[j] = grad->ptr.p_double[j]+v;
            }
            for(j=n1; j<=n2; j++)
            {
                v = network->weights.ptr.p_double[w1+j-n1];
                ae_v_addd(&network->chunks.ptr.pp_double[iderror+j][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1), v);
            }
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * Special neuron type, no back-propagation required
                 */
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPInternalCalculateGradient: unknown neuron type!", _state);
        }
    }
}


/*************************************************************************
Returns T*Ln(T/Z), guarded against overflow/underflow.
Internal subroutine.
*************************************************************************/
static double mlpbase_safecrossentropy(double t,
     double z,
     ae_state *_state)
{
    double r;
    double result;


    if( ae_fp_eq(t,0) )
    {
        result = 0;
    }
    else
    {
        if( ae_fp_greater(ae_fabs(z, _state),1) )
        {
            
            /*
             * Shouldn't be the case with softmax,
             * but we just want to be sure.
             */
            if( ae_fp_eq(t/z,0) )
            {
                r = ae_minrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        else
        {
            
            /*
             * Normal case
             */
            if( ae_fp_eq(z,0)||ae_fp_greater_eq(ae_fabs(t, _state),ae_maxrealnumber*ae_fabs(z, _state)) )
            {
                r = ae_maxrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        result = t*ae_log(r, _state);
    }
    return result;
}


ae_bool _multilayerperceptron_init(multilayerperceptron* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->structinfo, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->weights, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnmeans, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnsigmas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->neurons, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dfdnet, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->derror, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->chunks, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->nwbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _multilayerperceptron_init_copy(multilayerperceptron* dst, multilayerperceptron* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->structinfo, &src->structinfo, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->weights, &src->weights, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnmeans, &src->columnmeans, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnsigmas, &src->columnsigmas, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->neurons, &src->neurons, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dfdnet, &src->dfdnet, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->derror, &src->derror, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->chunks, &src->chunks, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->nwbuf, &src->nwbuf, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _multilayerperceptron_clear(multilayerperceptron* p)
{
    ae_vector_clear(&p->structinfo);
    ae_vector_clear(&p->weights);
    ae_vector_clear(&p->columnmeans);
    ae_vector_clear(&p->columnsigmas);
    ae_vector_clear(&p->neurons);
    ae_vector_clear(&p->dfdnet);
    ae_vector_clear(&p->derror);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->y);
    ae_matrix_clear(&p->chunks);
    ae_vector_clear(&p->nwbuf);
}


/*$ End $*/
