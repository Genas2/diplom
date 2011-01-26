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
#include "mlpe.h"


/*$ Declarations $*/
static ae_int_t mlpe_mlpntotaloffset = 3;
static ae_int_t mlpe_mlpevnum = 9;
static void mlpe_mlpeallerrors(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state);
static void mlpe_mlpebagginginternal(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_bool lmalgorithm,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate0(nin, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate1(nin, nhid, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate2(nin, nhid1, nhid2, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb0(nin, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb1(nin, nhid, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb2(nin, nhid1, nhid2, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater0(nin, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater1(nin, nhid, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater2(nin, nhid1, nhid2, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec0(nin, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec1(nin, nhid, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec2(nin, nhid1, nhid2, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(multilayerperceptron* network,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ccount;

    _mlpensemble_clear(ensemble);

    ae_assert(ensemblesize>0, "MLPECreate: incorrect ensemble size!", _state);
    
    /*
     * network properties
     */
    mlpproperties(network, &ensemble->nin, &ensemble->nout, &ensemble->wcount, _state);
    if( mlpissoftmax(network, _state) )
    {
        ccount = ensemble->nin;
    }
    else
    {
        ccount = ensemble->nin+ensemble->nout;
    }
    ensemble->postprocessing = ae_false;
    ensemble->issoftmax = mlpissoftmax(network, _state);
    ensemble->ensemblesize = ensemblesize;
    
    /*
     * structure information
     */
    ae_vector_set_length(&ensemble->structinfo, network->structinfo.ptr.p_int[0]-1+1, _state);
    for(i=0; i<=network->structinfo.ptr.p_int[0]-1; i++)
    {
        ensemble->structinfo.ptr.p_int[i] = network->structinfo.ptr.p_int[i];
    }
    
    /*
     * weights, means, sigmas
     */
    ae_vector_set_length(&ensemble->weights, ensemblesize*ensemble->wcount-1+1, _state);
    ae_vector_set_length(&ensemble->columnmeans, ensemblesize*ccount-1+1, _state);
    ae_vector_set_length(&ensemble->columnsigmas, ensemblesize*ccount-1+1, _state);
    for(i=0; i<=ensemblesize*ensemble->wcount-1; i++)
    {
        ensemble->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=ensemblesize-1; i++)
    {
        ae_v_move(&ensemble->columnmeans.ptr.p_double[i*ccount], 1, &network->columnmeans.ptr.p_double[0], 1, ae_v_len(i*ccount,(i+1)*ccount-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[i*ccount], 1, &network->columnsigmas.ptr.p_double[0], 1, ae_v_len(i*ccount,(i+1)*ccount-1));
    }
    
    /*
     * serialized part
     */
    mlpserialize(network, &ensemble->serializedmlp, &ensemble->serializedlen, _state);
    
    /*
     * temporaries, internal buffers
     */
    ae_vector_set_length(&ensemble->tmpweights, ensemble->wcount-1+1, _state);
    ae_vector_set_length(&ensemble->tmpmeans, ccount-1+1, _state);
    ae_vector_set_length(&ensemble->tmpsigmas, ccount-1+1, _state);
    ae_vector_set_length(&ensemble->neurons, ensemble->structinfo.ptr.p_int[mlpe_mlpntotaloffset]-1+1, _state);
    ae_vector_set_length(&ensemble->dfdnet, ensemble->structinfo.ptr.p_int[mlpe_mlpntotaloffset]-1+1, _state);
    ae_vector_set_length(&ensemble->y, ensemble->nout-1+1, _state);
}


/*************************************************************************
Copying of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble1 -   original

OUTPUT PARAMETERS:
    Ensemble2 -   copy

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecopy(mlpensemble* ensemble1,
     mlpensemble* ensemble2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ccount;
    ae_int_t ntotal;

    _mlpensemble_clear(ensemble2);

    
    /*
     * Unload info
     */
    ssize = ensemble1->structinfo.ptr.p_int[0];
    if( ensemble1->issoftmax )
    {
        ccount = ensemble1->nin;
    }
    else
    {
        ccount = ensemble1->nin+ensemble1->nout;
    }
    ntotal = ensemble1->structinfo.ptr.p_int[mlpe_mlpntotaloffset];
    
    /*
     * Allocate space
     */
    ae_vector_set_length(&ensemble2->structinfo, ssize-1+1, _state);
    ae_vector_set_length(&ensemble2->weights, ensemble1->ensemblesize*ensemble1->wcount-1+1, _state);
    ae_vector_set_length(&ensemble2->columnmeans, ensemble1->ensemblesize*ccount-1+1, _state);
    ae_vector_set_length(&ensemble2->columnsigmas, ensemble1->ensemblesize*ccount-1+1, _state);
    ae_vector_set_length(&ensemble2->tmpweights, ensemble1->wcount-1+1, _state);
    ae_vector_set_length(&ensemble2->tmpmeans, ccount-1+1, _state);
    ae_vector_set_length(&ensemble2->tmpsigmas, ccount-1+1, _state);
    ae_vector_set_length(&ensemble2->serializedmlp, ensemble1->serializedlen-1+1, _state);
    ae_vector_set_length(&ensemble2->neurons, ntotal-1+1, _state);
    ae_vector_set_length(&ensemble2->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&ensemble2->y, ensemble1->nout-1+1, _state);
    
    /*
     * Copy
     */
    ensemble2->nin = ensemble1->nin;
    ensemble2->nout = ensemble1->nout;
    ensemble2->wcount = ensemble1->wcount;
    ensemble2->ensemblesize = ensemble1->ensemblesize;
    ensemble2->issoftmax = ensemble1->issoftmax;
    ensemble2->postprocessing = ensemble1->postprocessing;
    ensemble2->serializedlen = ensemble1->serializedlen;
    for(i=0; i<=ssize-1; i++)
    {
        ensemble2->structinfo.ptr.p_int[i] = ensemble1->structinfo.ptr.p_int[i];
    }
    ae_v_move(&ensemble2->weights.ptr.p_double[0], 1, &ensemble1->weights.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*ensemble1->wcount-1));
    ae_v_move(&ensemble2->columnmeans.ptr.p_double[0], 1, &ensemble1->columnmeans.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*ccount-1));
    ae_v_move(&ensemble2->columnsigmas.ptr.p_double[0], 1, &ensemble1->columnsigmas.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*ccount-1));
    ae_v_move(&ensemble2->serializedmlp.ptr.p_double[0], 1, &ensemble1->serializedmlp.ptr.p_double[0], 1, ae_v_len(0,ensemble1->serializedlen-1));
}


/*************************************************************************
Serialization of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble-   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores ensemble,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeserialize(mlpensemble* ensemble,
     /* Real    */ ae_vector* ra,
     ae_int_t* rlen,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t ccount;
    ae_int_t hsize;
    ae_int_t offs;

    ae_vector_clear(ra);
    *rlen = 0;

    hsize = 13;
    ssize = ensemble->structinfo.ptr.p_int[0];
    if( ensemble->issoftmax )
    {
        ccount = ensemble->nin;
    }
    else
    {
        ccount = ensemble->nin+ensemble->nout;
    }
    ntotal = ensemble->structinfo.ptr.p_int[mlpe_mlpntotaloffset];
    *rlen = hsize+ssize+ensemble->ensemblesize*ensemble->wcount+2*ccount*ensemble->ensemblesize+ensemble->serializedlen;
    
    /*
     *  RA format:
     *  [0]     RLen
     *  [1]     Version (MLPEVNum)
     *  [2]     EnsembleSize
     *  [3]     NIn
     *  [4]     NOut
     *  [5]     WCount
     *  [6]     IsSoftmax 0/1
     *  [7]     PostProcessing 0/1
     *  [8]     sizeof(StructInfo)
     *  [9]     NTotal (sizeof(Neurons), sizeof(DFDNET))
     *  [10]    CCount (sizeof(ColumnMeans), sizeof(ColumnSigmas))
     *  [11]    data offset
     *  [12]    SerializedLen
     *
     *  [..]    StructInfo
     *  [..]    Weights
     *  [..]    ColumnMeans
     *  [..]    ColumnSigmas
     */
    ae_vector_set_length(ra, *rlen-1+1, _state);
    ra->ptr.p_double[0] = *rlen;
    ra->ptr.p_double[1] = mlpe_mlpevnum;
    ra->ptr.p_double[2] = ensemble->ensemblesize;
    ra->ptr.p_double[3] = ensemble->nin;
    ra->ptr.p_double[4] = ensemble->nout;
    ra->ptr.p_double[5] = ensemble->wcount;
    if( ensemble->issoftmax )
    {
        ra->ptr.p_double[6] = 1;
    }
    else
    {
        ra->ptr.p_double[6] = 0;
    }
    if( ensemble->postprocessing )
    {
        ra->ptr.p_double[7] = 1;
    }
    else
    {
        ra->ptr.p_double[7] = 9;
    }
    ra->ptr.p_double[8] = ssize;
    ra->ptr.p_double[9] = ntotal;
    ra->ptr.p_double[10] = ccount;
    ra->ptr.p_double[11] = hsize;
    ra->ptr.p_double[12] = ensemble->serializedlen;
    offs = hsize;
    for(i=offs; i<=offs+ssize-1; i++)
    {
        ra->ptr.p_double[i] = ensemble->structinfo.ptr.p_int[i-offs];
    }
    offs = offs+ssize;
    ae_v_move(&ra->ptr.p_double[offs], 1, &ensemble->weights.ptr.p_double[0], 1, ae_v_len(offs,offs+ensemble->ensemblesize*ensemble->wcount-1));
    offs = offs+ensemble->ensemblesize*ensemble->wcount;
    ae_v_move(&ra->ptr.p_double[offs], 1, &ensemble->columnmeans.ptr.p_double[0], 1, ae_v_len(offs,offs+ensemble->ensemblesize*ccount-1));
    offs = offs+ensemble->ensemblesize*ccount;
    ae_v_move(&ra->ptr.p_double[offs], 1, &ensemble->columnsigmas.ptr.p_double[0], 1, ae_v_len(offs,offs+ensemble->ensemblesize*ccount-1));
    offs = offs+ensemble->ensemblesize*ccount;
    ae_v_move(&ra->ptr.p_double[offs], 1, &ensemble->serializedmlp.ptr.p_double[0], 1, ae_v_len(offs,offs+ensemble->serializedlen-1));
    offs = offs+ensemble->serializedlen;
}


/*************************************************************************
Unserialization of MLPEnsemble strucure

INPUT PARAMETERS:
    RA      -   real array which stores ensemble

OUTPUT PARAMETERS:
    Ensemble-   restored structure

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeunserialize(/* Real    */ ae_vector* ra,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t ccount;
    ae_int_t hsize;
    ae_int_t offs;

    _mlpensemble_clear(ensemble);

    ae_assert(ae_round(ra->ptr.p_double[1], _state)==mlpe_mlpevnum, "MLPEUnserialize: incorrect array!", _state);
    
    /*
     * load info
     */
    hsize = 13;
    ensemble->ensemblesize = ae_round(ra->ptr.p_double[2], _state);
    ensemble->nin = ae_round(ra->ptr.p_double[3], _state);
    ensemble->nout = ae_round(ra->ptr.p_double[4], _state);
    ensemble->wcount = ae_round(ra->ptr.p_double[5], _state);
    ensemble->issoftmax = ae_round(ra->ptr.p_double[6], _state)==1;
    ensemble->postprocessing = ae_round(ra->ptr.p_double[7], _state)==1;
    ssize = ae_round(ra->ptr.p_double[8], _state);
    ntotal = ae_round(ra->ptr.p_double[9], _state);
    ccount = ae_round(ra->ptr.p_double[10], _state);
    offs = ae_round(ra->ptr.p_double[11], _state);
    ensemble->serializedlen = ae_round(ra->ptr.p_double[12], _state);
    
    /*
     *  Allocate arrays
     */
    ae_vector_set_length(&ensemble->structinfo, ssize-1+1, _state);
    ae_vector_set_length(&ensemble->weights, ensemble->ensemblesize*ensemble->wcount-1+1, _state);
    ae_vector_set_length(&ensemble->columnmeans, ensemble->ensemblesize*ccount-1+1, _state);
    ae_vector_set_length(&ensemble->columnsigmas, ensemble->ensemblesize*ccount-1+1, _state);
    ae_vector_set_length(&ensemble->tmpweights, ensemble->wcount-1+1, _state);
    ae_vector_set_length(&ensemble->tmpmeans, ccount-1+1, _state);
    ae_vector_set_length(&ensemble->tmpsigmas, ccount-1+1, _state);
    ae_vector_set_length(&ensemble->neurons, ntotal-1+1, _state);
    ae_vector_set_length(&ensemble->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&ensemble->serializedmlp, ensemble->serializedlen-1+1, _state);
    ae_vector_set_length(&ensemble->y, ensemble->nout-1+1, _state);
    
    /*
     * load data
     */
    for(i=offs; i<=offs+ssize-1; i++)
    {
        ensemble->structinfo.ptr.p_int[i-offs] = ae_round(ra->ptr.p_double[i], _state);
    }
    offs = offs+ssize;
    ae_v_move(&ensemble->weights.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,ensemble->ensemblesize*ensemble->wcount-1));
    offs = offs+ensemble->ensemblesize*ensemble->wcount;
    ae_v_move(&ensemble->columnmeans.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,ensemble->ensemblesize*ccount-1));
    offs = offs+ensemble->ensemblesize*ccount;
    ae_v_move(&ensemble->columnsigmas.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,ensemble->ensemblesize*ccount-1));
    offs = offs+ensemble->ensemblesize*ccount;
    ae_v_move(&ensemble->serializedmlp.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,ensemble->serializedlen-1));
    offs = offs+ensemble->serializedlen;
}


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(mlpensemble* ensemble, ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=ensemble->ensemblesize*ensemble->wcount-1; i++)
    {
        ensemble->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
}


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(mlpensemble* ensemble,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_state *_state)
{

    *nin = 0;
    *nout = 0;

    *nin = ensemble->nin;
    *nout = ensemble->nout;
}


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state)
{
    ae_bool result;


    result = ensemble->issoftmax;
    return result;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NOut, it will be reallocated. If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.


OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocess(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t es;
    ae_int_t wc;
    ae_int_t cc;
    double v;


    if( y->cnt<ensemble->nout )
    {
        ae_vector_set_length(y, ensemble->nout, _state);
    }
    es = ensemble->ensemblesize;
    wc = ensemble->wcount;
    if( ensemble->issoftmax )
    {
        cc = ensemble->nin;
    }
    else
    {
        cc = ensemble->nin+ensemble->nout;
    }
    v = (double)1/(double)es;
    for(i=0; i<=ensemble->nout-1; i++)
    {
        y->ptr.p_double[i] = 0;
    }
    for(i=0; i<=es-1; i++)
    {
        ae_v_move(&ensemble->tmpweights.ptr.p_double[0], 1, &ensemble->weights.ptr.p_double[i*wc], 1, ae_v_len(0,wc-1));
        ae_v_move(&ensemble->tmpmeans.ptr.p_double[0], 1, &ensemble->columnmeans.ptr.p_double[i*cc], 1, ae_v_len(0,cc-1));
        ae_v_move(&ensemble->tmpsigmas.ptr.p_double[0], 1, &ensemble->columnsigmas.ptr.p_double[i*cc], 1, ae_v_len(0,cc-1));
        mlpinternalprocessvector(&ensemble->structinfo, &ensemble->tmpweights, &ensemble->tmpmeans, &ensemble->tmpsigmas, &ensemble->neurons, &ensemble->dfdnet, x, &ensemble->y, _state);
        ae_v_addd(&y->ptr.p_double[0], 1, &ensemble->y.ptr.p_double[0], 1, ae_v_len(0,ensemble->nout-1), v);
    }
}


/*************************************************************************
'interactive'  variant  of  MLPEProcess  for  languages  like Python which
support constructs like "Y = MLPEProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocessi(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    mlpeprocess(ensemble, x, y, _state);
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlperelclserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = relcls;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgce(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avgce;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpermserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = rms;
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avg;
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgrelerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avgrel;
    return result;
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglm(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);

    mlpe_mlpebagginginternal(ensemble, xy, npoints, decay, restarts, 0.0, 0, ae_true, info, rep, ooberrors, _state);
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglbfgs(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);

    mlpe_mlpebagginginternal(ensemble, xy, npoints, decay, restarts, wstep, maxits, ae_false, info, rep, ooberrors, _state);
}


/*************************************************************************
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlpetraines(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t k;
    ae_int_t ccount;
    ae_int_t pcount;
    ae_matrix trnxy;
    ae_matrix valxy;
    ae_int_t trnsize;
    ae_int_t valsize;
    multilayerperceptron network;
    ae_int_t tmpinfo;
    mlpreport tmprep;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_matrix_init(&trnxy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&valxy, 0, 0, DT_REAL, _state, ae_true);
    _multilayerperceptron_init(&network, _state, ae_true);
    _mlpreport_init(&tmprep, _state, ae_true);

    if( (npoints<2||restarts<1)||ae_fp_less(decay,0) )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( ensemble->issoftmax )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][ensemble->nin], _state)<0||ae_round(xy->ptr.pp_double[i][ensemble->nin], _state)>=ensemble->nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    *info = 6;
    
    /*
     * allocate
     */
    if( ensemble->issoftmax )
    {
        ccount = ensemble->nin+1;
        pcount = ensemble->nin;
    }
    else
    {
        ccount = ensemble->nin+ensemble->nout;
        pcount = ensemble->nin+ensemble->nout;
    }
    ae_matrix_set_length(&trnxy, npoints-1+1, ccount-1+1, _state);
    ae_matrix_set_length(&valxy, npoints-1+1, ccount-1+1, _state);
    mlpunserialize(&ensemble->serializedmlp, &network, _state);
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    
    /*
     * train networks
     */
    for(k=0; k<=ensemble->ensemblesize-1; k++)
    {
        
        /*
         * Split set
         */
        do
        {
            trnsize = 0;
            valsize = 0;
            for(i=0; i<=npoints-1; i++)
            {
                if( ae_fp_less(ae_randomreal(_state),0.66) )
                {
                    
                    /*
                     * Assign sample to training set
                     */
                    ae_v_move(&trnxy.ptr.pp_double[trnsize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,ccount-1));
                    trnsize = trnsize+1;
                }
                else
                {
                    
                    /*
                     * Assign sample to validation set
                     */
                    ae_v_move(&valxy.ptr.pp_double[valsize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,ccount-1));
                    valsize = valsize+1;
                }
            }
        }
        while(!(trnsize!=0&&valsize!=0));
        
        /*
         * Train
         */
        mlptraines(&network, &trnxy, trnsize, &valxy, valsize, decay, restarts, &tmpinfo, &tmprep, _state);
        if( tmpinfo<0 )
        {
            *info = tmpinfo;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * save results
         */
        ae_v_move(&ensemble->weights.ptr.p_double[k*ensemble->wcount], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(k*ensemble->wcount,(k+1)*ensemble->wcount-1));
        ae_v_move(&ensemble->columnmeans.ptr.p_double[k*pcount], 1, &network.columnmeans.ptr.p_double[0], 1, ae_v_len(k*pcount,(k+1)*pcount-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[k*pcount], 1, &network.columnsigmas.ptr.p_double[0], 1, ae_v_len(k*pcount,(k+1)*pcount-1));
        rep->ngrad = rep->ngrad+tmprep.ngrad;
        rep->nhess = rep->nhess+tmprep.nhess;
        rep->ncholesky = rep->ncholesky+tmprep.ncholesky;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpe_mlpeallerrors(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector buf;
    ae_vector workx;
    ae_vector y;
    ae_vector dy;

    ae_frame_make(_state, &_frame_block);
    *relcls = 0;
    *avgce = 0;
    *rms = 0;
    *avg = 0;
    *avgrel = 0;
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dy, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&workx, ensemble->nin-1+1, _state);
    ae_vector_set_length(&y, ensemble->nout-1+1, _state);
    if( ensemble->issoftmax )
    {
        ae_vector_set_length(&dy, 0+1, _state);
        dserrallocate(ensemble->nout, &buf, _state);
    }
    else
    {
        ae_vector_set_length(&dy, ensemble->nout-1+1, _state);
        dserrallocate(-ensemble->nout, &buf, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,ensemble->nin-1));
        mlpeprocess(ensemble, &workx, &y, _state);
        if( ensemble->issoftmax )
        {
            dy.ptr.p_double[0] = xy->ptr.pp_double[i][ensemble->nin];
        }
        else
        {
            ae_v_move(&dy.ptr.p_double[0], 1, &xy->ptr.pp_double[i][ensemble->nin], 1, ae_v_len(0,ensemble->nout-1));
        }
        dserraccumulate(&buf, &y, &dy, _state);
    }
    dserrfinish(&buf, _state);
    *relcls = buf.ptr.p_double[0];
    *avgce = buf.ptr.p_double[1];
    *rms = buf.ptr.p_double[2];
    *avg = buf.ptr.p_double[3];
    *avgrel = buf.ptr.p_double[4];
    ae_frame_leave(_state);
}


/*************************************************************************
Internal bagging subroutine.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpe_mlpebagginginternal(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_bool lmalgorithm,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xys;
    ae_vector s;
    ae_matrix oobbuf;
    ae_vector oobcntbuf;
    ae_vector x;
    ae_vector y;
    ae_vector dy;
    ae_vector dsbuf;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t ccnt;
    ae_int_t pcnt;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double v;
    mlpreport tmprep;
    multilayerperceptron network;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);
    ae_matrix_init(&xys, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&s, 0, DT_BOOL, _state, ae_true);
    ae_matrix_init(&oobbuf, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&oobcntbuf, 0, DT_INT, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dsbuf, 0, DT_REAL, _state, ae_true);
    _mlpreport_init(&tmprep, _state, ae_true);
    _multilayerperceptron_init(&network, _state, ae_true);

    
    /*
     * Test for inputs
     */
    if( (!lmalgorithm&&ae_fp_eq(wstep,0))&&maxits==0 )
    {
        *info = -8;
        ae_frame_leave(_state);
        return;
    }
    if( ((npoints<=0||restarts<1)||ae_fp_less(wstep,0))||maxits<0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( ensemble->issoftmax )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][ensemble->nin], _state)<0||ae_round(xy->ptr.pp_double[i][ensemble->nin], _state)>=ensemble->nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    
    /*
     * allocate temporaries
     */
    *info = 2;
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    ooberrors->relclserror = 0;
    ooberrors->avgce = 0;
    ooberrors->rmserror = 0;
    ooberrors->avgerror = 0;
    ooberrors->avgrelerror = 0;
    nin = ensemble->nin;
    nout = ensemble->nout;
    if( ensemble->issoftmax )
    {
        ccnt = nin+1;
        pcnt = nin;
    }
    else
    {
        ccnt = nin+nout;
        pcnt = nin+nout;
    }
    ae_matrix_set_length(&xys, npoints-1+1, ccnt-1+1, _state);
    ae_vector_set_length(&s, npoints-1+1, _state);
    ae_matrix_set_length(&oobbuf, npoints-1+1, nout-1+1, _state);
    ae_vector_set_length(&oobcntbuf, npoints-1+1, _state);
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&y, nout-1+1, _state);
    if( ensemble->issoftmax )
    {
        ae_vector_set_length(&dy, 0+1, _state);
    }
    else
    {
        ae_vector_set_length(&dy, nout-1+1, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        for(j=0; j<=nout-1; j++)
        {
            oobbuf.ptr.pp_double[i][j] = 0;
        }
    }
    for(i=0; i<=npoints-1; i++)
    {
        oobcntbuf.ptr.p_int[i] = 0;
    }
    mlpunserialize(&ensemble->serializedmlp, &network, _state);
    
    /*
     * main bagging cycle
     */
    for(k=0; k<=ensemble->ensemblesize-1; k++)
    {
        
        /*
         * prepare dataset
         */
        for(i=0; i<=npoints-1; i++)
        {
            s.ptr.p_bool[i] = ae_false;
        }
        for(i=0; i<=npoints-1; i++)
        {
            j = ae_randominteger(npoints, _state);
            s.ptr.p_bool[j] = ae_true;
            ae_v_move(&xys.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[j][0], 1, ae_v_len(0,ccnt-1));
        }
        
        /*
         * train
         */
        if( lmalgorithm )
        {
            mlptrainlm(&network, &xys, npoints, decay, restarts, info, &tmprep, _state);
        }
        else
        {
            mlptrainlbfgs(&network, &xys, npoints, decay, restarts, wstep, maxits, info, &tmprep, _state);
        }
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * save results
         */
        rep->ngrad = rep->ngrad+tmprep.ngrad;
        rep->nhess = rep->nhess+tmprep.nhess;
        rep->ncholesky = rep->ncholesky+tmprep.ncholesky;
        ae_v_move(&ensemble->weights.ptr.p_double[k*ensemble->wcount], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(k*ensemble->wcount,(k+1)*ensemble->wcount-1));
        ae_v_move(&ensemble->columnmeans.ptr.p_double[k*pcnt], 1, &network.columnmeans.ptr.p_double[0], 1, ae_v_len(k*pcnt,(k+1)*pcnt-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[k*pcnt], 1, &network.columnsigmas.ptr.p_double[0], 1, ae_v_len(k*pcnt,(k+1)*pcnt-1));
        
        /*
         * OOB estimates
         */
        for(i=0; i<=npoints-1; i++)
        {
            if( !s.ptr.p_bool[i] )
            {
                ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
                mlpprocess(&network, &x, &y, _state);
                ae_v_add(&oobbuf.ptr.pp_double[i][0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
                oobcntbuf.ptr.p_int[i] = oobcntbuf.ptr.p_int[i]+1;
            }
        }
    }
    
    /*
     * OOB estimates
     */
    if( ensemble->issoftmax )
    {
        dserrallocate(nout, &dsbuf, _state);
    }
    else
    {
        dserrallocate(-nout, &dsbuf, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( oobcntbuf.ptr.p_int[i]!=0 )
        {
            v = (double)1/(double)oobcntbuf.ptr.p_int[i];
            ae_v_moved(&y.ptr.p_double[0], 1, &oobbuf.ptr.pp_double[i][0], 1, ae_v_len(0,nout-1), v);
            if( ensemble->issoftmax )
            {
                dy.ptr.p_double[0] = xy->ptr.pp_double[i][nin];
            }
            else
            {
                ae_v_moved(&dy.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1), v);
            }
            dserraccumulate(&dsbuf, &y, &dy, _state);
        }
    }
    dserrfinish(&dsbuf, _state);
    ooberrors->relclserror = dsbuf.ptr.p_double[0];
    ooberrors->avgce = dsbuf.ptr.p_double[1];
    ooberrors->rmserror = dsbuf.ptr.p_double[2];
    ooberrors->avgerror = dsbuf.ptr.p_double[3];
    ooberrors->avgrelerror = dsbuf.ptr.p_double[4];
    ae_frame_leave(_state);
}


ae_bool _mlpensemble_init(mlpensemble* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->structinfo, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->weights, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnmeans, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnsigmas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->serializedmlp, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpweights, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpmeans, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpsigmas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->neurons, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dfdnet, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _mlpensemble_init_copy(mlpensemble* dst, mlpensemble* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->structinfo, &src->structinfo, _state, make_automatic) )
        return ae_false;
    dst->ensemblesize = src->ensemblesize;
    dst->nin = src->nin;
    dst->nout = src->nout;
    dst->wcount = src->wcount;
    dst->issoftmax = src->issoftmax;
    dst->postprocessing = src->postprocessing;
    if( !ae_vector_init_copy(&dst->weights, &src->weights, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnmeans, &src->columnmeans, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnsigmas, &src->columnsigmas, _state, make_automatic) )
        return ae_false;
    dst->serializedlen = src->serializedlen;
    if( !ae_vector_init_copy(&dst->serializedmlp, &src->serializedmlp, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpweights, &src->tmpweights, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpmeans, &src->tmpmeans, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpsigmas, &src->tmpsigmas, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->neurons, &src->neurons, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dfdnet, &src->dfdnet, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _mlpensemble_clear(mlpensemble* p)
{
    ae_vector_clear(&p->structinfo);
    ae_vector_clear(&p->weights);
    ae_vector_clear(&p->columnmeans);
    ae_vector_clear(&p->columnsigmas);
    ae_vector_clear(&p->serializedmlp);
    ae_vector_clear(&p->tmpweights);
    ae_vector_clear(&p->tmpmeans);
    ae_vector_clear(&p->tmpsigmas);
    ae_vector_clear(&p->neurons);
    ae_vector_clear(&p->dfdnet);
    ae_vector_clear(&p->y);
}


/*$ End $*/
