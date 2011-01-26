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
#include "mlptrain.h"


/*$ Declarations $*/
static double mlptrain_mindecay = 0.001;
static void mlptrain_mlpkfoldcvgeneral(multilayerperceptron* n,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_bool lmalgorithm,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state);
static void mlptrain_mlpkfoldsplit(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nclasses,
     ae_int_t foldscount,
     ae_bool stratifiedsplits,
     /* Integer */ ae_vector* folds,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptrainlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double lmftol;
    double lmsteptol;
    ae_int_t i;
    ae_int_t k;
    double v;
    double e;
    double enew;
    double xnorm2;
    double stepnorm;
    ae_vector g;
    ae_vector d;
    ae_matrix h;
    ae_matrix hmod;
    ae_matrix z;
    ae_bool spd;
    double nu;
    double lambdav;
    double lambdaup;
    double lambdadown;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    ae_vector x;
    ae_vector y;
    ae_vector wbase;
    ae_vector wdir;
    ae_vector wt;
    ae_vector wx;
    ae_int_t pass;
    ae_vector wbest;
    double ebest;
    ae_int_t invinfo;
    matinvreport invrep;
    ae_int_t solverinfo;
    densesolverreport solverrep;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&g, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&h, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&hmod, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z, 0, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbase, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wdir, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wt, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    _matinvreport_init(&invrep, _state, ae_true);
    _densesolverreport_init(&solverrep, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    lambdaup = 10;
    lambdadown = 0.3;
    lmftol = 0.001;
    lmsteptol = 0.001;
    
    /*
     * Test for inputs
     */
    if( npoints<=0||restarts<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    decay = ae_maxreal(decay, mlptrain_mindecay, _state);
    *info = 2;
    
    /*
     * Initialize data
     */
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    
    /*
     * General case.
     * Prepare task and network. Allocate space.
     */
    mlpinitpreprocessor(network, xy, npoints, _state);
    ae_vector_set_length(&g, wcount-1+1, _state);
    ae_matrix_set_length(&h, wcount-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&hmod, wcount-1+1, wcount-1+1, _state);
    ae_vector_set_length(&wbase, wcount-1+1, _state);
    ae_vector_set_length(&wdir, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ae_vector_set_length(&wt, wcount-1+1, _state);
    ae_vector_set_length(&wx, wcount-1+1, _state);
    ebest = ae_maxrealnumber;
    
    /*
     * Multiple passes
     */
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Initialize weights
         */
        mlprandomize(network, _state);
        
        /*
         * First stage of the hybrid algorithm: LBFGS
         */
        ae_v_move(&wbase.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 5, _state), &wbase, &state, _state);
        minlbfgssetcond(&state, 0, 0, 0, ae_maxint(25, wcount, _state), _state);
        while(minlbfgsiteration(&state, _state))
        {
            
            /*
             * gradient
             */
            ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            mlpgradbatch(network, xy, npoints, &state.f, &state.g, _state);
            
            /*
             * weight decay
             */
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            
            /*
             * next iteration
             */
            rep->ngrad = rep->ngrad+1;
        }
        minlbfgsresults(&state, &wbase, &internalrep, _state);
        ae_v_move(&network->weights.ptr.p_double[0], 1, &wbase.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Second stage of the hybrid algorithm: LM
         *
         * Initialize H with identity matrix,
         * G with gradient,
         * E with regularized error.
         */
        mlphessianbatch(network, xy, npoints, &e, &g, &h, _state);
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = e+0.5*decay*v;
        ae_v_addd(&g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
        for(k=0; k<=wcount-1; k++)
        {
            h.ptr.pp_double[k][k] = h.ptr.pp_double[k][k]+decay;
        }
        rep->nhess = rep->nhess+1;
        lambdav = 0.001;
        nu = 2;
        for(;;)
        {
            
            /*
             * 1. HMod = H+lambda*I
             * 2. Try to solve (H+Lambda*I)*dx = -g.
             *    Increase lambda if left part is not positive definite.
             */
            for(i=0; i<=wcount-1; i++)
            {
                ae_v_move(&hmod.ptr.pp_double[i][0], 1, &h.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                hmod.ptr.pp_double[i][i] = hmod.ptr.pp_double[i][i]+lambdav;
            }
            spd = spdmatrixcholesky(&hmod, wcount, ae_true, _state);
            rep->ncholesky = rep->ncholesky+1;
            if( !spd )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            spdmatrixcholeskysolve(&hmod, wcount, ae_true, &g, &solverinfo, &solverrep, &wdir, _state);
            if( solverinfo<0 )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            ae_v_muld(&wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1), -1);
            
            /*
             * Lambda found.
             * 1. Save old w in WBase
             * 1. Test some stopping criterions
             * 2. If error(w+wdir)>error(w), increase lambda
             */
            ae_v_add(&network->weights.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            xnorm2 = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            stepnorm = ae_v_dotproduct(&wdir.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            stepnorm = ae_sqrt(stepnorm, _state);
            enew = mlperror(network, xy, npoints, _state)+0.5*decay*xnorm2;
            if( ae_fp_less(stepnorm,lmsteptol*(1+ae_sqrt(xnorm2, _state))) )
            {
                break;
            }
            if( ae_fp_greater(enew,e) )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            
            /*
             * Optimize using inv(cholesky(H)) as preconditioner
             */
            rmatrixtrinverse(&hmod, wcount, ae_true, ae_false, &invinfo, &invrep, _state);
            if( invinfo<=0 )
            {
                
                /*
                 * if matrix can't be inverted then exit with errors
                 * TODO: make WCount steps in direction suggested by HMod
                 */
                *info = -9;
                ae_frame_leave(_state);
                return;
            }
            ae_v_move(&wbase.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            for(i=0; i<=wcount-1; i++)
            {
                wt.ptr.p_double[i] = 0;
            }
            minlbfgscreatex(wcount, wcount, &wt, 1, &state, _state);
            minlbfgssetcond(&state, 0, 0, 0, 5, _state);
            while(minlbfgsiteration(&state, _state))
            {
                
                /*
                 * gradient
                 */
                for(i=0; i<=wcount-1; i++)
                {
                    v = ae_v_dotproduct(&state.x.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1));
                    network->weights.ptr.p_double[i] = wbase.ptr.p_double[i]+v;
                }
                mlpgradbatch(network, xy, npoints, &state.f, &g, _state);
                for(i=0; i<=wcount-1; i++)
                {
                    state.g.ptr.p_double[i] = 0;
                }
                for(i=0; i<=wcount-1; i++)
                {
                    v = g.ptr.p_double[i];
                    ae_v_addd(&state.g.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1), v);
                }
                
                /*
                 * weight decay
                 * grad(x'*x) = A'*(x0+A*t)
                 */
                v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                state.f = state.f+0.5*decay*v;
                for(i=0; i<=wcount-1; i++)
                {
                    v = decay*network->weights.ptr.p_double[i];
                    ae_v_addd(&state.g.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1), v);
                }
                
                /*
                 * next iteration
                 */
                rep->ngrad = rep->ngrad+1;
            }
            minlbfgsresults(&state, &wt, &internalrep, _state);
            
            /*
             * Accept new position.
             * Calculate Hessian
             */
            for(i=0; i<=wcount-1; i++)
            {
                v = ae_v_dotproduct(&wt.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1));
                network->weights.ptr.p_double[i] = wbase.ptr.p_double[i]+v;
            }
            mlphessianbatch(network, xy, npoints, &e, &g, &h, _state);
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            e = e+0.5*decay*v;
            ae_v_addd(&g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            for(k=0; k<=wcount-1; k++)
            {
                h.ptr.pp_double[k][k] = h.ptr.pp_double[k][k]+decay;
            }
            rep->nhess = rep->nhess+1;
            
            /*
             * Update lambda
             */
            lambdav = lambdav*lambdadown;
            nu = 2;
        }
        
        /*
         * update WBest
         */
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = 0.5*decay*v+mlperror(network, xy, npoints, _state);
        if( ae_fp_less(e,ebest) )
        {
            ebest = e;
            ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
    }
    
    /*
     * copy WBest to output
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlptrainlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector w;
    ae_vector wbest;
    double e;
    double v;
    double ebest;
    minlbfgsreport internalrep;
    minlbfgsstate state;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);

    
    /*
     * Test inputs, parse flags, read network geometry
     */
    if( ae_fp_eq(wstep,0)&&maxits==0 )
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
    mlpproperties(network, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    decay = ae_maxreal(decay, mlptrain_mindecay, _state);
    *info = 2;
    
    /*
     * Prepare
     */
    mlpinitpreprocessor(network, xy, npoints, _state);
    ae_vector_set_length(&w, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ebest = ae_maxrealnumber;
    
    /*
     * Multiple starts
     */
    rep->ncholesky = 0;
    rep->nhess = 0;
    rep->ngrad = 0;
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Process
         */
        mlprandomize(network, _state);
        ae_v_move(&w.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 10, _state), &w, &state, _state);
        minlbfgssetcond(&state, 0.0, 0.0, wstep, maxits, _state);
        while(minlbfgsiteration(&state, _state))
        {
            ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            mlpgradnbatch(network, xy, npoints, &state.f, &state.g, _state);
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            rep->ngrad = rep->ngrad+1;
        }
        minlbfgsresults(&state, &w, &internalrep, _state);
        ae_v_move(&network->weights.ptr.p_double[0], 1, &w.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Compare with best
         */
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = mlperrorn(network, xy, npoints, _state)+0.5*decay*v;
        if( ae_fp_less(e,ebest) )
        {
            ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            ebest = e;
        }
    }
    
    /*
     * The best network
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size
    ValXY       -   validation set
    ValSize     -   validation set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptraines(multilayerperceptron* network,
     /* Real    */ ae_matrix* trnxy,
     ae_int_t trnsize,
     /* Real    */ ae_matrix* valxy,
     ae_int_t valsize,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector w;
    ae_vector wbest;
    double e;
    double v;
    double ebest;
    ae_vector wfinal;
    double efinal;
    ae_int_t itbest;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    double wstep;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wfinal, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);

    wstep = 0.001;
    
    /*
     * Test inputs, parse flags, read network geometry
     */
    if( ((trnsize<=0||valsize<=0)||restarts<1)||ae_fp_less(decay,0) )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    mlpproperties(network, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=trnsize-1; i++)
        {
            if( ae_round(trnxy->ptr.pp_double[i][nin], _state)<0||ae_round(trnxy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
        for(i=0; i<=valsize-1; i++)
        {
            if( ae_round(valxy->ptr.pp_double[i][nin], _state)<0||ae_round(valxy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    *info = 2;
    
    /*
     * Prepare
     */
    mlpinitpreprocessor(network, trnxy, trnsize, _state);
    ae_vector_set_length(&w, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ae_vector_set_length(&wfinal, wcount-1+1, _state);
    efinal = ae_maxrealnumber;
    for(i=0; i<=wcount-1; i++)
    {
        wfinal.ptr.p_double[i] = 0;
    }
    
    /*
     * Multiple starts
     */
    rep->ncholesky = 0;
    rep->nhess = 0;
    rep->ngrad = 0;
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Process
         */
        mlprandomize(network, _state);
        ebest = mlperror(network, valxy, valsize, _state);
        ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        itbest = 0;
        ae_v_move(&w.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 10, _state), &w, &state, _state);
        minlbfgssetcond(&state, 0.0, 0.0, wstep, 0, _state);
        minlbfgssetxrep(&state, ae_true, _state);
        while(minlbfgsiteration(&state, _state))
        {
            
            /*
             * Calculate gradient
             */
            ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            mlpgradnbatch(network, trnxy, trnsize, &state.f, &state.g, _state);
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            rep->ngrad = rep->ngrad+1;
            
            /*
             * Validation set
             */
            if( state.xupdated )
            {
                ae_v_move(&network->weights.ptr.p_double[0], 1, &w.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                e = mlperror(network, valxy, valsize, _state);
                if( ae_fp_less(e,ebest) )
                {
                    ebest = e;
                    ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                    itbest = internalrep.iterationscount;
                }
                if( internalrep.iterationscount>30&&ae_fp_greater(internalrep.iterationscount,1.5*itbest) )
                {
                    *info = 6;
                    break;
                }
            }
        }
        minlbfgsresults(&state, &w, &internalrep, _state);
        
        /*
         * Compare with final answer
         */
        if( ae_fp_less(ebest,efinal) )
        {
            ae_v_move(&wfinal.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            efinal = ebest;
        }
    }
    
    /*
     * The best network
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wfinal.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);

    mlptrain_mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, ae_false, wstep, maxits, info, rep, cvrep, _state);
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);

    mlptrain_mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, ae_true, 0.0, 0, info, rep, cvrep, _state);
}


/*************************************************************************
Internal cross-validation subroutine
*************************************************************************/
static void mlptrain_mlpkfoldcvgeneral(multilayerperceptron* n,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_bool lmalgorithm,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t fold;
    ae_int_t j;
    ae_int_t k;
    multilayerperceptron network;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t rowlen;
    ae_int_t wcount;
    ae_int_t nclasses;
    ae_int_t tssize;
    ae_int_t cvssize;
    ae_matrix cvset;
    ae_matrix testset;
    ae_vector folds;
    ae_int_t relcnt;
    mlpreport internalrep;
    ae_vector x;
    ae_vector y;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);
    _multilayerperceptron_init(&network, _state, ae_true);
    ae_matrix_init(&cvset, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&testset, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&folds, 0, DT_INT, _state, ae_true);
    _mlpreport_init(&internalrep, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    
    /*
     * Read network geometry, test parameters
     */
    mlpproperties(n, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(n, _state) )
    {
        nclasses = nout;
        rowlen = nin+1;
    }
    else
    {
        nclasses = -nout;
        rowlen = nin+nout;
    }
    if( (npoints<=0||foldscount<2)||foldscount>npoints )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    mlpcopy(n, &network, _state);
    
    /*
     * K-fold out cross-validation.
     * First, estimate generalization error
     */
    ae_matrix_set_length(&testset, npoints-1+1, rowlen-1+1, _state);
    ae_matrix_set_length(&cvset, npoints-1+1, rowlen-1+1, _state);
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&y, nout-1+1, _state);
    mlptrain_mlpkfoldsplit(xy, npoints, nclasses, foldscount, ae_false, &folds, _state);
    cvrep->relclserror = 0;
    cvrep->avgce = 0;
    cvrep->rmserror = 0;
    cvrep->avgerror = 0;
    cvrep->avgrelerror = 0;
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    relcnt = 0;
    for(fold=0; fold<=foldscount-1; fold++)
    {
        
        /*
         * Separate set
         */
        tssize = 0;
        cvssize = 0;
        for(i=0; i<=npoints-1; i++)
        {
            if( folds.ptr.p_int[i]==fold )
            {
                ae_v_move(&testset.ptr.pp_double[tssize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,rowlen-1));
                tssize = tssize+1;
            }
            else
            {
                ae_v_move(&cvset.ptr.pp_double[cvssize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,rowlen-1));
                cvssize = cvssize+1;
            }
        }
        
        /*
         * Train on CV training set
         */
        if( lmalgorithm )
        {
            mlptrainlm(&network, &cvset, cvssize, decay, restarts, info, &internalrep, _state);
        }
        else
        {
            mlptrainlbfgs(&network, &cvset, cvssize, decay, restarts, wstep, maxits, info, &internalrep, _state);
        }
        if( *info<0 )
        {
            cvrep->relclserror = 0;
            cvrep->avgce = 0;
            cvrep->rmserror = 0;
            cvrep->avgerror = 0;
            cvrep->avgrelerror = 0;
            ae_frame_leave(_state);
            return;
        }
        rep->ngrad = rep->ngrad+internalrep.ngrad;
        rep->nhess = rep->nhess+internalrep.nhess;
        rep->ncholesky = rep->ncholesky+internalrep.ncholesky;
        
        /*
         * Estimate error using CV test set
         */
        if( mlpissoftmax(&network, _state) )
        {
            
            /*
             * classification-only code
             */
            cvrep->relclserror = cvrep->relclserror+mlpclserror(&network, &testset, tssize, _state);
            cvrep->avgce = cvrep->avgce+mlperrorn(&network, &testset, tssize, _state);
        }
        for(i=0; i<=tssize-1; i++)
        {
            ae_v_move(&x.ptr.p_double[0], 1, &testset.ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
            mlpprocess(&network, &x, &y, _state);
            if( mlpissoftmax(&network, _state) )
            {
                
                /*
                 * Classification-specific code
                 */
                k = ae_round(testset.ptr.pp_double[i][nin], _state);
                for(j=0; j<=nout-1; j++)
                {
                    if( j==k )
                    {
                        cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j]-1, _state);
                        cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j]-1, _state);
                        cvrep->avgrelerror = cvrep->avgrelerror+ae_fabs(y.ptr.p_double[j]-1, _state);
                        relcnt = relcnt+1;
                    }
                    else
                    {
                        cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j], _state);
                        cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j], _state);
                    }
                }
            }
            else
            {
                
                /*
                 * Regression-specific code
                 */
                for(j=0; j<=nout-1; j++)
                {
                    cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j], _state);
                    cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j], _state);
                    if( ae_fp_neq(testset.ptr.pp_double[i][nin+j],0) )
                    {
                        cvrep->avgrelerror = cvrep->avgrelerror+ae_fabs((y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j])/testset.ptr.pp_double[i][nin+j], _state);
                        relcnt = relcnt+1;
                    }
                }
            }
        }
    }
    if( mlpissoftmax(&network, _state) )
    {
        cvrep->relclserror = cvrep->relclserror/npoints;
        cvrep->avgce = cvrep->avgce/(ae_log(2, _state)*npoints);
    }
    cvrep->rmserror = ae_sqrt(cvrep->rmserror/(npoints*nout), _state);
    cvrep->avgerror = cvrep->avgerror/(npoints*nout);
    cvrep->avgrelerror = cvrep->avgrelerror/relcnt;
    *info = 1;
    ae_frame_leave(_state);
}


/*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.
*************************************************************************/
static void mlptrain_mlpkfoldsplit(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nclasses,
     ae_int_t foldscount,
     ae_bool stratifiedsplits,
     /* Integer */ ae_vector* folds,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;

    ae_vector_clear(folds);

    
    /*
     * test parameters
     */
    ae_assert(npoints>0, "MLPKFoldSplit: wrong NPoints!", _state);
    ae_assert(nclasses>1||nclasses<0, "MLPKFoldSplit: wrong NClasses!", _state);
    ae_assert(foldscount>=2&&foldscount<=npoints, "MLPKFoldSplit: wrong FoldsCount!", _state);
    ae_assert(!stratifiedsplits, "MLPKFoldSplit: stratified splits are not supported!", _state);
    
    /*
     * Folds
     */
    ae_vector_set_length(folds, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        folds->ptr.p_int[i] = i*foldscount/npoints;
    }
    for(i=0; i<=npoints-2; i++)
    {
        j = i+ae_randominteger(npoints-i, _state);
        if( j!=i )
        {
            k = folds->ptr.p_int[i];
            folds->ptr.p_int[i] = folds->ptr.p_int[j];
            folds->ptr.p_int[j] = k;
        }
    }
}


ae_bool _mlpreport_init(mlpreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _mlpreport_init_copy(mlpreport* dst, mlpreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
    return ae_true;
}


void _mlpreport_clear(mlpreport* p)
{
}


ae_bool _mlpcvreport_init(mlpcvreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _mlpcvreport_init_copy(mlpcvreport* dst, mlpcvreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    return ae_true;
}


void _mlpcvreport_clear(mlpcvreport* p)
{
}


/*$ End $*/
