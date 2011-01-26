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

#ifndef _mlpe_h
#define _mlpe_h

#include "aenv.h"
#include "ialglib.h"
#include "mlpbase.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "apserv.h"
#include "matinv.h"
#include "linmin.h"
#include "fbls.h"
#include "minlbfgs.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "mlptrain.h"
#include "tsort.h"
#include "basicstatops.h"
#include "basestat.h"
#include "bdss.h"


/*$ Declarations $*/


/*************************************************************************
Neural networks ensemble
*************************************************************************/
typedef struct
{
    ae_vector structinfo;
    ae_int_t ensemblesize;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_bool issoftmax;
    ae_bool postprocessing;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    ae_int_t serializedlen;
    ae_vector serializedmlp;
    ae_vector tmpweights;
    ae_vector tmpmeans;
    ae_vector tmpsigmas;
    ae_vector neurons;
    ae_vector dfdnet;
    ae_vector y;
} mlpensemble;


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(multilayerperceptron* network,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(mlpensemble* ensemble, ae_state *_state);


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(mlpensemble* ensemble,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_state *_state);


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);


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
     ae_state *_state);
ae_bool _mlpensemble_init(mlpensemble* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mlpensemble_init_copy(mlpensemble* dst, mlpensemble* src, ae_state *_state, ae_bool make_automatic);
void _mlpensemble_clear(mlpensemble* p);


/*$ End $*/
#endif

