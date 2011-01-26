/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _basestat_h
#define _basestat_h

#include "aenv.h"
#include "ialglib.h"
#include "tsort.h"
#include "apserv.h"
#include "basicstatops.h"
#include "ablasf.h"
#include "ablas.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Calculation of the distribution moments: mean, variance, skewness, kurtosis.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X
    
OUTPUT PARAMETERS
    Mean    -   mean.
    Variance-   variance.
    Skewness-   skewness (if variance<>0; zero otherwise).
    Kurtosis-   kurtosis (if variance<>0; zero otherwise).


  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void samplemoments(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* mean,
     double* variance,
     double* skewness,
     double* kurtosis,
     ae_state *_state);


/*************************************************************************
ADev

Input parameters:
    X   -   sample
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X
    
Output parameters:
    ADev-   ADev

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void sampleadev(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* adev,
     ae_state *_state);


/*************************************************************************
Median calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X

Output parameters:
    Median

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void samplemedian(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* median,
     ae_state *_state);


/*************************************************************************
Percentile calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X
    P   -   percentile (0<=P<=1)

Output parameters:
    V   -   percentile

  -- ALGLIB --
     Copyright 01.03.2008 by Bochkanov Sergey
*************************************************************************/
void samplepercentile(/* Real    */ ae_vector* x,
     ae_int_t n,
     double p,
     double* v,
     ae_state *_state);


/*************************************************************************
2-sample covariance

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    covariance (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
double cov2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
Pearson product-moment correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    Pearson product-moment correlation coefficient
    (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
double pearsoncorr2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
Spearman's rank correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    Spearman's rank correlation coefficient
    (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double spearmancorr2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
Covariance matrix

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], covariance matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void covm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Pearson product-moment correlation matrix

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void pearsoncorrm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Spearman's rank correlation matrix

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void spearmancorrm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Cross-covariance matrix

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-covariance matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void covm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Pearson product-moment cross-correlation matrix

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void pearsoncorrm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Spearman's rank cross-correlation matrix

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void spearmancorrm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);


/*************************************************************************
Obsolete function, we recommend to use PearsonCorr2().

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double pearsoncorrelation(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);


/*************************************************************************
Obsolete function, we recommend to use SpearmanCorr2().

    -- ALGLIB --
    Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double spearmanrankcorrelation(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);


/*$ End $*/
#endif

