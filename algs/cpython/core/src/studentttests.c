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


#include <stdafx.h>
#include "studentttests.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
One-sample t-test

This test checks three hypotheses about the mean of the given sample.  The
following tests are performed:
    * two-tailed test (null hypothesis - the mean is equal  to  the  given
      value)
    * left-tailed test (null hypothesis - the  mean  is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis - the mean is less than or  equal
      to the given value).

The test is based on the assumption that  a  given  sample  has  a  normal
distribution and  an  unknown  dispersion.  If  the  distribution  sharply
differs from normal, the test will work incorrectly.

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of sample.
    Mean    -   assumed value of the mean.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************/
void studentttest1(/* Real    */ ae_vector* x,
     ae_int_t n,
     double mean,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state)
{
    ae_int_t i;
    double xmean;
    double xvariance;
    double xstddev;
    double v1;
    double v2;
    double stat;
    double s;

    *bothtails = 0;
    *lefttail = 0;
    *righttail = 0;

    if( n<=1 )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Mean
     */
    xmean = 0;
    for(i=0; i<=n-1; i++)
    {
        xmean = xmean+x->ptr.p_double[i];
    }
    xmean = xmean/n;
    
    /*
     * Variance (using corrected two-pass algorithm)
     */
    xvariance = 0;
    xstddev = 0;
    if( n!=1 )
    {
        v1 = 0;
        for(i=0; i<=n-1; i++)
        {
            v1 = v1+ae_sqr(x->ptr.p_double[i]-xmean, _state);
        }
        v2 = 0;
        for(i=0; i<=n-1; i++)
        {
            v2 = v2+(x->ptr.p_double[i]-xmean);
        }
        v2 = ae_sqr(v2, _state)/n;
        xvariance = (v1-v2)/(n-1);
        if( ae_fp_less(xvariance,0) )
        {
            xvariance = 0;
        }
        xstddev = ae_sqrt(xvariance, _state);
    }
    if( ae_fp_eq(xstddev,0) )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Statistic
     */
    stat = (xmean-mean)/(xstddev/ae_sqrt(n, _state));
    s = studenttdistribution(n-1, stat, _state);
    *bothtails = 2*ae_minreal(s, 1-s, _state);
    *lefttail = s;
    *righttail = 1-s;
}


/*************************************************************************
Two-sample pooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * dispersions are equal
    * samples are independent.

Input parameters:
    X       -   sample 1. Array whose index goes from 0 to N-1.
    N       -   size of sample.
    Y       -   sample 2. Array whose index goes from 0 to M-1.
    M       -   size of sample.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************/
void studentttest2(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state)
{
    ae_int_t i;
    double xmean;
    double ymean;
    double stat;
    double s;
    double p;

    *bothtails = 0;
    *lefttail = 0;
    *righttail = 0;

    if( n<=1||m<=1 )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Mean
     */
    xmean = 0;
    for(i=0; i<=n-1; i++)
    {
        xmean = xmean+x->ptr.p_double[i];
    }
    xmean = xmean/n;
    ymean = 0;
    for(i=0; i<=m-1; i++)
    {
        ymean = ymean+y->ptr.p_double[i];
    }
    ymean = ymean/m;
    
    /*
     * S
     */
    s = 0;
    for(i=0; i<=n-1; i++)
    {
        s = s+ae_sqr(x->ptr.p_double[i]-xmean, _state);
    }
    for(i=0; i<=m-1; i++)
    {
        s = s+ae_sqr(y->ptr.p_double[i]-ymean, _state);
    }
    s = ae_sqrt(s*((double)1/(double)n+(double)1/(double)m)/(n+m-2), _state);
    if( ae_fp_eq(s,0) )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Statistic
     */
    stat = (xmean-ymean)/s;
    p = studenttdistribution(n+m-2, stat, _state);
    *bothtails = 2*ae_minreal(p, 1-p, _state);
    *lefttail = p;
    *righttail = 1-p;
}


/*************************************************************************
Two-sample unpooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * samples are independent.
Dispersion equality is not required

Input parameters:
    X - sample 1. Array whose index goes from 0 to N-1.
    N - size of the sample.
    Y - sample 2. Array whose index goes from 0 to M-1.
    M - size of the sample.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************/
void unequalvariancettest(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state)
{
    ae_int_t i;
    double xmean;
    double ymean;
    double xvar;
    double yvar;
    double df;
    double p;
    double stat;
    double c;

    *bothtails = 0;
    *lefttail = 0;
    *righttail = 0;

    if( n<=1||m<=1 )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Mean
     */
    xmean = 0;
    for(i=0; i<=n-1; i++)
    {
        xmean = xmean+x->ptr.p_double[i];
    }
    xmean = xmean/n;
    ymean = 0;
    for(i=0; i<=m-1; i++)
    {
        ymean = ymean+y->ptr.p_double[i];
    }
    ymean = ymean/m;
    
    /*
     * Variance (using corrected two-pass algorithm)
     */
    xvar = 0;
    for(i=0; i<=n-1; i++)
    {
        xvar = xvar+ae_sqr(x->ptr.p_double[i]-xmean, _state);
    }
    xvar = xvar/(n-1);
    yvar = 0;
    for(i=0; i<=m-1; i++)
    {
        yvar = yvar+ae_sqr(y->ptr.p_double[i]-ymean, _state);
    }
    yvar = yvar/(m-1);
    if( ae_fp_eq(xvar,0)||ae_fp_eq(yvar,0) )
    {
        *bothtails = 1.0;
        *lefttail = 1.0;
        *righttail = 1.0;
        return;
    }
    
    /*
     * Statistic
     */
    stat = (xmean-ymean)/ae_sqrt(xvar/n+yvar/m, _state);
    c = xvar/n/(xvar/n+yvar/m);
    df = (n-1)*(m-1)/((m-1)*ae_sqr(c, _state)+(n-1)*ae_sqr(1-c, _state));
    if( ae_fp_greater(stat,0) )
    {
        p = 1-0.5*incompletebeta(df/2, 0.5, df/(df+ae_sqr(stat, _state)), _state);
    }
    else
    {
        p = 0.5*incompletebeta(df/2, 0.5, df/(df+ae_sqr(stat, _state)), _state);
    }
    *bothtails = 2*ae_minreal(p, 1-p, _state);
    *lefttail = p;
    *righttail = 1-p;
}


/*$ End $*/
