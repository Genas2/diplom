/*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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
#include "lsfit.h"


/*$ Declarations $*/
static ae_int_t lsfit_rfsmax = 10;
static void lsfit_spline1dfitinternal(ae_int_t st,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state);
static void lsfit_lsfitlinearinternal(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);
static void lsfit_lsfitclearrequestfields(lsfitstate* state,
     ae_state *_state);
static void lsfit_barycentriccalcbasis(barycentricinterpolant* b,
     double t,
     /* Real    */ ae_vector* y,
     ae_state *_state);
static void lsfit_internalchebyshevfit(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state);
static void lsfit_barycentricfitwcfixedd(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t d,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Fitting by polynomials in barycentric form. This function provides  simple
unterface for unconstrained unweighted fitting. See  PolynomialFitWC()  if
you need constrained fitting.

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO:
    PolynomialFitWC()

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0
            * if given, only leading N elements of X/Y are used
            * if not given, automatically determined from sizes of X/Y
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
    P   -   interpolant in barycentric form.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

NOTES:
    you can convert P from barycentric form  to  the  power  or  Chebyshev
    basis with PolynomialBar2Pow() or PolynomialBar2Cheb() functions  from
    POLINT subpackage.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************/
void polynomialfit(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* p,
     polynomialfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector w;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _barycentricinterpolant_clear(p);
    _polynomialfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);

    ae_assert(n>0, "PolynomialFit: N<=0!", _state);
    ae_assert(m>0, "PolynomialFit: M<=0!", _state);
    ae_assert(x->cnt>=n, "PolynomialFit: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "PolynomialFit: Length(Y)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "PolynomialFit: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "PolynomialFit: Y contains infinite or NaN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    polynomialfitwc(x, y, &w, n, &xc, &yc, &dc, 0, m, info, p, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Weighted  fitting by polynomials in barycentric form, with constraints  on
function values or first derivatives.

Small regularizing term is used when solving constrained tasks (to improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO:
    PolynomialFit()

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
            * if given, only leading N elements of X/Y/W are used
            * if not given, automatically determined from sizes of X/Y/W
    XC  -   points where polynomial values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that P(XC[i])=YC[i]
            * DC[i]=1   means that P'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    P   -   interpolant in barycentric form.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

NOTES:
    you can convert P from barycentric form  to  the  power  or  Chebyshev
    basis with PolynomialBar2Pow() or PolynomialBar2Cheb() functions  from
    POLINT subpackage.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* even simple constraints can be inconsistent, see  Wikipedia  article  on
  this subject: http://en.wikipedia.org/wiki/Birkhoff_interpolation
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the one special cases, however, we can  guarantee  consistency.  This
  case  is:  M>1  and constraints on the function values (NOT DERIVATIVES)

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************/
void polynomialfitwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* p,
     polynomialfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    double xa;
    double xb;
    double sa;
    double sb;
    ae_vector xoriginal;
    ae_vector yoriginal;
    ae_vector y2;
    ae_vector w2;
    ae_vector tmp;
    ae_vector tmp2;
    ae_vector bx;
    ae_vector by;
    ae_vector bw;
    ae_int_t i;
    ae_int_t j;
    double u;
    double v;
    double s;
    ae_int_t relcnt;
    lsfitreport lrep;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_vector_init_copy(&_w, w, _state, ae_true);
    w = &_w;
    ae_vector_init_copy(&_xc, xc, _state, ae_true);
    xc = &_xc;
    ae_vector_init_copy(&_yc, yc, _state, ae_true);
    yc = &_yc;
    *info = 0;
    _barycentricinterpolant_clear(p);
    _polynomialfitreport_clear(rep);
    ae_vector_init(&xoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&by, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bw, 0, DT_REAL, _state, ae_true);
    _lsfitreport_init(&lrep, _state, ae_true);

    ae_assert(n>0, "PolynomialFitWC: N<=0!", _state);
    ae_assert(m>0, "PolynomialFitWC: M<=0!", _state);
    ae_assert(k>=0, "PolynomialFitWC: K<0!", _state);
    ae_assert(k<m, "PolynomialFitWC: K>=M!", _state);
    ae_assert(x->cnt>=n, "PolynomialFitWC: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "PolynomialFitWC: Length(Y)<N!", _state);
    ae_assert(w->cnt>=n, "PolynomialFitWC: Length(W)<N!", _state);
    ae_assert(xc->cnt>=k, "PolynomialFitWC: Length(XC)<K!", _state);
    ae_assert(yc->cnt>=k, "PolynomialFitWC: Length(YC)<K!", _state);
    ae_assert(dc->cnt>=k, "PolynomialFitWC: Length(DC)<K!", _state);
    ae_assert(isfinitevector(x, n, _state), "PolynomialFitWC: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "PolynomialFitWC: Y contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(w, n, _state), "PolynomialFitWC: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(xc, k, _state), "PolynomialFitWC: XC contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(yc, k, _state), "PolynomialFitWC: YC contains infinite or NaN values!", _state);
    for(i=0; i<=k-1; i++)
    {
        ae_assert(dc->ptr.p_int[i]==0||dc->ptr.p_int[i]==1, "PolynomialFitWC: one of DC[] is not 0 or 1!", _state);
    }
    
    /*
     * Scale X, Y, XC, YC.
     * Solve scaled problem using internal Chebyshev fitting function.
     */
    lsfitscalexy(x, y, w, n, xc, yc, dc, k, &xa, &xb, &sa, &sb, &xoriginal, &yoriginal, _state);
    lsfit_internalchebyshevfit(x, y, w, n, xc, yc, dc, k, m, info, &tmp, &lrep, _state);
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Generate barycentric model and scale it
     * * BX, BY store barycentric model nodes
     * * FMatrix is reused (remember - it is at least MxM, what we need)
     *
     * Model intialization is done in O(M^2). In principle, it can be
     * done in O(M*log(M)), but before it we solved task with O(N*M^2)
     * complexity, so it is only a small amount of total time spent.
     */
    ae_vector_set_length(&bx, m, _state);
    ae_vector_set_length(&by, m, _state);
    ae_vector_set_length(&bw, m, _state);
    ae_vector_set_length(&tmp2, m, _state);
    s = 1;
    for(i=0; i<=m-1; i++)
    {
        if( m!=1 )
        {
            u = ae_cos(ae_pi*i/(m-1), _state);
        }
        else
        {
            u = 0;
        }
        v = 0;
        for(j=0; j<=m-1; j++)
        {
            if( j==0 )
            {
                tmp2.ptr.p_double[j] = 1;
            }
            else
            {
                if( j==1 )
                {
                    tmp2.ptr.p_double[j] = u;
                }
                else
                {
                    tmp2.ptr.p_double[j] = 2*u*tmp2.ptr.p_double[j-1]-tmp2.ptr.p_double[j-2];
                }
            }
            v = v+tmp.ptr.p_double[j]*tmp2.ptr.p_double[j];
        }
        bx.ptr.p_double[i] = u;
        by.ptr.p_double[i] = v;
        bw.ptr.p_double[i] = s;
        if( i==0||i==m-1 )
        {
            bw.ptr.p_double[i] = 0.5*bw.ptr.p_double[i];
        }
        s = -s;
    }
    barycentricbuildxyw(&bx, &by, &bw, m, p, _state);
    barycentriclintransx(p, 2/(xb-xa), -(xa+xb)/(xb-xa), _state);
    barycentriclintransy(p, sb-sa, sa, _state);
    
    /*
     * Scale absolute errors obtained from LSFitLinearW.
     * Relative error should be calculated separately
     * (because of shifting/scaling of the task)
     */
    rep->taskrcond = lrep.taskrcond;
    rep->rmserror = lrep.rmserror*(sb-sa);
    rep->avgerror = lrep.avgerror*(sb-sa);
    rep->maxerror = lrep.maxerror*(sb-sa);
    rep->avgrelerror = 0;
    relcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(yoriginal.ptr.p_double[i],0) )
        {
            rep->avgrelerror = rep->avgrelerror+ae_fabs(barycentriccalc(p, xoriginal.ptr.p_double[i], _state)-yoriginal.ptr.p_double[i], _state)/ae_fabs(yoriginal.ptr.p_double[i], _state);
            relcnt = relcnt+1;
        }
    }
    if( relcnt!=0 )
    {
        rep->avgrelerror = rep->avgrelerror/relcnt;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Weghted rational least  squares  fitting  using  Floater-Hormann  rational
functions  with  optimal  D  chosen  from  [0,9],  with  constraints   and
individual weights.

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal D (least WEIGHTED root
mean square error) is chosen.  Task  is  linear,  so  linear least squares
solver  is  used.  Complexity  of  this  computational  scheme is O(N*M^2)
(mostly dominated by the least squares solver).

SEE ALSO
* BarycentricFitFloaterHormann(), "lightweight" fitting without invididual
  weights and constraints.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where function values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroutine doesn't calculate task's condition number for K<>0.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained barycentric interpolants:
* excessive  constraints  can  be  inconsistent.   Floater-Hormann   basis
  functions aren't as flexible as splines (although they are very smooth).
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints on the function  VALUES at the interval
  boundaries. Note that consustency of the  constraints  on  the  function
  DERIVATIVES is NOT guaranteed (you can use in such cases  cubic  splines
  which are more flexible).
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormannwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t d;
    ae_int_t i;
    double wrmscur;
    double wrmsbest;
    barycentricinterpolant locb;
    barycentricfitreport locrep;
    ae_int_t locinfo;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _barycentricinterpolant_clear(b);
    _barycentricfitreport_clear(rep);
    _barycentricinterpolant_init(&locb, _state, ae_true);
    _barycentricfitreport_init(&locrep, _state, ae_true);

    ae_assert(n>0, "BarycentricFitFloaterHormannWC: N<=0!", _state);
    ae_assert(m>0, "BarycentricFitFloaterHormannWC: M<=0!", _state);
    ae_assert(k>=0, "BarycentricFitFloaterHormannWC: K<0!", _state);
    ae_assert(k<m, "BarycentricFitFloaterHormannWC: K>=M!", _state);
    ae_assert(x->cnt>=n, "BarycentricFitFloaterHormannWC: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "BarycentricFitFloaterHormannWC: Length(Y)<N!", _state);
    ae_assert(w->cnt>=n, "BarycentricFitFloaterHormannWC: Length(W)<N!", _state);
    ae_assert(xc->cnt>=k, "BarycentricFitFloaterHormannWC: Length(XC)<K!", _state);
    ae_assert(yc->cnt>=k, "BarycentricFitFloaterHormannWC: Length(YC)<K!", _state);
    ae_assert(dc->cnt>=k, "BarycentricFitFloaterHormannWC: Length(DC)<K!", _state);
    ae_assert(isfinitevector(x, n, _state), "BarycentricFitFloaterHormannWC: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "BarycentricFitFloaterHormannWC: Y contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(w, n, _state), "BarycentricFitFloaterHormannWC: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(xc, k, _state), "BarycentricFitFloaterHormannWC: XC contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(yc, k, _state), "BarycentricFitFloaterHormannWC: YC contains infinite or NaN values!", _state);
    for(i=0; i<=k-1; i++)
    {
        ae_assert(dc->ptr.p_int[i]==0||dc->ptr.p_int[i]==1, "BarycentricFitFloaterHormannWC: one of DC[] is not 0 or 1!", _state);
    }
    
    /*
     * Find optimal D
     *
     * Info is -3 by default (degenerate constraints).
     * If LocInfo will always be equal to -3, Info will remain equal to -3.
     * If at least once LocInfo will be -4, Info will be -4.
     */
    wrmsbest = ae_maxrealnumber;
    rep->dbest = -1;
    *info = -3;
    for(d=0; d<=ae_minint(9, n-1, _state); d++)
    {
        lsfit_barycentricfitwcfixedd(x, y, w, n, xc, yc, dc, k, m, d, &locinfo, &locb, &locrep, _state);
        ae_assert((locinfo==-4||locinfo==-3)||locinfo>0, "BarycentricFitFloaterHormannWC: unexpected result from BarycentricFitWCFixedD!", _state);
        if( locinfo>0 )
        {
            
            /*
             * Calculate weghted RMS
             */
            wrmscur = 0;
            for(i=0; i<=n-1; i++)
            {
                wrmscur = wrmscur+ae_sqr(w->ptr.p_double[i]*(y->ptr.p_double[i]-barycentriccalc(&locb, x->ptr.p_double[i], _state)), _state);
            }
            wrmscur = ae_sqrt(wrmscur/n, _state);
            if( ae_fp_less(wrmscur,wrmsbest)||rep->dbest<0 )
            {
                barycentriccopy(&locb, b, _state);
                rep->dbest = d;
                *info = 1;
                rep->rmserror = locrep.rmserror;
                rep->avgerror = locrep.avgerror;
                rep->avgrelerror = locrep.avgrelerror;
                rep->maxerror = locrep.maxerror;
                rep->taskrcond = locrep.taskrcond;
                wrmsbest = wrmscur;
            }
        }
        else
        {
            if( locinfo!=-3&&*info<0 )
            {
                *info = locinfo;
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Rational least squares fitting using  Floater-Hormann  rational  functions
with optimal D chosen from [0,9].

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal  D  (least  root  mean
square error) is chosen.  Task  is  linear, so linear least squares solver
is used. Complexity  of  this  computational  scheme is  O(N*M^2)  (mostly
dominated by the least squares solver).

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0.
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormann(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector w;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _barycentricinterpolant_clear(b);
    _barycentricfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);

    ae_assert(n>0, "BarycentricFitFloaterHormann: N<=0!", _state);
    ae_assert(m>0, "BarycentricFitFloaterHormann: M<=0!", _state);
    ae_assert(x->cnt>=n, "BarycentricFitFloaterHormann: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "BarycentricFitFloaterHormann: Length(Y)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "BarycentricFitFloaterHormann: X contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "BarycentricFitFloaterHormann: Y contains infinite or NaN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    barycentricfitfloaterhormannwc(x, y, &w, n, &xc, &yc, &dc, 0, m, info, b, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Rational least squares fitting using  Floater-Hormann  rational  functions
with optimal D chosen from [0,9].

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal  D  (least  root  mean
square error) is chosen.  Task  is  linear, so linear least squares solver
is used. Complexity  of  this  computational  scheme is  O(N*M^2)  (mostly
dominated by the least squares solver).

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of points, N>0.
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitpenalized(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     double rho,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector w;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=1, "Spline1DFitPenalized: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitPenalized: M<4!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitPenalized: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitPenalized: Length(Y)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitPenalized: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitPenalized: Y contains infinite or NAN values!", _state);
    ae_assert(ae_isfinite(rho, _state), "Spline1DFitPenalized: Rho is infinite!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    spline1dfitpenalizedw(x, y, &w, n, m, rho, info, s, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Weighted fitting by penalized cubic spline.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are cubic splines with  natural  boundary
conditions. Problem is regularized by  adding non-linearity penalty to the
usual least squares penalty function:

    S(x) = arg min { LS + P }, where
    LS   = SUM { w[i]^2*(y[i] - S(x[i]))^2 } - least squares penalty
    P    = C*10^rho*integral{ S''(x)^2*dx } - non-linearity penalty
    rho  - tunable constant given by user
    C    - automatically determined scale parameter,
           makes penalty invariant with respect to scaling of X, Y, W.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            problem.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    M   -   number of basis functions ( = number_of_nodes), M>=4.
    Rho -   regularization  constant  passed   by   user.   It   penalizes
            nonlinearity in the regression spline. It  is  logarithmically
            scaled,  i.e.  actual  value  of  regularization  constant  is
            calculated as 10^Rho. It is automatically scaled so that:
            * Rho=2.0 corresponds to moderate amount of nonlinearity
            * generally, it should be somewhere in the [-8.0,+8.0]
            If you do not want to penalize nonlineary,
            pass small Rho. Values as low as -15 should work.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD or
                           Cholesky decomposition; problem may be
                           too ill-conditioned (very rare)
    S   -   spline interpolant.
    Rep -   Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

NOTE 1: additional nodes are added to the spline outside  of  the  fitting
interval to force linearity when x<min(x,xc) or x>max(x,xc).  It  is  done
for consistency - we penalize non-linearity  at [min(x,xc),max(x,xc)],  so
it is natural to force linearity outside of this interval.

NOTE 2: function automatically sorts points,  so  caller may pass unsorted
array.

  -- ALGLIB PROJECT --
     Copyright 19.10.2010 by Bochkanov Sergey
*************************************************************************/
void spline1dfitpenalizedw(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     ae_int_t m,
     double rho,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_int_t i;
    ae_int_t j;
    ae_int_t b;
    double v;
    double relcnt;
    double xa;
    double xb;
    double sa;
    double sb;
    ae_vector xoriginal;
    ae_vector yoriginal;
    double pdecay;
    double tdecay;
    ae_matrix fmatrix;
    ae_vector fcolumn;
    ae_vector y2;
    ae_vector w2;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;
    double fdmax;
    double admax;
    ae_matrix amatrix;
    ae_matrix d2matrix;
    double fa;
    double ga;
    double fb;
    double gb;
    double lambdav;
    ae_vector bx;
    ae_vector by;
    ae_vector bd1;
    ae_vector bd2;
    ae_vector tx;
    ae_vector ty;
    ae_vector td;
    spline1dinterpolant bs;
    ae_matrix nmatrix;
    ae_vector rightpart;
    fblslincgstate cgstate;
    ae_vector c;
    ae_vector tmp0;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_vector_init_copy(&_w, w, _state, ae_true);
    w = &_w;
    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);
    ae_vector_init(&xoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yoriginal, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&fmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&fcolumn, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&amatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&d2matrix, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&by, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bd1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bd2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ty, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&td, 0, DT_REAL, _state, ae_true);
    _spline1dinterpolant_init(&bs, _state, ae_true);
    ae_matrix_init(&nmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&rightpart, 0, DT_REAL, _state, ae_true);
    _fblslincgstate_init(&cgstate, _state, ae_true);
    ae_vector_init(&c, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp0, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=1, "Spline1DFitPenalizedW: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitPenalizedW: M<4!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitPenalizedW: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitPenalizedW: Length(Y)<N!", _state);
    ae_assert(w->cnt>=n, "Spline1DFitPenalizedW: Length(W)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitPenalizedW: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitPenalizedW: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(w, n, _state), "Spline1DFitPenalizedW: Y contains infinite or NAN values!", _state);
    ae_assert(ae_isfinite(rho, _state), "Spline1DFitPenalizedW: Rho is infinite!", _state);
    
    /*
     * Prepare LambdaV
     */
    v = -ae_log(ae_machineepsilon, _state)/ae_log(10, _state);
    if( ae_fp_less(rho,-v) )
    {
        rho = -v;
    }
    if( ae_fp_greater(rho,v) )
    {
        rho = v;
    }
    lambdav = ae_pow(10, rho, _state);
    
    /*
     * Sort X, Y, W
     */
    heapsortdpoints(x, y, w, n, _state);
    
    /*
     * Scale X, Y, XC, YC
     */
    lsfitscalexy(x, y, w, n, &xc, &yc, &dc, 0, &xa, &xb, &sa, &sb, &xoriginal, &yoriginal, _state);
    
    /*
     * Allocate space
     */
    ae_matrix_set_length(&fmatrix, n, m, _state);
    ae_matrix_set_length(&amatrix, m, m, _state);
    ae_matrix_set_length(&d2matrix, m, m, _state);
    ae_vector_set_length(&bx, m, _state);
    ae_vector_set_length(&by, m, _state);
    ae_vector_set_length(&fcolumn, n, _state);
    ae_matrix_set_length(&nmatrix, m, m, _state);
    ae_vector_set_length(&rightpart, m, _state);
    ae_vector_set_length(&tmp0, ae_maxint(m, n, _state), _state);
    ae_vector_set_length(&c, m, _state);
    
    /*
     * Fill:
     * * FMatrix by values of basis functions
     * * TmpAMatrix by second derivatives of I-th function at J-th point
     * * CMatrix by constraints
     */
    fdmax = 0;
    for(b=0; b<=m-1; b++)
    {
        
        /*
         * Prepare I-th basis function
         */
        for(j=0; j<=m-1; j++)
        {
            bx.ptr.p_double[j] = (double)(2*j)/(double)(m-1)-1;
            by.ptr.p_double[j] = 0;
        }
        by.ptr.p_double[b] = 1;
        spline1dgriddiff2cubic(&bx, &by, m, 2, 0.0, 2, 0.0, &bd1, &bd2, _state);
        spline1dbuildcubic(&bx, &by, m, 2, 0.0, 2, 0.0, &bs, _state);
        
        /*
         * Calculate B-th column of FMatrix
         * Update FDMax (maximum column norm)
         */
        spline1dconvcubic(&bx, &by, m, 2, 0.0, 2, 0.0, x, n, &fcolumn, _state);
        ae_v_move(&fmatrix.ptr.pp_double[0][b], fmatrix.stride, &fcolumn.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = 0;
        for(i=0; i<=n-1; i++)
        {
            v = v+ae_sqr(w->ptr.p_double[i]*fcolumn.ptr.p_double[i], _state);
        }
        fdmax = ae_maxreal(fdmax, v, _state);
        
        /*
         * Fill temporary with second derivatives of basis function
         */
        ae_v_move(&d2matrix.ptr.pp_double[b][0], 1, &bd2.ptr.p_double[0], 1, ae_v_len(0,m-1));
    }
    
    /*
     * * calculate penalty matrix A
     * * calculate max of diagonal elements of A
     * * calculate PDecay - coefficient before penalty matrix
     */
    for(i=0; i<=m-1; i++)
    {
        for(j=i; j<=m-1; j++)
        {
            
            /*
             * calculate integral(B_i''*B_j'') where B_i and B_j are
             * i-th and j-th basis splines.
             * B_i and B_j are piecewise linear functions.
             */
            v = 0;
            for(b=0; b<=m-2; b++)
            {
                fa = d2matrix.ptr.pp_double[i][b];
                fb = d2matrix.ptr.pp_double[i][b+1];
                ga = d2matrix.ptr.pp_double[j][b];
                gb = d2matrix.ptr.pp_double[j][b+1];
                v = v+(bx.ptr.p_double[b+1]-bx.ptr.p_double[b])*(fa*ga+(fa*(gb-ga)+ga*(fb-fa))/2+(fb-fa)*(gb-ga)/3);
            }
            amatrix.ptr.pp_double[i][j] = v;
            amatrix.ptr.pp_double[j][i] = v;
        }
    }
    admax = 0;
    for(i=0; i<=m-1; i++)
    {
        admax = ae_maxreal(admax, ae_fabs(amatrix.ptr.pp_double[i][i], _state), _state);
    }
    pdecay = lambdav*fdmax/admax;
    
    /*
     * Calculate TDecay for Tikhonov regularization
     */
    tdecay = fdmax*(1+pdecay)*10*ae_machineepsilon;
    
    /*
     * Prepare system
     *
     * NOTE: FMatrix is spoiled during this process
     */
    for(i=0; i<=n-1; i++)
    {
        v = w->ptr.p_double[i];
        ae_v_muld(&fmatrix.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
    }
    rmatrixgemm(m, m, n, 1.0, &fmatrix, 0, 0, 1, &fmatrix, 0, 0, 0, 0.0, &nmatrix, 0, 0, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=m-1; j++)
        {
            nmatrix.ptr.pp_double[i][j] = nmatrix.ptr.pp_double[i][j]+pdecay*amatrix.ptr.pp_double[i][j];
        }
    }
    for(i=0; i<=m-1; i++)
    {
        nmatrix.ptr.pp_double[i][i] = nmatrix.ptr.pp_double[i][i]+tdecay;
    }
    for(i=0; i<=m-1; i++)
    {
        rightpart.ptr.p_double[i] = 0;
    }
    for(i=0; i<=n-1; i++)
    {
        v = y->ptr.p_double[i]*w->ptr.p_double[i];
        ae_v_addd(&rightpart.ptr.p_double[0], 1, &fmatrix.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
    }
    
    /*
     * Solve system
     */
    if( !spdmatrixcholesky(&nmatrix, m, ae_true, _state) )
    {
        *info = -4;
        ae_frame_leave(_state);
        return;
    }
    fblscholeskysolve(&nmatrix, 1.0, m, ae_true, &rightpart, &tmp0, _state);
    ae_v_move(&c.ptr.p_double[0], 1, &rightpart.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * add nodes to force linearity outside of the fitting interval
     */
    spline1dgriddiffcubic(&bx, &c, m, 2, 0.0, 2, 0.0, &bd1, _state);
    ae_vector_set_length(&tx, m+2, _state);
    ae_vector_set_length(&ty, m+2, _state);
    ae_vector_set_length(&td, m+2, _state);
    ae_v_move(&tx.ptr.p_double[1], 1, &bx.ptr.p_double[0], 1, ae_v_len(1,m));
    ae_v_move(&ty.ptr.p_double[1], 1, &rightpart.ptr.p_double[0], 1, ae_v_len(1,m));
    ae_v_move(&td.ptr.p_double[1], 1, &bd1.ptr.p_double[0], 1, ae_v_len(1,m));
    tx.ptr.p_double[0] = tx.ptr.p_double[1]-(tx.ptr.p_double[2]-tx.ptr.p_double[1]);
    ty.ptr.p_double[0] = ty.ptr.p_double[1]-td.ptr.p_double[1]*(tx.ptr.p_double[2]-tx.ptr.p_double[1]);
    td.ptr.p_double[0] = td.ptr.p_double[1];
    tx.ptr.p_double[m+1] = tx.ptr.p_double[m]+(tx.ptr.p_double[m]-tx.ptr.p_double[m-1]);
    ty.ptr.p_double[m+1] = ty.ptr.p_double[m]+td.ptr.p_double[m]*(tx.ptr.p_double[m]-tx.ptr.p_double[m-1]);
    td.ptr.p_double[m+1] = td.ptr.p_double[m];
    spline1dbuildhermite(&tx, &ty, &td, m+2, s, _state);
    spline1dlintransx(s, 2/(xb-xa), -(xa+xb)/(xb-xa), _state);
    spline1dlintransy(s, sb-sa, sa, _state);
    *info = 1;
    
    /*
     * Fill report
     */
    rep->rmserror = 0;
    rep->avgerror = 0;
    rep->avgrelerror = 0;
    rep->maxerror = 0;
    relcnt = 0;
    spline1dconvcubic(&bx, &rightpart, m, 2, 0.0, 2, 0.0, x, n, &fcolumn, _state);
    for(i=0; i<=n-1; i++)
    {
        v = (sb-sa)*fcolumn.ptr.p_double[i]+sa;
        rep->rmserror = rep->rmserror+ae_sqr(v-yoriginal.ptr.p_double[i], _state);
        rep->avgerror = rep->avgerror+ae_fabs(v-yoriginal.ptr.p_double[i], _state);
        if( ae_fp_neq(yoriginal.ptr.p_double[i],0) )
        {
            rep->avgrelerror = rep->avgrelerror+ae_fabs(v-yoriginal.ptr.p_double[i], _state)/ae_fabs(yoriginal.ptr.p_double[i], _state);
            relcnt = relcnt+1;
        }
        rep->maxerror = ae_maxreal(rep->maxerror, ae_fabs(v-yoriginal.ptr.p_double[i], _state), _state);
    }
    rep->rmserror = ae_sqrt(rep->rmserror/n, _state);
    rep->avgerror = rep->avgerror/n;
    if( ae_fp_neq(relcnt,0) )
    {
        rep->avgrelerror = rep->avgrelerror/relcnt;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Weighted fitting by cubic  spline,  with constraints on function values or
derivatives.

Equidistant grid with M-2 nodes on [min(x,xc),max(x,xc)] is  used to build
basis functions. Basis functions are cubic splines with continuous  second
derivatives  and  non-fixed first  derivatives  at  interval  ends.  Small
regularizing term is used  when  solving  constrained  tasks  (to  improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitHermiteWC()  -   fitting by Hermite splines (more flexible,
                                less smooth)
    Spline1DFitCubic()      -   "lightweight" fitting  by  cubic  splines,
                                without invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints (optional):
            * 0<=K<M.
            * K=0 means no constraints (XC/YC/DC are not used)
            * if given, only first K elements of XC/YC/DC are used
            * if not given, automatically determined from XC/YC/DC
    M   -   number of basis functions ( = number_of_nodes+2), M>=4.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints  on  the  function  values  AND/OR  its
  derivatives at the interval boundaries.
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.


  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubicwc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_int_t i;

    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);

    ae_assert(n>=1, "Spline1DFitCubicWC: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitCubicWC: M<4!", _state);
    ae_assert(k>=0, "Spline1DFitCubicWC: K<0!", _state);
    ae_assert(k<m, "Spline1DFitCubicWC: K>=M!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitCubicWC: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitCubicWC: Length(Y)<N!", _state);
    ae_assert(w->cnt>=n, "Spline1DFitCubicWC: Length(W)<N!", _state);
    ae_assert(xc->cnt>=k, "Spline1DFitCubicWC: Length(XC)<K!", _state);
    ae_assert(yc->cnt>=k, "Spline1DFitCubicWC: Length(YC)<K!", _state);
    ae_assert(dc->cnt>=k, "Spline1DFitCubicWC: Length(DC)<K!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitCubicWC: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitCubicWC: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(w, n, _state), "Spline1DFitCubicWC: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(xc, k, _state), "Spline1DFitCubicWC: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(yc, k, _state), "Spline1DFitCubicWC: Y contains infinite or NAN values!", _state);
    for(i=0; i<=k-1; i++)
    {
        ae_assert(dc->ptr.p_int[i]==0||dc->ptr.p_int[i]==1, "Spline1DFitCubicWC: DC[i] is neither 0 or 1!", _state);
    }
    lsfit_spline1dfitinternal(0, x, y, w, n, xc, yc, dc, k, m, info, s, rep, _state);
}


/*************************************************************************
Weighted  fitting  by Hermite spline,  with constraints on function values
or first derivatives.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are Hermite splines.  Small  regularizing
term is used when solving constrained tasks (to improve stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitCubicWC()    -   fitting by Cubic splines (less flexible,
                                more smooth)
    Spline1DFitHermite()    -   "lightweight" Hermite fitting, without
                                invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points (optional):
            * N>0
            * if given, only first N elements of X/Y/W are processed
            * if not given, automatically determined from X/Y/W sizes
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints (optional):
            * 0<=K<M.
            * K=0 means no constraints (XC/YC/DC are not used)
            * if given, only first K elements of XC/YC/DC are used
            * if not given, automatically determined from XC/YC/DC
    M   -   number of basis functions (= 2 * number of nodes),
            M>=4,
            M IS EVEN!

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -2 means odd M was passed (which is not supported)
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

IMPORTANT:
    this subroitine supports only even M's


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the several special cases, however, we can guarantee consistency.
* one of this cases is  M>=4  and   constraints  on   the  function  value
  (AND/OR its derivative) at the interval boundaries.
* another special case is M>=4  and  ONE  constraint on the function value
  (OR, BUT NOT AND, derivative) anywhere in [min(x),max(x)]

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermitewc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_int_t i;

    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);

    ae_assert(n>=1, "Spline1DFitHermiteWC: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitHermiteWC: M<4!", _state);
    ae_assert(m%2==0, "Spline1DFitHermiteWC: M is odd!", _state);
    ae_assert(k>=0, "Spline1DFitHermiteWC: K<0!", _state);
    ae_assert(k<m, "Spline1DFitHermiteWC: K>=M!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitHermiteWC: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitHermiteWC: Length(Y)<N!", _state);
    ae_assert(w->cnt>=n, "Spline1DFitHermiteWC: Length(W)<N!", _state);
    ae_assert(xc->cnt>=k, "Spline1DFitHermiteWC: Length(XC)<K!", _state);
    ae_assert(yc->cnt>=k, "Spline1DFitHermiteWC: Length(YC)<K!", _state);
    ae_assert(dc->cnt>=k, "Spline1DFitHermiteWC: Length(DC)<K!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitHermiteWC: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitHermiteWC: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(w, n, _state), "Spline1DFitHermiteWC: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(xc, k, _state), "Spline1DFitHermiteWC: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(yc, k, _state), "Spline1DFitHermiteWC: Y contains infinite or NAN values!", _state);
    for(i=0; i<=k-1; i++)
    {
        ae_assert(dc->ptr.p_int[i]==0||dc->ptr.p_int[i]==1, "Spline1DFitHermiteWC: DC[i] is neither 0 or 1!", _state);
    }
    lsfit_spline1dfitinternal(1, x, y, w, n, xc, yc, dc, k, m, info, s, rep, _state);
}


/*************************************************************************
Least squares fitting by cubic spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitCubicWC().  See  Spline1DFitCubicWC() for more information
about subroutine parameters (we don't duplicate it here because of length)

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubic(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector w;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);

    ae_assert(n>=1, "Spline1DFitCubic: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitCubic: M<4!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitCubic: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitCubic: Length(Y)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitCubic: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitCubic: Y contains infinite or NAN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    spline1dfitcubicwc(x, y, &w, n, &xc, &yc, &dc, 0, m, info, s, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Least squares fitting by Hermite spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitHermiteWC().  See Spline1DFitHermiteWC()  description  for
more information about subroutine parameters (we don't duplicate  it  here
because of length).

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermite(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector w;
    ae_vector xc;
    ae_vector yc;
    ae_vector dc;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yc, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dc, 0, DT_INT, _state, ae_true);

    ae_assert(n>=1, "Spline1DFitHermite: N<1!", _state);
    ae_assert(m>=4, "Spline1DFitHermite: M<4!", _state);
    ae_assert(m%2==0, "Spline1DFitHermite: M is odd!", _state);
    ae_assert(x->cnt>=n, "Spline1DFitHermite: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DFitHermite: Length(Y)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "Spline1DFitHermite: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, n, _state), "Spline1DFitHermite: Y contains infinite or NAN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    spline1dfithermitewc(x, y, &w, n, &xc, &yc, &dc, 0, m, info, s, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Weighted linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -1    incorrect N/M were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearw(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{

    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);

    ae_assert(n>=1, "LSFitLinearW: N<1!", _state);
    ae_assert(m>=1, "LSFitLinearW: M<1!", _state);
    ae_assert(y->cnt>=n, "LSFitLinearW: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitLinearW: Y contains infinite or NaN values!", _state);
    ae_assert(w->cnt>=n, "LSFitLinearW: length(W)<N!", _state);
    ae_assert(isfinitevector(w, n, _state), "LSFitLinearW: W contains infinite or NaN values!", _state);
    ae_assert(fmatrix->rows>=n, "LSFitLinearW: rows(FMatrix)<N!", _state);
    ae_assert(fmatrix->cols>=m, "LSFitLinearW: cols(FMatrix)<M!", _state);
    ae_assert(apservisfinitematrix(fmatrix, n, m, _state), "LSFitLinearW: FMatrix contains infinite or NaN values!", _state);
    lsfit_lsfitlinearinternal(y, w, fmatrix, n, m, info, c, rep, _state);
}


/*************************************************************************
Weighted constained linear least squares fitting.

This  is  variation  of LSFitLinearW(), which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then LSFitLinearW()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearwc(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _y;
    ae_matrix _cmatrix;
    ae_int_t i;
    ae_int_t j;
    ae_vector tau;
    ae_matrix q;
    ae_matrix f2;
    ae_vector tmp;
    ae_vector c0;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_matrix_init_copy(&_cmatrix, cmatrix, _state, ae_true);
    cmatrix = &_cmatrix;
    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);
    ae_vector_init(&tau, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&q, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&f2, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c0, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=1, "LSFitLinearWC: N<1!", _state);
    ae_assert(m>=1, "LSFitLinearWC: M<1!", _state);
    ae_assert(k>=0, "LSFitLinearWC: K<0!", _state);
    ae_assert(y->cnt>=n, "LSFitLinearWC: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitLinearWC: Y contains infinite or NaN values!", _state);
    ae_assert(w->cnt>=n, "LSFitLinearWC: length(W)<N!", _state);
    ae_assert(isfinitevector(w, n, _state), "LSFitLinearWC: W contains infinite or NaN values!", _state);
    ae_assert(fmatrix->rows>=n, "LSFitLinearWC: rows(FMatrix)<N!", _state);
    ae_assert(fmatrix->cols>=m, "LSFitLinearWC: cols(FMatrix)<M!", _state);
    ae_assert(apservisfinitematrix(fmatrix, n, m, _state), "LSFitLinearWC: FMatrix contains infinite or NaN values!", _state);
    ae_assert(cmatrix->rows>=k, "LSFitLinearWC: rows(CMatrix)<K!", _state);
    ae_assert(cmatrix->cols>=m+1||k==0, "LSFitLinearWC: cols(CMatrix)<M+1!", _state);
    ae_assert(apservisfinitematrix(cmatrix, k, m+1, _state), "LSFitLinearWC: CMatrix contains infinite or NaN values!", _state);
    if( k>=m )
    {
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Solve
     */
    if( k==0 )
    {
        
        /*
         * no constraints
         */
        lsfit_lsfitlinearinternal(y, w, fmatrix, n, m, info, c, rep, _state);
    }
    else
    {
        
        /*
         * First, find general form solution of constraints system:
         * * factorize C = L*Q
         * * unpack Q
         * * fill upper part of C with zeros (for RCond)
         *
         * We got C=C0+Q2'*y where Q2 is lower M-K rows of Q.
         */
        rmatrixlq(cmatrix, k, m, &tau, _state);
        rmatrixlqunpackq(cmatrix, k, m, &tau, m, &q, _state);
        for(i=0; i<=k-1; i++)
        {
            for(j=i+1; j<=m-1; j++)
            {
                cmatrix->ptr.pp_double[i][j] = 0.0;
            }
        }
        if( ae_fp_less(rmatrixlurcondinf(cmatrix, k, _state),1000*ae_machineepsilon) )
        {
            *info = -3;
            ae_frame_leave(_state);
            return;
        }
        ae_vector_set_length(&tmp, k, _state);
        for(i=0; i<=k-1; i++)
        {
            if( i>0 )
            {
                v = ae_v_dotproduct(&cmatrix->ptr.pp_double[i][0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,i-1));
            }
            else
            {
                v = 0;
            }
            tmp.ptr.p_double[i] = (cmatrix->ptr.pp_double[i][m]-v)/cmatrix->ptr.pp_double[i][i];
        }
        ae_vector_set_length(&c0, m, _state);
        for(i=0; i<=m-1; i++)
        {
            c0.ptr.p_double[i] = 0;
        }
        for(i=0; i<=k-1; i++)
        {
            v = tmp.ptr.p_double[i];
            ae_v_addd(&c0.ptr.p_double[0], 1, &q.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
        }
        
        /*
         * Second, prepare modified matrix F2 = F*Q2' and solve modified task
         */
        ae_vector_set_length(&tmp, ae_maxint(n, m, _state)+1, _state);
        ae_matrix_set_length(&f2, n, m-k, _state);
        matrixvectormultiply(fmatrix, 0, n-1, 0, m-1, ae_false, &c0, 0, m-1, -1.0, y, 0, n-1, 1.0, _state);
        matrixmatrixmultiply(fmatrix, 0, n-1, 0, m-1, ae_false, &q, k, m-1, 0, m-1, ae_true, 1.0, &f2, 0, n-1, 0, m-k-1, 0.0, &tmp, _state);
        lsfit_lsfitlinearinternal(y, w, &f2, n, m-k, info, &tmp, rep, _state);
        rep->taskrcond = -1;
        if( *info<=0 )
        {
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * then, convert back to original answer: C = C0 + Q2'*Y0
         */
        ae_vector_set_length(c, m, _state);
        ae_v_move(&c->ptr.p_double[0], 1, &c0.ptr.p_double[0], 1, ae_v_len(0,m-1));
        matrixvectormultiply(&q, k, m-1, 0, m-1, ae_true, &tmp, 0, m-k-1, 1.0, c, 0, m-1, 1.0, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinear(/* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector w;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=1, "LSFitLinear: N<1!", _state);
    ae_assert(m>=1, "LSFitLinear: M<1!", _state);
    ae_assert(y->cnt>=n, "LSFitLinear: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitLinear: Y contains infinite or NaN values!", _state);
    ae_assert(fmatrix->rows>=n, "LSFitLinear: rows(FMatrix)<N!", _state);
    ae_assert(fmatrix->cols>=m, "LSFitLinear: cols(FMatrix)<M!", _state);
    ae_assert(apservisfinitematrix(fmatrix, n, m, _state), "LSFitLinear: FMatrix contains infinite or NaN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    lsfit_lsfitlinearinternal(y, &w, fmatrix, n, m, info, c, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Constained linear least squares fitting.

This  is  variation  of LSFitLinear(),  which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then  LSFitLinear()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearc(/* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* fmatrix,
     /* Real    */ ae_matrix* cmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _y;
    ae_vector w;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=1, "LSFitLinearC: N<1!", _state);
    ae_assert(m>=1, "LSFitLinearC: M<1!", _state);
    ae_assert(k>=0, "LSFitLinearC: K<0!", _state);
    ae_assert(y->cnt>=n, "LSFitLinearC: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitLinearC: Y contains infinite or NaN values!", _state);
    ae_assert(fmatrix->rows>=n, "LSFitLinearC: rows(FMatrix)<N!", _state);
    ae_assert(fmatrix->cols>=m, "LSFitLinearC: cols(FMatrix)<M!", _state);
    ae_assert(apservisfinitematrix(fmatrix, n, m, _state), "LSFitLinearC: FMatrix contains infinite or NaN values!", _state);
    ae_assert(cmatrix->rows>=k, "LSFitLinearC: rows(CMatrix)<K!", _state);
    ae_assert(cmatrix->cols>=m+1||k==0, "LSFitLinearC: cols(CMatrix)<M+1!", _state);
    ae_assert(apservisfinitematrix(cmatrix, k, m+1, _state), "LSFitLinearC: CMatrix contains infinite or NaN values!", _state);
    ae_vector_set_length(&w, n, _state);
    for(i=0; i<=n-1; i++)
    {
        w.ptr.p_double[i] = 1;
    }
    lsfitlinearwc(y, &w, fmatrix, cmatrix, n, m, k, info, c, rep, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Weighted nonlinear least squares fitting using function values only.

Combination of numerical differentiation and secant updates is used to
obtain function Jacobian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]).

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    DiffStep-   numerical differentiation step;
                should not be very small or large;
                large = loss of accuracy
                small = growth of round-off errors

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 18.10.2008 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewf(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     double diffstep,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateWF: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateWF: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateWF: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateWF: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateWF: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateWF: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateWF: Y contains infinite or NaN values!", _state);
    ae_assert(w->cnt>=n, "LSFitCreateWF: length(W)<N!", _state);
    ae_assert(isfinitevector(w, n, _state), "LSFitCreateWF: W contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateWF: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateWF: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateWF: X contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(diffstep, _state), "LSFitCreateWF: DiffStep is not finite!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "LSFitCreateWF: DiffStep<=0!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->w.ptr.p_double[0], 1, &w->ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
    }
    minlmcreatev(k, n, &state->c, diffstep, &state->optstate, _state);
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Nonlinear least squares fitting using function values only.

Combination of numerical differentiation and secant updates is used to
obtain function Jacobian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (f(c,x[0])-y[0])^2 + ... + (f(c,x[n-1])-y[n-1])^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]).

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    DiffStep-   numerical differentiation step;
                should not be very small or large;
                large = loss of accuracy
                small = growth of round-off errors

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 18.10.2008 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatef(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     double diffstep,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateFG: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateFG: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateFG: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateFG: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateFG: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateFG: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateFG: Y contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateFG: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateFG: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateFG: X contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateFG: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateFG: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateFG: X contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(diffstep, _state), "LSFitCreateWF: DiffStep is not finite!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "LSFitCreateWF: DiffStep<=0!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
        state->w.ptr.p_double[i] = 1;
    }
    minlmcreatev(k, n, &state->c, diffstep, &state->optstate, _state);
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Weighted nonlinear least squares fitting using gradient only.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,
    
    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted
    
This subroutine uses only f(c,x[i]) and its gradient.
    
INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also:
    LSFitResults
    LSFitCreateFG (fitting without weights)
    LSFitCreateWFGH (fitting using Hessian)
    LSFitCreateFGH (fitting using Hessian, without weights)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewfg(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_bool cheapfg,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateWFG: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateWFG: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateWFG: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateWFG: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateWFG: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateWFG: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateWFG: Y contains infinite or NaN values!", _state);
    ae_assert(w->cnt>=n, "LSFitCreateWFG: length(W)<N!", _state);
    ae_assert(isfinitevector(w, n, _state), "LSFitCreateWFG: W contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateWFG: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateWFG: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateWFG: X contains infinite or NaN values!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_vector_set_length(&state->g, k, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->w.ptr.p_double[0], 1, &w->ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
    }
    if( cheapfg )
    {
        minlmcreatevgj(k, n, &state->c, &state->optstate, _state);
    }
    else
    {
        minlmcreatevj(k, n, &state->c, &state->optstate, _state);
    }
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Nonlinear least squares fitting using gradient only, without individual
weights.

Nonlinear task min(F(c)) is solved, where

    F(c) = ((f(c,x[0])-y[0]))^2 + ... + ((f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses only f(c,x[i]) and its gradient.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatefg(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     ae_bool cheapfg,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateFG: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateFG: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateFG: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateFG: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateFG: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateFG: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateFG: Y contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateFG: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateFG: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateFG: X contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateFG: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateFG: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateFG: X contains infinite or NaN values!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_vector_set_length(&state->g, k, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
        state->w.ptr.p_double[i] = 1;
    }
    if( cheapfg )
    {
        minlmcreatevgj(k, n, &state->c, &state->optstate, _state);
    }
    else
    {
        minlmcreatevj(k, n, &state->c, &state->optstate, _state);
    }
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Weighted nonlinear least squares fitting using gradient/Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(c,x[i]), its gradient and its Hessian.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatewfgh(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateWFGH: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateWFGH: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateWFGH: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateWFGH: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateWFGH: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateWFGH: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateWFGH: Y contains infinite or NaN values!", _state);
    ae_assert(w->cnt>=n, "LSFitCreateWFGH: length(W)<N!", _state);
    ae_assert(isfinitevector(w, n, _state), "LSFitCreateWFGH: W contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateWFGH: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateWFGH: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateWFGH: X contains infinite or NaN values!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_matrix_set_length(&state->h, k, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_vector_set_length(&state->g, k, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->w.ptr.p_double[0], 1, &w->ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
    }
    minlmcreatefgh(k, &state->c, &state->optstate, _state);
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Nonlinear least squares fitting using gradient/Hessian, without individial
weights.

Nonlinear task min(F(c)) is solved, where

    F(c) = ((f(c,x[0])-y[0]))^2 + ... + ((f(c,x[n-1])-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(c,x[i]), its gradient and its Hessian.

INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitcreatefgh(/* Real    */ ae_matrix* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* c,
     ae_int_t n,
     ae_int_t m,
     ae_int_t k,
     lsfitstate* state,
     ae_state *_state)
{
    ae_int_t i;

    _lsfitstate_clear(state);

    ae_assert(n>=1, "LSFitCreateFGH: N<1!", _state);
    ae_assert(m>=1, "LSFitCreateFGH: M<1!", _state);
    ae_assert(k>=1, "LSFitCreateFGH: K<1!", _state);
    ae_assert(c->cnt>=k, "LSFitCreateFGH: length(C)<K!", _state);
    ae_assert(isfinitevector(c, k, _state), "LSFitCreateFGH: C contains infinite or NaN values!", _state);
    ae_assert(y->cnt>=n, "LSFitCreateFGH: length(Y)<N!", _state);
    ae_assert(isfinitevector(y, n, _state), "LSFitCreateFGH: Y contains infinite or NaN values!", _state);
    ae_assert(x->rows>=n, "LSFitCreateFGH: rows(X)<N!", _state);
    ae_assert(x->cols>=m, "LSFitCreateFGH: cols(X)<M!", _state);
    ae_assert(apservisfinitematrix(x, n, m, _state), "LSFitCreateFGH: X contains infinite or NaN values!", _state);
    state->n = n;
    state->m = m;
    state->k = k;
    lsfitsetcond(state, 0.0, 0.0, 0, _state);
    lsfitsetstpmax(state, 0.0, _state);
    lsfitsetxrep(state, ae_false, _state);
    ae_matrix_set_length(&state->taskx, n, m, _state);
    ae_vector_set_length(&state->tasky, n, _state);
    ae_vector_set_length(&state->w, n, _state);
    ae_vector_set_length(&state->c, k, _state);
    ae_matrix_set_length(&state->h, k, k, _state);
    ae_vector_set_length(&state->x, m, _state);
    ae_vector_set_length(&state->g, k, _state);
    ae_v_move(&state->c.ptr.p_double[0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,k-1));
    for(i=0; i<=n-1; i++)
    {
        ae_v_move(&state->taskx.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
        state->tasky.ptr.p_double[i] = y->ptr.p_double[i];
        state->w.ptr.p_double[i] = 1;
    }
    minlmcreatefgh(k, &state->c, &state->optstate, _state);
    lsfit_lsfitclearrequestfields(state, _state);
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
Stopping conditions for nonlinear least squares fitting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (according to the scheme used by MINLM unit).


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitsetcond(lsfitstate* state,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsf, _state), "LSFitSetCond: EpsF is not finite!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "LSFitSetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "LSFitSetCond: EpsX is not finite!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "LSFitSetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "LSFitSetCond: negative MaxIts!", _state);
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void lsfitsetstpmax(lsfitstate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_fp_greater_eq(stpmax,0), "LSFitSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not
    
When reports are needed, State.C (current parameters) and State.F (current
value of fitting function) are reported.


  -- ALGLIB --
     Copyright 15.08.2010 by Bochkanov Sergey
*************************************************************************/
void lsfitsetxrep(lsfitstate* state, ae_bool needxrep, ae_state *_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
NOTES:

1. this algorithm is somewhat unusual because it works with  parameterized
   function f(C,X), where X is a function argument (we  have  many  points
   which are characterized by different  argument  values),  and  C  is  a
   parameter to fit.

   For example, if we want to do linear fit by f(c0,c1,x) = c0*x+c1,  then
   x will be argument, and {c0,c1} will be parameters.
   
   It is important to understand that this algorithm finds minimum in  the
   space of function PARAMETERS (not arguments), so it  needs  derivatives
   of f() with respect to C, not X.
   
   In the example above it will need f=c0*x+c1 and {df/dc0,df/dc1} = {x,1}
   instead of {df/dx} = {c0}.

2. Callback functions accept C as the first parameter, and X as the second

3. If  state  was  created  with  LSFitCreateFG(),  algorithm  needs  just
   function   and   its   gradient,   but   if   state   was  created with
   LSFitCreateFGH(), algorithm will need function, gradient and Hessian.
   
   According  to  the  said  above,  there  ase  several  versions of this
   function, which accept different sets of callbacks.
   
   This flexibility opens way to subtle errors - you may create state with
   LSFitCreateFGH() (optimization using Hessian), but call function  which
   does not accept Hessian. So when algorithm will request Hessian,  there
   will be no callback to call. In this case exception will be thrown.
   
   Be careful to avoid such errors because there is no way to find them at
   compile time - you can see them at runtime only.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool lsfititeration(lsfitstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t k;
    ae_int_t i;
    ae_int_t j;
    double v;
    double relcnt;
    ae_bool result;


    
    /*
     * Reverse communication preparations
     * I know it looks ugly, but it works the same way
     * anywhere from C++ to Python.
     *
     * This code initializes locals by:
     * * random values determined during code
     *   generation - on first subroutine call
     * * values from previous call - on subsequent calls
     */
    if( state->rstate.stage>=0 )
    {
        n = state->rstate.ia.ptr.p_int[0];
        m = state->rstate.ia.ptr.p_int[1];
        k = state->rstate.ia.ptr.p_int[2];
        i = state->rstate.ia.ptr.p_int[3];
        j = state->rstate.ia.ptr.p_int[4];
        v = state->rstate.ra.ptr.p_double[0];
        relcnt = state->rstate.ra.ptr.p_double[1];
    }
    else
    {
        n = -983;
        m = -989;
        k = -834;
        i = 900;
        j = -287;
        v = 364;
        relcnt = 214;
    }
    if( state->rstate.stage==0 )
    {
        goto lbl_0;
    }
    if( state->rstate.stage==1 )
    {
        goto lbl_1;
    }
    if( state->rstate.stage==2 )
    {
        goto lbl_2;
    }
    if( state->rstate.stage==3 )
    {
        goto lbl_3;
    }
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    
    /*
     * Routine body
     */
    
    /*
     * init
     */
    n = state->n;
    m = state->m;
    k = state->k;
    minlmsetcond(&state->optstate, 0.0, state->epsf, state->epsx, state->maxits, _state);
    minlmsetstpmax(&state->optstate, state->stpmax, _state);
    minlmsetxrep(&state->optstate, state->xrep, _state);
    
    /*
     * Optimize
     */
lbl_7:
    if( !minlmiteration(&state->optstate, _state) )
    {
        goto lbl_8;
    }
    if( !state->optstate.needfi )
    {
        goto lbl_9;
    }
    
    /*
     * calculate f[] = wi*(f(xi,c)-yi)
     */
    i = 0;
lbl_11:
    if( i>n-1 )
    {
        goto lbl_13;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needf = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needf = ae_false;
    state->optstate.fi.ptr.p_double[i] = state->w.ptr.p_double[i]*(state->f-state->tasky.ptr.p_double[i]);
    i = i+1;
    goto lbl_11;
lbl_13:
    goto lbl_7;
lbl_9:
    if( !state->optstate.needf )
    {
        goto lbl_14;
    }
    
    /*
     * calculate F = sum (wi*(f(xi,c)-yi))^2
     */
    state->optstate.f = 0;
    i = 0;
lbl_16:
    if( i>n-1 )
    {
        goto lbl_18;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needf = ae_true;
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->needf = ae_false;
    state->optstate.f = state->optstate.f+ae_sqr(state->w.ptr.p_double[i]*(state->f-state->tasky.ptr.p_double[i]), _state);
    i = i+1;
    goto lbl_16;
lbl_18:
    goto lbl_7;
lbl_14:
    if( !state->optstate.needfg )
    {
        goto lbl_19;
    }
    
    /*
     * calculate F/gradF
     */
    state->optstate.f = 0;
    for(i=0; i<=k-1; i++)
    {
        state->optstate.g.ptr.p_double[i] = 0;
    }
    i = 0;
lbl_21:
    if( i>n-1 )
    {
        goto lbl_23;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->needfg = ae_false;
    state->optstate.f = state->optstate.f+ae_sqr(state->w.ptr.p_double[i]*(state->f-state->tasky.ptr.p_double[i]), _state);
    v = ae_sqr(state->w.ptr.p_double[i], _state)*2*(state->f-state->tasky.ptr.p_double[i]);
    ae_v_addd(&state->optstate.g.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,k-1), v);
    i = i+1;
    goto lbl_21;
lbl_23:
    goto lbl_7;
lbl_19:
    if( !state->optstate.needfij )
    {
        goto lbl_24;
    }
    
    /*
     * calculate Fi/jac(Fi)
     */
    i = 0;
lbl_26:
    if( i>n-1 )
    {
        goto lbl_28;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    state->optstate.fi.ptr.p_double[i] = state->w.ptr.p_double[i]*(state->f-state->tasky.ptr.p_double[i]);
    v = state->w.ptr.p_double[i];
    ae_v_moved(&state->optstate.j.ptr.pp_double[i][0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,k-1), v);
    i = i+1;
    goto lbl_26;
lbl_28:
    goto lbl_7;
lbl_24:
    if( !state->optstate.needfgh )
    {
        goto lbl_29;
    }
    
    /*
     * calculate F/grad(F)/hess(F)
     */
    state->optstate.f = 0;
    for(i=0; i<=k-1; i++)
    {
        state->optstate.g.ptr.p_double[i] = 0;
    }
    for(i=0; i<=k-1; i++)
    {
        for(j=0; j<=k-1; j++)
        {
            state->optstate.h.ptr.pp_double[i][j] = 0;
        }
    }
    i = 0;
lbl_31:
    if( i>n-1 )
    {
        goto lbl_33;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needfgh = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->needfgh = ae_false;
    state->optstate.f = state->optstate.f+ae_sqr(state->w.ptr.p_double[i]*(state->f-state->tasky.ptr.p_double[i]), _state);
    v = ae_sqr(state->w.ptr.p_double[i], _state)*2*(state->f-state->tasky.ptr.p_double[i]);
    ae_v_addd(&state->optstate.g.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,k-1), v);
    for(j=0; j<=k-1; j++)
    {
        v = 2*ae_sqr(state->w.ptr.p_double[i], _state)*state->g.ptr.p_double[j];
        ae_v_addd(&state->optstate.h.ptr.pp_double[j][0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,k-1), v);
        v = 2*ae_sqr(state->w.ptr.p_double[i], _state)*(state->f-state->tasky.ptr.p_double[i]);
        ae_v_addd(&state->optstate.h.ptr.pp_double[j][0], 1, &state->h.ptr.pp_double[j][0], 1, ae_v_len(0,k-1), v);
    }
    i = i+1;
    goto lbl_31;
lbl_33:
    goto lbl_7;
lbl_29:
    if( !state->optstate.xupdated )
    {
        goto lbl_34;
    }
    
    /*
     * Report new iteration
     */
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->optstate.x.ptr.p_double[0], 1, ae_v_len(0,k-1));
    state->f = state->optstate.f;
    lsfit_lsfitclearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->xupdated = ae_false;
    goto lbl_7;
lbl_34:
    goto lbl_7;
lbl_8:
    minlmresults(&state->optstate, &state->c, &state->optrep, _state);
    state->repterminationtype = state->optrep.terminationtype;
    
    /*
     * calculate errors
     */
    if( state->repterminationtype<=0 )
    {
        goto lbl_36;
    }
    state->reprmserror = 0;
    state->repavgerror = 0;
    state->repavgrelerror = 0;
    state->repmaxerror = 0;
    relcnt = 0;
    i = 0;
lbl_38:
    if( i>n-1 )
    {
        goto lbl_40;
    }
    ae_v_move(&state->c.ptr.p_double[0], 1, &state->c.ptr.p_double[0], 1, ae_v_len(0,k-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->taskx.ptr.pp_double[i][0], 1, ae_v_len(0,m-1));
    state->pointindex = i;
    lsfit_lsfitclearrequestfields(state, _state);
    state->needf = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->needf = ae_false;
    v = state->f;
    state->reprmserror = state->reprmserror+ae_sqr(v-state->tasky.ptr.p_double[i], _state);
    state->repavgerror = state->repavgerror+ae_fabs(v-state->tasky.ptr.p_double[i], _state);
    if( ae_fp_neq(state->tasky.ptr.p_double[i],0) )
    {
        state->repavgrelerror = state->repavgrelerror+ae_fabs(v-state->tasky.ptr.p_double[i], _state)/ae_fabs(state->tasky.ptr.p_double[i], _state);
        relcnt = relcnt+1;
    }
    state->repmaxerror = ae_maxreal(state->repmaxerror, ae_fabs(v-state->tasky.ptr.p_double[i], _state), _state);
    i = i+1;
    goto lbl_38;
lbl_40:
    state->reprmserror = ae_sqrt(state->reprmserror/n, _state);
    state->repavgerror = state->repavgerror/n;
    if( ae_fp_neq(relcnt,0) )
    {
        state->repavgrelerror = state->repavgrelerror/relcnt;
    }
lbl_36:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = k;
    state->rstate.ia.ptr.p_int[3] = i;
    state->rstate.ia.ptr.p_int[4] = j;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = relcnt;
    return result;
}


/*************************************************************************
Nonlinear least squares fitting results.

Called after return from LSFitFit().

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    Info    -   completetion code:
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
    C       -   array[0..K-1], solution
    Rep     -   optimization report. Following fields are set:
                * Rep.TerminationType completetion code:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitresults(lsfitstate* state,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{

    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);

    *info = state->repterminationtype;
    if( *info>0 )
    {
        ae_vector_set_length(c, state->k, _state);
        ae_v_move(&c->ptr.p_double[0], 1, &state->c.ptr.p_double[0], 1, ae_v_len(0,state->k-1));
        rep->rmserror = state->reprmserror;
        rep->avgerror = state->repavgerror;
        rep->avgrelerror = state->repavgrelerror;
        rep->maxerror = state->repmaxerror;
    }
}


/*************************************************************************
Internal subroutine: automatic scaling for LLS tasks.
NEVER CALL IT DIRECTLY!

Maps abscissas to [-1,1], standartizes ordinates and correspondingly scales
constraints. It also scales weights so that max(W[i])=1

Transformations performed:
* X, XC         [XA,XB] => [-1,+1]
                transformation makes min(X)=-1, max(X)=+1

* Y             [SA,SB] => [0,1]
                transformation makes mean(Y)=0, stddev(Y)=1
                
* YC            transformed accordingly to SA, SB, DC[I]

  -- ALGLIB PROJECT --
     Copyright 08.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitscalexy(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     double* xa,
     double* xb,
     double* sa,
     double* sb,
     /* Real    */ ae_vector* xoriginal,
     /* Real    */ ae_vector* yoriginal,
     ae_state *_state)
{
    double xmin;
    double xmax;
    ae_int_t i;
    double mx;

    *xa = 0;
    *xb = 0;
    *sa = 0;
    *sb = 0;
    ae_vector_clear(xoriginal);
    ae_vector_clear(yoriginal);

    ae_assert(n>=1, "LSFitScaleXY: incorrect N", _state);
    ae_assert(k>=0, "LSFitScaleXY: incorrect K", _state);
    
    /*
     * Calculate xmin/xmax.
     * Force xmin<>xmax.
     */
    xmin = x->ptr.p_double[0];
    xmax = x->ptr.p_double[0];
    for(i=1; i<=n-1; i++)
    {
        xmin = ae_minreal(xmin, x->ptr.p_double[i], _state);
        xmax = ae_maxreal(xmax, x->ptr.p_double[i], _state);
    }
    for(i=0; i<=k-1; i++)
    {
        xmin = ae_minreal(xmin, xc->ptr.p_double[i], _state);
        xmax = ae_maxreal(xmax, xc->ptr.p_double[i], _state);
    }
    if( ae_fp_eq(xmin,xmax) )
    {
        if( ae_fp_eq(xmin,0) )
        {
            xmin = -1;
            xmax = 1;
        }
        else
        {
            if( ae_fp_greater(xmin,0) )
            {
                xmin = 0.5*xmin;
            }
            else
            {
                xmax = 0.5*xmax;
            }
        }
    }
    
    /*
     * Transform abscissas: map [XA,XB] to [0,1]
     *
     * Store old X[] in XOriginal[] (it will be used
     * to calculate relative error).
     */
    ae_vector_set_length(xoriginal, n, _state);
    ae_v_move(&xoriginal->ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    *xa = xmin;
    *xb = xmax;
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[i] = 2*(x->ptr.p_double[i]-0.5*(*xa+(*xb)))/(*xb-(*xa));
    }
    for(i=0; i<=k-1; i++)
    {
        ae_assert(dc->ptr.p_int[i]>=0, "LSFitScaleXY: internal error!", _state);
        xc->ptr.p_double[i] = 2*(xc->ptr.p_double[i]-0.5*(*xa+(*xb)))/(*xb-(*xa));
        yc->ptr.p_double[i] = yc->ptr.p_double[i]*ae_pow(0.5*(*xb-(*xa)), dc->ptr.p_int[i], _state);
    }
    
    /*
     * Transform function values: map [SA,SB] to [0,1]
     * SA = mean(Y),
     * SB = SA+stddev(Y).
     *
     * Store old Y[] in YOriginal[] (it will be used
     * to calculate relative error).
     */
    ae_vector_set_length(yoriginal, n, _state);
    ae_v_move(&yoriginal->ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
    *sa = 0;
    for(i=0; i<=n-1; i++)
    {
        *sa = *sa+y->ptr.p_double[i];
    }
    *sa = *sa/n;
    *sb = 0;
    for(i=0; i<=n-1; i++)
    {
        *sb = *sb+ae_sqr(y->ptr.p_double[i]-(*sa), _state);
    }
    *sb = ae_sqrt(*sb/n, _state)+(*sa);
    if( ae_fp_eq(*sb,*sa) )
    {
        *sb = 2*(*sa);
    }
    if( ae_fp_eq(*sb,*sa) )
    {
        *sb = *sa+1;
    }
    for(i=0; i<=n-1; i++)
    {
        y->ptr.p_double[i] = (y->ptr.p_double[i]-(*sa))/(*sb-(*sa));
    }
    for(i=0; i<=k-1; i++)
    {
        if( dc->ptr.p_int[i]==0 )
        {
            yc->ptr.p_double[i] = (yc->ptr.p_double[i]-(*sa))/(*sb-(*sa));
        }
        else
        {
            yc->ptr.p_double[i] = yc->ptr.p_double[i]/(*sb-(*sa));
        }
    }
    
    /*
     * Scale weights
     */
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        mx = ae_maxreal(mx, ae_fabs(w->ptr.p_double[i], _state), _state);
    }
    if( ae_fp_neq(mx,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            w->ptr.p_double[i] = w->ptr.p_double[i]/mx;
        }
    }
}


/*************************************************************************
Internal spline fitting subroutine

  -- ALGLIB PROJECT --
     Copyright 08.09.2009 by Bochkanov Sergey
*************************************************************************/
static void lsfit_spline1dfitinternal(ae_int_t st,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     spline1dinterpolant* s,
     spline1dfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_matrix fmatrix;
    ae_matrix cmatrix;
    ae_vector y2;
    ae_vector w2;
    ae_vector sx;
    ae_vector sy;
    ae_vector sd;
    ae_vector tmp;
    ae_vector xoriginal;
    ae_vector yoriginal;
    lsfitreport lrep;
    double v0;
    double v1;
    double v2;
    double mx;
    spline1dinterpolant s2;
    ae_int_t i;
    ae_int_t j;
    ae_int_t relcnt;
    double xa;
    double xb;
    double sa;
    double sb;
    double bl;
    double br;
    double decay;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_vector_init_copy(&_w, w, _state, ae_true);
    w = &_w;
    ae_vector_init_copy(&_xc, xc, _state, ae_true);
    xc = &_xc;
    ae_vector_init_copy(&_yc, yc, _state, ae_true);
    yc = &_yc;
    *info = 0;
    _spline1dinterpolant_clear(s);
    _spline1dfitreport_clear(rep);
    ae_matrix_init(&fmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sd, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yoriginal, 0, DT_REAL, _state, ae_true);
    _lsfitreport_init(&lrep, _state, ae_true);
    _spline1dinterpolant_init(&s2, _state, ae_true);

    ae_assert(st==0||st==1, "Spline1DFit: internal error!", _state);
    if( st==0&&m<4 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( st==1&&m<4 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( (n<1||k<0)||k>=m )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=k-1; i++)
    {
        *info = 0;
        if( dc->ptr.p_int[i]<0 )
        {
            *info = -1;
        }
        if( dc->ptr.p_int[i]>1 )
        {
            *info = -1;
        }
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
    }
    if( st==1&&m%2!=0 )
    {
        
        /*
         * Hermite fitter must have even number of basis functions
         */
        *info = -2;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * weight decay for correct handling of task which becomes
     * degenerate after constraints are applied
     */
    decay = 10000*ae_machineepsilon;
    
    /*
     * Scale X, Y, XC, YC
     */
    lsfitscalexy(x, y, w, n, xc, yc, dc, k, &xa, &xb, &sa, &sb, &xoriginal, &yoriginal, _state);
    
    /*
     * allocate space, initialize:
     * * SX     -   grid for basis functions
     * * SY     -   values of basis functions at grid points
     * * FMatrix-   values of basis functions at X[]
     * * CMatrix-   values (derivatives) of basis functions at XC[]
     */
    ae_vector_set_length(&y2, n+m, _state);
    ae_vector_set_length(&w2, n+m, _state);
    ae_matrix_set_length(&fmatrix, n+m, m, _state);
    if( k>0 )
    {
        ae_matrix_set_length(&cmatrix, k, m+1, _state);
    }
    if( st==0 )
    {
        
        /*
         * allocate space for cubic spline
         */
        ae_vector_set_length(&sx, m-2, _state);
        ae_vector_set_length(&sy, m-2, _state);
        for(j=0; j<=m-2-1; j++)
        {
            sx.ptr.p_double[j] = (double)(2*j)/(double)(m-2-1)-1;
        }
    }
    if( st==1 )
    {
        
        /*
         * allocate space for Hermite spline
         */
        ae_vector_set_length(&sx, m/2, _state);
        ae_vector_set_length(&sy, m/2, _state);
        ae_vector_set_length(&sd, m/2, _state);
        for(j=0; j<=m/2-1; j++)
        {
            sx.ptr.p_double[j] = (double)(2*j)/(double)(m/2-1)-1;
        }
    }
    
    /*
     * Prepare design and constraints matrices:
     * * fill constraints matrix
     * * fill first N rows of design matrix with values
     * * fill next M rows of design matrix with regularizing term
     * * append M zeros to Y
     * * append M elements, mean(abs(W)) each, to W
     */
    for(j=0; j<=m-1; j++)
    {
        
        /*
         * prepare Jth basis function
         */
        if( st==0 )
        {
            
            /*
             * cubic spline basis
             */
            for(i=0; i<=m-2-1; i++)
            {
                sy.ptr.p_double[i] = 0;
            }
            bl = 0;
            br = 0;
            if( j<m-2 )
            {
                sy.ptr.p_double[j] = 1;
            }
            if( j==m-2 )
            {
                bl = 1;
            }
            if( j==m-1 )
            {
                br = 1;
            }
            spline1dbuildcubic(&sx, &sy, m-2, 1, bl, 1, br, &s2, _state);
        }
        if( st==1 )
        {
            
            /*
             * Hermite basis
             */
            for(i=0; i<=m/2-1; i++)
            {
                sy.ptr.p_double[i] = 0;
                sd.ptr.p_double[i] = 0;
            }
            if( j%2==0 )
            {
                sy.ptr.p_double[j/2] = 1;
            }
            else
            {
                sd.ptr.p_double[j/2] = 1;
            }
            spline1dbuildhermite(&sx, &sy, &sd, m/2, &s2, _state);
        }
        
        /*
         * values at X[], XC[]
         */
        for(i=0; i<=n-1; i++)
        {
            fmatrix.ptr.pp_double[i][j] = spline1dcalc(&s2, x->ptr.p_double[i], _state);
        }
        for(i=0; i<=k-1; i++)
        {
            ae_assert(dc->ptr.p_int[i]>=0&&dc->ptr.p_int[i]<=2, "Spline1DFit: internal error!", _state);
            spline1ddiff(&s2, xc->ptr.p_double[i], &v0, &v1, &v2, _state);
            if( dc->ptr.p_int[i]==0 )
            {
                cmatrix.ptr.pp_double[i][j] = v0;
            }
            if( dc->ptr.p_int[i]==1 )
            {
                cmatrix.ptr.pp_double[i][j] = v1;
            }
            if( dc->ptr.p_int[i]==2 )
            {
                cmatrix.ptr.pp_double[i][j] = v2;
            }
        }
    }
    for(i=0; i<=k-1; i++)
    {
        cmatrix.ptr.pp_double[i][m] = yc->ptr.p_double[i];
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=m-1; j++)
        {
            if( i==j )
            {
                fmatrix.ptr.pp_double[n+i][j] = decay;
            }
            else
            {
                fmatrix.ptr.pp_double[n+i][j] = 0;
            }
        }
    }
    ae_vector_set_length(&y2, n+m, _state);
    ae_vector_set_length(&w2, n+m, _state);
    ae_v_move(&y2.ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&w2.ptr.p_double[0], 1, &w->ptr.p_double[0], 1, ae_v_len(0,n-1));
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        mx = mx+ae_fabs(w->ptr.p_double[i], _state);
    }
    mx = mx/n;
    for(i=0; i<=m-1; i++)
    {
        y2.ptr.p_double[n+i] = 0;
        w2.ptr.p_double[n+i] = mx;
    }
    
    /*
     * Solve constrained task
     */
    if( k>0 )
    {
        
        /*
         * solve using regularization
         */
        lsfitlinearwc(&y2, &w2, &fmatrix, &cmatrix, n+m, m, k, info, &tmp, &lrep, _state);
    }
    else
    {
        
        /*
         * no constraints, no regularization needed
         */
        lsfitlinearwc(y, w, &fmatrix, &cmatrix, n, m, k, info, &tmp, &lrep, _state);
    }
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Generate spline and scale it
     */
    if( st==0 )
    {
        
        /*
         * cubic spline basis
         */
        ae_v_move(&sy.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,m-2-1));
        spline1dbuildcubic(&sx, &sy, m-2, 1, tmp.ptr.p_double[m-2], 1, tmp.ptr.p_double[m-1], s, _state);
    }
    if( st==1 )
    {
        
        /*
         * Hermite basis
         */
        for(i=0; i<=m/2-1; i++)
        {
            sy.ptr.p_double[i] = tmp.ptr.p_double[2*i];
            sd.ptr.p_double[i] = tmp.ptr.p_double[2*i+1];
        }
        spline1dbuildhermite(&sx, &sy, &sd, m/2, s, _state);
    }
    spline1dlintransx(s, 2/(xb-xa), -(xa+xb)/(xb-xa), _state);
    spline1dlintransy(s, sb-sa, sa, _state);
    
    /*
     * Scale absolute errors obtained from LSFitLinearW.
     * Relative error should be calculated separately
     * (because of shifting/scaling of the task)
     */
    rep->taskrcond = lrep.taskrcond;
    rep->rmserror = lrep.rmserror*(sb-sa);
    rep->avgerror = lrep.avgerror*(sb-sa);
    rep->maxerror = lrep.maxerror*(sb-sa);
    rep->avgrelerror = 0;
    relcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(yoriginal.ptr.p_double[i],0) )
        {
            rep->avgrelerror = rep->avgrelerror+ae_fabs(spline1dcalc(s, xoriginal.ptr.p_double[i], _state)-yoriginal.ptr.p_double[i], _state)/ae_fabs(yoriginal.ptr.p_double[i], _state);
            relcnt = relcnt+1;
        }
    }
    if( relcnt!=0 )
    {
        rep->avgrelerror = rep->avgrelerror/relcnt;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal fitting subroutine
*************************************************************************/
static void lsfit_lsfitlinearinternal(/* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     /* Real    */ ae_matrix* fmatrix,
     ae_int_t n,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    double threshold;
    ae_matrix ft;
    ae_matrix q;
    ae_matrix l;
    ae_matrix r;
    ae_vector b;
    ae_vector wmod;
    ae_vector tau;
    ae_int_t i;
    ae_int_t j;
    double v;
    ae_vector sv;
    ae_matrix u;
    ae_matrix vt;
    ae_vector tmp;
    ae_vector utb;
    ae_vector sutb;
    ae_int_t relcnt;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);
    ae_matrix_init(&ft, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&q, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&l, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&r, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wmod, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tau, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sv, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&utb, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sutb, 0, DT_REAL, _state, ae_true);

    if( n<1||m<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    threshold = ae_sqrt(ae_machineepsilon, _state);
    
    /*
     * Degenerate case, needs special handling
     */
    if( n<m )
    {
        
        /*
         * Create design matrix.
         */
        ae_matrix_set_length(&ft, n, m, _state);
        ae_vector_set_length(&b, n, _state);
        ae_vector_set_length(&wmod, n, _state);
        for(j=0; j<=n-1; j++)
        {
            v = w->ptr.p_double[j];
            ae_v_moved(&ft.ptr.pp_double[j][0], 1, &fmatrix->ptr.pp_double[j][0], 1, ae_v_len(0,m-1), v);
            b.ptr.p_double[j] = w->ptr.p_double[j]*y->ptr.p_double[j];
            wmod.ptr.p_double[j] = 1;
        }
        
        /*
         * LQ decomposition and reduction to M=N
         */
        ae_vector_set_length(c, m, _state);
        for(i=0; i<=m-1; i++)
        {
            c->ptr.p_double[i] = 0;
        }
        rep->taskrcond = 0;
        rmatrixlq(&ft, n, m, &tau, _state);
        rmatrixlqunpackq(&ft, n, m, &tau, n, &q, _state);
        rmatrixlqunpackl(&ft, n, m, &l, _state);
        lsfit_lsfitlinearinternal(&b, &wmod, &l, n, n, info, &tmp, rep, _state);
        if( *info<=0 )
        {
            ae_frame_leave(_state);
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            v = tmp.ptr.p_double[i];
            ae_v_addd(&c->ptr.p_double[0], 1, &q.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * N>=M. Generate design matrix and reduce to N=M using
     * QR decomposition.
     */
    ae_matrix_set_length(&ft, n, m, _state);
    ae_vector_set_length(&b, n, _state);
    for(j=0; j<=n-1; j++)
    {
        v = w->ptr.p_double[j];
        ae_v_moved(&ft.ptr.pp_double[j][0], 1, &fmatrix->ptr.pp_double[j][0], 1, ae_v_len(0,m-1), v);
        b.ptr.p_double[j] = w->ptr.p_double[j]*y->ptr.p_double[j];
    }
    rmatrixqr(&ft, n, m, &tau, _state);
    rmatrixqrunpackq(&ft, n, m, &tau, m, &q, _state);
    rmatrixqrunpackr(&ft, n, m, &r, _state);
    ae_vector_set_length(&tmp, m, _state);
    for(i=0; i<=m-1; i++)
    {
        tmp.ptr.p_double[i] = 0;
    }
    for(i=0; i<=n-1; i++)
    {
        v = b.ptr.p_double[i];
        ae_v_addd(&tmp.ptr.p_double[0], 1, &q.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
    }
    ae_vector_set_length(&b, m, _state);
    ae_v_move(&b.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * R contains reduced MxM design upper triangular matrix,
     * B contains reduced Mx1 right part.
     *
     * Determine system condition number and decide
     * should we use triangular solver (faster) or
     * SVD-based solver (more stable).
     *
     * We can use LU-based RCond estimator for this task.
     */
    rep->taskrcond = rmatrixlurcondinf(&r, m, _state);
    if( ae_fp_greater(rep->taskrcond,threshold) )
    {
        
        /*
         * use QR-based solver
         */
        ae_vector_set_length(c, m, _state);
        c->ptr.p_double[m-1] = b.ptr.p_double[m-1]/r.ptr.pp_double[m-1][m-1];
        for(i=m-2; i>=0; i--)
        {
            v = ae_v_dotproduct(&r.ptr.pp_double[i][i+1], 1, &c->ptr.p_double[i+1], 1, ae_v_len(i+1,m-1));
            c->ptr.p_double[i] = (b.ptr.p_double[i]-v)/r.ptr.pp_double[i][i];
        }
    }
    else
    {
        
        /*
         * use SVD-based solver
         */
        if( !rmatrixsvd(&r, m, m, 1, 1, 2, &sv, &u, &vt, _state) )
        {
            *info = -4;
            ae_frame_leave(_state);
            return;
        }
        ae_vector_set_length(&utb, m, _state);
        ae_vector_set_length(&sutb, m, _state);
        for(i=0; i<=m-1; i++)
        {
            utb.ptr.p_double[i] = 0;
        }
        for(i=0; i<=m-1; i++)
        {
            v = b.ptr.p_double[i];
            ae_v_addd(&utb.ptr.p_double[0], 1, &u.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
        }
        if( ae_fp_greater(sv.ptr.p_double[0],0) )
        {
            rep->taskrcond = sv.ptr.p_double[m-1]/sv.ptr.p_double[0];
            for(i=0; i<=m-1; i++)
            {
                if( ae_fp_greater(sv.ptr.p_double[i],threshold*sv.ptr.p_double[0]) )
                {
                    sutb.ptr.p_double[i] = utb.ptr.p_double[i]/sv.ptr.p_double[i];
                }
                else
                {
                    sutb.ptr.p_double[i] = 0;
                }
            }
        }
        else
        {
            rep->taskrcond = 0;
            for(i=0; i<=m-1; i++)
            {
                sutb.ptr.p_double[i] = 0;
            }
        }
        ae_vector_set_length(c, m, _state);
        for(i=0; i<=m-1; i++)
        {
            c->ptr.p_double[i] = 0;
        }
        for(i=0; i<=m-1; i++)
        {
            v = sutb.ptr.p_double[i];
            ae_v_addd(&c->ptr.p_double[0], 1, &vt.ptr.pp_double[i][0], 1, ae_v_len(0,m-1), v);
        }
    }
    
    /*
     * calculate errors
     */
    rep->rmserror = 0;
    rep->avgerror = 0;
    rep->avgrelerror = 0;
    rep->maxerror = 0;
    relcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&fmatrix->ptr.pp_double[i][0], 1, &c->ptr.p_double[0], 1, ae_v_len(0,m-1));
        rep->rmserror = rep->rmserror+ae_sqr(v-y->ptr.p_double[i], _state);
        rep->avgerror = rep->avgerror+ae_fabs(v-y->ptr.p_double[i], _state);
        if( ae_fp_neq(y->ptr.p_double[i],0) )
        {
            rep->avgrelerror = rep->avgrelerror+ae_fabs(v-y->ptr.p_double[i], _state)/ae_fabs(y->ptr.p_double[i], _state);
            relcnt = relcnt+1;
        }
        rep->maxerror = ae_maxreal(rep->maxerror, ae_fabs(v-y->ptr.p_double[i], _state), _state);
    }
    rep->rmserror = ae_sqrt(rep->rmserror/n, _state);
    rep->avgerror = rep->avgerror/n;
    if( relcnt!=0 )
    {
        rep->avgrelerror = rep->avgrelerror/relcnt;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine
*************************************************************************/
static void lsfit_lsfitclearrequestfields(lsfitstate* state,
     ae_state *_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->needfgh = ae_false;
    state->xupdated = ae_false;
}


/*************************************************************************
Internal subroutine, calculates barycentric basis functions.
Used for efficient simultaneous calculation of N basis functions.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
static void lsfit_barycentriccalcbasis(barycentricinterpolant* b,
     double t,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    double s2;
    double s;
    double v;
    ae_int_t i;
    ae_int_t j;


    
    /*
     * special case: N=1
     */
    if( b->n==1 )
    {
        y->ptr.p_double[0] = 1;
        return;
    }
    
    /*
     * Here we assume that task is normalized, i.e.:
     * 1. abs(Y[i])<=1
     * 2. abs(W[i])<=1
     * 3. X[] is ordered
     *
     * First, we decide: should we use "safe" formula (guarded
     * against overflow) or fast one?
     */
    s = ae_fabs(t-b->x.ptr.p_double[0], _state);
    for(i=0; i<=b->n-1; i++)
    {
        v = b->x.ptr.p_double[i];
        if( ae_fp_eq(v,t) )
        {
            for(j=0; j<=b->n-1; j++)
            {
                y->ptr.p_double[j] = 0;
            }
            y->ptr.p_double[i] = 1;
            return;
        }
        v = ae_fabs(t-v, _state);
        if( ae_fp_less(v,s) )
        {
            s = v;
        }
    }
    s2 = 0;
    for(i=0; i<=b->n-1; i++)
    {
        v = s/(t-b->x.ptr.p_double[i]);
        v = v*b->w.ptr.p_double[i];
        y->ptr.p_double[i] = v;
        s2 = s2+v;
    }
    v = 1/s2;
    ae_v_muld(&y->ptr.p_double[0], 1, ae_v_len(0,b->n-1), v);
}


/*************************************************************************
This is internal function for Chebyshev fitting.

It assumes that input data are normalized:
* X/XC belong to [-1,+1],
* mean(Y)=0, stddev(Y)=1.

It does not checks inputs for errors.

This function is used to fit general (shifted) Chebyshev models, power
basis models or barycentric models.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
    N   -   number of points, N>0.
    XC  -   points where polynomial values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that P(XC[i])=YC[i]
            * DC[i]=1   means that P'(XC[i])=YC[i]
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
    C   -   interpolant in Chebyshev form; [-1,+1] is used as base interval
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************/
static void lsfit_internalchebyshevfit(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t* info,
     /* Real    */ ae_vector* c,
     lsfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _xc;
    ae_vector _yc;
    ae_vector y2;
    ae_vector w2;
    ae_vector tmp;
    ae_vector tmp2;
    ae_vector tmpdiff;
    ae_vector bx;
    ae_vector by;
    ae_vector bw;
    ae_matrix fmatrix;
    ae_matrix cmatrix;
    ae_int_t i;
    ae_int_t j;
    double mx;
    double decay;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_xc, xc, _state, ae_true);
    xc = &_xc;
    ae_vector_init_copy(&_yc, yc, _state, ae_true);
    yc = &_yc;
    *info = 0;
    ae_vector_clear(c);
    _lsfitreport_clear(rep);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpdiff, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&by, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bw, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&fmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cmatrix, 0, 0, DT_REAL, _state, ae_true);

    
    /*
     * weight decay for correct handling of task which becomes
     * degenerate after constraints are applied
     */
    decay = 10000*ae_machineepsilon;
    
    /*
     * allocate space, initialize/fill:
     * * FMatrix-   values of basis functions at X[]
     * * CMatrix-   values (derivatives) of basis functions at XC[]
     * * fill constraints matrix
     * * fill first N rows of design matrix with values
     * * fill next M rows of design matrix with regularizing term
     * * append M zeros to Y
     * * append M elements, mean(abs(W)) each, to W
     */
    ae_vector_set_length(&y2, n+m, _state);
    ae_vector_set_length(&w2, n+m, _state);
    ae_vector_set_length(&tmp, m, _state);
    ae_vector_set_length(&tmpdiff, m, _state);
    ae_matrix_set_length(&fmatrix, n+m, m, _state);
    if( k>0 )
    {
        ae_matrix_set_length(&cmatrix, k, m+1, _state);
    }
    
    /*
     * Fill design matrix, Y2, W2:
     * * first N rows with basis functions for original points
     * * next M rows with decay terms
     */
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * prepare Ith row
         * use Tmp for calculations to avoid multidimensional arrays overhead
         */
        for(j=0; j<=m-1; j++)
        {
            if( j==0 )
            {
                tmp.ptr.p_double[j] = 1;
            }
            else
            {
                if( j==1 )
                {
                    tmp.ptr.p_double[j] = x->ptr.p_double[i];
                }
                else
                {
                    tmp.ptr.p_double[j] = 2*x->ptr.p_double[i]*tmp.ptr.p_double[j-1]-tmp.ptr.p_double[j-2];
                }
            }
        }
        ae_v_move(&fmatrix.ptr.pp_double[i][0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,m-1));
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=m-1; j++)
        {
            if( i==j )
            {
                fmatrix.ptr.pp_double[n+i][j] = decay;
            }
            else
            {
                fmatrix.ptr.pp_double[n+i][j] = 0;
            }
        }
    }
    ae_v_move(&y2.ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&w2.ptr.p_double[0], 1, &w->ptr.p_double[0], 1, ae_v_len(0,n-1));
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        mx = mx+ae_fabs(w->ptr.p_double[i], _state);
    }
    mx = mx/n;
    for(i=0; i<=m-1; i++)
    {
        y2.ptr.p_double[n+i] = 0;
        w2.ptr.p_double[n+i] = mx;
    }
    
    /*
     * fill constraints matrix
     */
    for(i=0; i<=k-1; i++)
    {
        
        /*
         * prepare Ith row
         * use Tmp for basis function values,
         * TmpDiff for basos function derivatives
         */
        for(j=0; j<=m-1; j++)
        {
            if( j==0 )
            {
                tmp.ptr.p_double[j] = 1;
                tmpdiff.ptr.p_double[j] = 0;
            }
            else
            {
                if( j==1 )
                {
                    tmp.ptr.p_double[j] = xc->ptr.p_double[i];
                    tmpdiff.ptr.p_double[j] = 1;
                }
                else
                {
                    tmp.ptr.p_double[j] = 2*xc->ptr.p_double[i]*tmp.ptr.p_double[j-1]-tmp.ptr.p_double[j-2];
                    tmpdiff.ptr.p_double[j] = 2*(tmp.ptr.p_double[j-1]+xc->ptr.p_double[i]*tmpdiff.ptr.p_double[j-1])-tmpdiff.ptr.p_double[j-2];
                }
            }
        }
        if( dc->ptr.p_int[i]==0 )
        {
            ae_v_move(&cmatrix.ptr.pp_double[i][0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,m-1));
        }
        if( dc->ptr.p_int[i]==1 )
        {
            ae_v_move(&cmatrix.ptr.pp_double[i][0], 1, &tmpdiff.ptr.p_double[0], 1, ae_v_len(0,m-1));
        }
        cmatrix.ptr.pp_double[i][m] = yc->ptr.p_double[i];
    }
    
    /*
     * Solve constrained task
     */
    if( k>0 )
    {
        
        /*
         * solve using regularization
         */
        lsfitlinearwc(&y2, &w2, &fmatrix, &cmatrix, n+m, m, k, info, c, rep, _state);
    }
    else
    {
        
        /*
         * no constraints, no regularization needed
         */
        lsfitlinearwc(y, w, &fmatrix, &cmatrix, n, m, 0, info, c, rep, _state);
    }
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal Floater-Hormann fitting subroutine for fixed D
*************************************************************************/
static void lsfit_barycentricfitwcfixedd(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* w,
     ae_int_t n,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* yc,
     /* Integer */ ae_vector* dc,
     ae_int_t k,
     ae_int_t m,
     ae_int_t d,
     ae_int_t* info,
     barycentricinterpolant* b,
     barycentricfitreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_matrix fmatrix;
    ae_matrix cmatrix;
    ae_vector y2;
    ae_vector w2;
    ae_vector sx;
    ae_vector sy;
    ae_vector sbf;
    ae_vector xoriginal;
    ae_vector yoriginal;
    ae_vector tmp;
    lsfitreport lrep;
    double v0;
    double v1;
    double mx;
    barycentricinterpolant b2;
    ae_int_t i;
    ae_int_t j;
    ae_int_t relcnt;
    double xa;
    double xb;
    double sa;
    double sb;
    double decay;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_vector_init_copy(&_w, w, _state, ae_true);
    w = &_w;
    ae_vector_init_copy(&_xc, xc, _state, ae_true);
    xc = &_xc;
    ae_vector_init_copy(&_yc, yc, _state, ae_true);
    yc = &_yc;
    *info = 0;
    _barycentricinterpolant_clear(b);
    _barycentricfitreport_clear(rep);
    ae_matrix_init(&fmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&cmatrix, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&w2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sbf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&yoriginal, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    _lsfitreport_init(&lrep, _state, ae_true);
    _barycentricinterpolant_init(&b2, _state, ae_true);

    if( ((n<1||m<2)||k<0)||k>=m )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=k-1; i++)
    {
        *info = 0;
        if( dc->ptr.p_int[i]<0 )
        {
            *info = -1;
        }
        if( dc->ptr.p_int[i]>1 )
        {
            *info = -1;
        }
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
    }
    
    /*
     * weight decay for correct handling of task which becomes
     * degenerate after constraints are applied
     */
    decay = 10000*ae_machineepsilon;
    
    /*
     * Scale X, Y, XC, YC
     */
    lsfitscalexy(x, y, w, n, xc, yc, dc, k, &xa, &xb, &sa, &sb, &xoriginal, &yoriginal, _state);
    
    /*
     * allocate space, initialize:
     * * FMatrix-   values of basis functions at X[]
     * * CMatrix-   values (derivatives) of basis functions at XC[]
     */
    ae_vector_set_length(&y2, n+m, _state);
    ae_vector_set_length(&w2, n+m, _state);
    ae_matrix_set_length(&fmatrix, n+m, m, _state);
    if( k>0 )
    {
        ae_matrix_set_length(&cmatrix, k, m+1, _state);
    }
    ae_vector_set_length(&y2, n+m, _state);
    ae_vector_set_length(&w2, n+m, _state);
    
    /*
     * Prepare design and constraints matrices:
     * * fill constraints matrix
     * * fill first N rows of design matrix with values
     * * fill next M rows of design matrix with regularizing term
     * * append M zeros to Y
     * * append M elements, mean(abs(W)) each, to W
     */
    ae_vector_set_length(&sx, m, _state);
    ae_vector_set_length(&sy, m, _state);
    ae_vector_set_length(&sbf, m, _state);
    for(j=0; j<=m-1; j++)
    {
        sx.ptr.p_double[j] = (double)(2*j)/(double)(m-1)-1;
    }
    for(i=0; i<=m-1; i++)
    {
        sy.ptr.p_double[i] = 1;
    }
    barycentricbuildfloaterhormann(&sx, &sy, m, d, &b2, _state);
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        lsfit_barycentriccalcbasis(&b2, x->ptr.p_double[i], &sbf, _state);
        ae_v_move(&fmatrix.ptr.pp_double[i][0], 1, &sbf.ptr.p_double[0], 1, ae_v_len(0,m-1));
        y2.ptr.p_double[i] = y->ptr.p_double[i];
        w2.ptr.p_double[i] = w->ptr.p_double[i];
        mx = mx+ae_fabs(w->ptr.p_double[i], _state)/n;
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=m-1; j++)
        {
            if( i==j )
            {
                fmatrix.ptr.pp_double[n+i][j] = decay;
            }
            else
            {
                fmatrix.ptr.pp_double[n+i][j] = 0;
            }
        }
        y2.ptr.p_double[n+i] = 0;
        w2.ptr.p_double[n+i] = mx;
    }
    if( k>0 )
    {
        for(j=0; j<=m-1; j++)
        {
            for(i=0; i<=m-1; i++)
            {
                sy.ptr.p_double[i] = 0;
            }
            sy.ptr.p_double[j] = 1;
            barycentricbuildfloaterhormann(&sx, &sy, m, d, &b2, _state);
            for(i=0; i<=k-1; i++)
            {
                ae_assert(dc->ptr.p_int[i]>=0&&dc->ptr.p_int[i]<=1, "BarycentricFit: internal error!", _state);
                barycentricdiff1(&b2, xc->ptr.p_double[i], &v0, &v1, _state);
                if( dc->ptr.p_int[i]==0 )
                {
                    cmatrix.ptr.pp_double[i][j] = v0;
                }
                if( dc->ptr.p_int[i]==1 )
                {
                    cmatrix.ptr.pp_double[i][j] = v1;
                }
            }
        }
        for(i=0; i<=k-1; i++)
        {
            cmatrix.ptr.pp_double[i][m] = yc->ptr.p_double[i];
        }
    }
    
    /*
     * Solve constrained task
     */
    if( k>0 )
    {
        
        /*
         * solve using regularization
         */
        lsfitlinearwc(&y2, &w2, &fmatrix, &cmatrix, n+m, m, k, info, &tmp, &lrep, _state);
    }
    else
    {
        
        /*
         * no constraints, no regularization needed
         */
        lsfitlinearwc(y, w, &fmatrix, &cmatrix, n, m, k, info, &tmp, &lrep, _state);
    }
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Generate interpolant and scale it
     */
    ae_v_move(&sy.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,m-1));
    barycentricbuildfloaterhormann(&sx, &sy, m, d, b, _state);
    barycentriclintransx(b, 2/(xb-xa), -(xa+xb)/(xb-xa), _state);
    barycentriclintransy(b, sb-sa, sa, _state);
    
    /*
     * Scale absolute errors obtained from LSFitLinearW.
     * Relative error should be calculated separately
     * (because of shifting/scaling of the task)
     */
    rep->taskrcond = lrep.taskrcond;
    rep->rmserror = lrep.rmserror*(sb-sa);
    rep->avgerror = lrep.avgerror*(sb-sa);
    rep->maxerror = lrep.maxerror*(sb-sa);
    rep->avgrelerror = 0;
    relcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(yoriginal.ptr.p_double[i],0) )
        {
            rep->avgrelerror = rep->avgrelerror+ae_fabs(barycentriccalc(b, xoriginal.ptr.p_double[i], _state)-yoriginal.ptr.p_double[i], _state)/ae_fabs(yoriginal.ptr.p_double[i], _state);
            relcnt = relcnt+1;
        }
    }
    if( relcnt!=0 )
    {
        rep->avgrelerror = rep->avgrelerror/relcnt;
    }
    ae_frame_leave(_state);
}


ae_bool _polynomialfitreport_init(polynomialfitreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _polynomialfitreport_init_copy(polynomialfitreport* dst, polynomialfitreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
    return ae_true;
}


void _polynomialfitreport_clear(polynomialfitreport* p)
{
}


ae_bool _barycentricfitreport_init(barycentricfitreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _barycentricfitreport_init_copy(barycentricfitreport* dst, barycentricfitreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->dbest = src->dbest;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
    return ae_true;
}


void _barycentricfitreport_clear(barycentricfitreport* p)
{
}


ae_bool _spline1dfitreport_init(spline1dfitreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _spline1dfitreport_init_copy(spline1dfitreport* dst, spline1dfitreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
    return ae_true;
}


void _spline1dfitreport_clear(spline1dfitreport* p)
{
}


ae_bool _lsfitreport_init(lsfitreport* p, ae_state *_state, ae_bool make_automatic)
{
    return ae_true;
}


ae_bool _lsfitreport_init_copy(lsfitreport* dst, lsfitreport* src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
    return ae_true;
}


void _lsfitreport_clear(lsfitreport* p)
{
}


ae_bool _lsfitstate_init(lsfitstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_matrix_init(&p->taskx, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tasky, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->w, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->c, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->h, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_minlmstate_init(&p->optstate, _state, make_automatic) )
        return ae_false;
    if( !_minlmreport_init(&p->optrep, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _lsfitstate_init_copy(lsfitstate* dst, lsfitstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->m = src->m;
    dst->k = src->k;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->stpmax = src->stpmax;
    dst->xrep = src->xrep;
    if( !ae_matrix_init_copy(&dst->taskx, &src->taskx, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tasky, &src->tasky, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->w, &src->w, _state, make_automatic) )
        return ae_false;
    dst->xupdated = src->xupdated;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->needfgh = src->needfgh;
    dst->pointindex = src->pointindex;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->c, &src->c, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->h, &src->h, _state, make_automatic) )
        return ae_false;
    dst->repterminationtype = src->repterminationtype;
    dst->reprmserror = src->reprmserror;
    dst->repavgerror = src->repavgerror;
    dst->repavgrelerror = src->repavgrelerror;
    dst->repmaxerror = src->repmaxerror;
    if( !_minlmstate_init_copy(&dst->optstate, &src->optstate, _state, make_automatic) )
        return ae_false;
    if( !_minlmreport_init_copy(&dst->optrep, &src->optrep, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _lsfitstate_clear(lsfitstate* p)
{
    ae_matrix_clear(&p->taskx);
    ae_vector_clear(&p->tasky);
    ae_vector_clear(&p->w);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->c);
    ae_vector_clear(&p->g);
    ae_matrix_clear(&p->h);
    _minlmstate_clear(&p->optstate);
    _minlmreport_clear(&p->optrep);
    _rcommstate_clear(&p->rstate);
}


/*$ End $*/
