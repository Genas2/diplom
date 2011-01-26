/*************************************************************************
Copyright (c) 2006-2010, Sergey Bochkanov (ALGLIB project).

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
#include "pspline.h"


/*$ Declarations $*/
static void pspline_pspline2par(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t pt,
     /* Real    */ ae_vector* p,
     ae_state *_state);
static void pspline_pspline3par(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t pt,
     /* Real    */ ae_vector* p,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
This function  builds  non-periodic 2-dimensional parametric spline  which
starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            Order of points is important!
    N   -   points count, N>=5 for Akima splines, N>=2 for other types  of
            splines.
    ST  -   spline type:
            * 0     Akima spline
            * 1     parabolically terminated Catmull-Rom spline (Tension=0)
            * 2     parabolically terminated cubic spline
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that  there all consequent points  are distinct.
  I.e. (x0,y0)<>(x1,y1),  (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2build(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t st,
     ae_int_t pt,
     pspline2interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _xy;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_xy, xy, _state, ae_true);
    xy = &_xy;
    _pspline2interpolant_clear(p);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_assert(st>=0&&st<=2, "PSpline2Build: incorrect spline type!", _state);
    ae_assert(pt>=0&&pt<=2, "PSpline2Build: incorrect parameterization type!", _state);
    if( st==0 )
    {
        ae_assert(n>=5, "PSpline2Build: N<5 (minimum value for Akima splines)!", _state);
    }
    else
    {
        ae_assert(n>=2, "PSpline2Build: N<2!", _state);
    }
    
    /*
     * Prepare
     */
    p->n = n;
    p->periodic = ae_false;
    ae_vector_set_length(&tmp, n, _state);
    
    /*
     * Build parameterization, check that all parameters are distinct
     */
    pspline_pspline2par(xy, n, pt, &p->p, _state);
    ae_assert(aredistinct(&p->p, n, _state), "PSpline2Build: consequent points are too close!", _state);
    
    /*
     * Build splines
     */
    if( st==0 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildakima(&p->p, &tmp, n, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildakima(&p->p, &tmp, n, &p->y, _state);
    }
    if( st==1 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcatmullrom(&p->p, &tmp, n, 0, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcatmullrom(&p->p, &tmp, n, 0, 0.0, &p->y, _state);
    }
    if( st==2 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcubic(&p->p, &tmp, n, 0, 0.0, 0, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcubic(&p->p, &tmp, n, 0, 0.0, 0, 0.0, &p->y, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function  builds  non-periodic 3-dimensional parametric spline  which
starts at (X[0],Y[0],Z[0]) and ends at (X[N-1],Y[N-1],Z[N-1]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3build(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t st,
     ae_int_t pt,
     pspline3interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _xy;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_xy, xy, _state, ae_true);
    xy = &_xy;
    _pspline3interpolant_clear(p);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_assert(st>=0&&st<=2, "PSpline3Build: incorrect spline type!", _state);
    ae_assert(pt>=0&&pt<=2, "PSpline3Build: incorrect parameterization type!", _state);
    if( st==0 )
    {
        ae_assert(n>=5, "PSpline3Build: N<5 (minimum value for Akima splines)!", _state);
    }
    else
    {
        ae_assert(n>=2, "PSpline3Build: N<2!", _state);
    }
    
    /*
     * Prepare
     */
    p->n = n;
    p->periodic = ae_false;
    ae_vector_set_length(&tmp, n, _state);
    
    /*
     * Build parameterization, check that all parameters are distinct
     */
    pspline_pspline3par(xy, n, pt, &p->p, _state);
    ae_assert(aredistinct(&p->p, n, _state), "PSpline3Build: consequent points are too close!", _state);
    
    /*
     * Build splines
     */
    if( st==0 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildakima(&p->p, &tmp, n, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildakima(&p->p, &tmp, n, &p->y, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][2], xy->stride, ae_v_len(0,n-1));
        spline1dbuildakima(&p->p, &tmp, n, &p->z, _state);
    }
    if( st==1 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcatmullrom(&p->p, &tmp, n, 0, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcatmullrom(&p->p, &tmp, n, 0, 0.0, &p->y, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][2], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcatmullrom(&p->p, &tmp, n, 0, 0.0, &p->z, _state);
    }
    if( st==2 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcubic(&p->p, &tmp, n, 0, 0.0, 0, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcubic(&p->p, &tmp, n, 0, 0.0, 0, 0.0, &p->y, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][2], xy->stride, ae_v_len(0,n-1));
        spline1dbuildcubic(&p->p, &tmp, n, 0, 0.0, 0, 0.0, &p->z, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This  function  builds  periodic  2-dimensional  parametric  spline  which
starts at (X[0],Y[0]), goes through all points to (X[N-1],Y[N-1]) and then
back to (X[0],Y[0]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            XY[N-1,0:1] must be different from XY[0,0:1].
            Order of points is important!
    N   -   points count, N>=3 for other types of splines.
    ST  -   spline type:
            * 1     Catmull-Rom spline (Tension=0) with cyclic boundary conditions
            * 2     cubic spline with cyclic boundary conditions
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that there all consequent points  are  distinct.
  I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so  on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).
* last point of sequence is NOT equal to the first  point.  You  shouldn't
  make curve "explicitly periodic" by making them equal.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2buildperiodic(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t st,
     ae_int_t pt,
     pspline2interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _xy;
    ae_matrix xyp;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_xy, xy, _state, ae_true);
    xy = &_xy;
    _pspline2interpolant_clear(p);
    ae_matrix_init(&xyp, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_assert(st>=1&&st<=2, "PSpline2BuildPeriodic: incorrect spline type!", _state);
    ae_assert(pt>=0&&pt<=2, "PSpline2BuildPeriodic: incorrect parameterization type!", _state);
    ae_assert(n>=3, "PSpline2BuildPeriodic: N<3!", _state);
    
    /*
     * Prepare
     */
    p->n = n;
    p->periodic = ae_true;
    ae_vector_set_length(&tmp, n+1, _state);
    ae_matrix_set_length(&xyp, n+1, 2, _state);
    ae_v_move(&xyp.ptr.pp_double[0][0], xyp.stride, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
    ae_v_move(&xyp.ptr.pp_double[0][1], xyp.stride, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
    ae_v_move(&xyp.ptr.pp_double[n][0], 1, &xy->ptr.pp_double[0][0], 1, ae_v_len(0,1));
    
    /*
     * Build parameterization, check that all parameters are distinct
     */
    pspline_pspline2par(&xyp, n+1, pt, &p->p, _state);
    ae_assert(aredistinct(&p->p, n+1, _state), "PSpline2BuildPeriodic: consequent (or first and last) points are too close!", _state);
    
    /*
     * Build splines
     */
    if( st==1 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][0], xyp.stride, ae_v_len(0,n));
        spline1dbuildcatmullrom(&p->p, &tmp, n+1, -1, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][1], xyp.stride, ae_v_len(0,n));
        spline1dbuildcatmullrom(&p->p, &tmp, n+1, -1, 0.0, &p->y, _state);
    }
    if( st==2 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][0], xyp.stride, ae_v_len(0,n));
        spline1dbuildcubic(&p->p, &tmp, n+1, -1, 0.0, -1, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][1], xyp.stride, ae_v_len(0,n));
        spline1dbuildcubic(&p->p, &tmp, n+1, -1, 0.0, -1, 0.0, &p->y, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This  function  builds  periodic  3-dimensional  parametric  spline  which
starts at (X[0],Y[0],Z[0]), goes through all points to (X[N-1],Y[N-1],Z[N-1])
and then back to (X[0],Y[0],Z[0]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3buildperiodic(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t st,
     ae_int_t pt,
     pspline3interpolant* p,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _xy;
    ae_matrix xyp;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_xy, xy, _state, ae_true);
    xy = &_xy;
    _pspline3interpolant_clear(p);
    ae_matrix_init(&xyp, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_assert(st>=1&&st<=2, "PSpline3BuildPeriodic: incorrect spline type!", _state);
    ae_assert(pt>=0&&pt<=2, "PSpline3BuildPeriodic: incorrect parameterization type!", _state);
    ae_assert(n>=3, "PSpline3BuildPeriodic: N<3!", _state);
    
    /*
     * Prepare
     */
    p->n = n;
    p->periodic = ae_true;
    ae_vector_set_length(&tmp, n+1, _state);
    ae_matrix_set_length(&xyp, n+1, 3, _state);
    ae_v_move(&xyp.ptr.pp_double[0][0], xyp.stride, &xy->ptr.pp_double[0][0], xy->stride, ae_v_len(0,n-1));
    ae_v_move(&xyp.ptr.pp_double[0][1], xyp.stride, &xy->ptr.pp_double[0][1], xy->stride, ae_v_len(0,n-1));
    ae_v_move(&xyp.ptr.pp_double[0][2], xyp.stride, &xy->ptr.pp_double[0][2], xy->stride, ae_v_len(0,n-1));
    ae_v_move(&xyp.ptr.pp_double[n][0], 1, &xy->ptr.pp_double[0][0], 1, ae_v_len(0,2));
    
    /*
     * Build parameterization, check that all parameters are distinct
     */
    pspline_pspline3par(&xyp, n+1, pt, &p->p, _state);
    ae_assert(aredistinct(&p->p, n+1, _state), "PSplineBuild2Periodic: consequent (or first and last) points are too close!", _state);
    
    /*
     * Build splines
     */
    if( st==1 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][0], xyp.stride, ae_v_len(0,n));
        spline1dbuildcatmullrom(&p->p, &tmp, n+1, -1, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][1], xyp.stride, ae_v_len(0,n));
        spline1dbuildcatmullrom(&p->p, &tmp, n+1, -1, 0.0, &p->y, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][2], xyp.stride, ae_v_len(0,n));
        spline1dbuildcatmullrom(&p->p, &tmp, n+1, -1, 0.0, &p->z, _state);
    }
    if( st==2 )
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][0], xyp.stride, ae_v_len(0,n));
        spline1dbuildcubic(&p->p, &tmp, n+1, -1, 0.0, -1, 0.0, &p->x, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][1], xyp.stride, ae_v_len(0,n));
        spline1dbuildcubic(&p->p, &tmp, n+1, -1, 0.0, -1, 0.0, &p->y, _state);
        ae_v_move(&tmp.ptr.p_double[0], 1, &xyp.ptr.pp_double[0][2], xyp.stride, ae_v_len(0,n));
        spline1dbuildcubic(&p->p, &tmp, n+1, -1, 0.0, -1, 0.0, &p->z, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function returns vector of parameter values correspoding to points.

I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and U=TValues(P)  we
have
    (X[0],Y[0]) = PSpline2Calc(P,U[0]),
    (X[1],Y[1]) = PSpline2Calc(P,U[1]),
    (X[2],Y[2]) = PSpline2Calc(P,U[2]),
    ...

INPUT PARAMETERS:
    P   -   parametric spline interpolant

OUTPUT PARAMETERS:
    N   -   array size
    T   -   array[0..N-1]


NOTES:
* for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1
* for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2parametervalues(pspline2interpolant* p,
     ae_int_t* n,
     /* Real    */ ae_vector* t,
     ae_state *_state)
{

    *n = 0;
    ae_vector_clear(t);

    ae_assert(p->n>=2, "PSpline2ParameterValues: internal error!", _state);
    *n = p->n;
    ae_vector_set_length(t, *n, _state);
    ae_v_move(&t->ptr.p_double[0], 1, &p->p.ptr.p_double[0], 1, ae_v_len(0,*n-1));
    t->ptr.p_double[0] = 0;
    if( !p->periodic )
    {
        t->ptr.p_double[*n-1] = 1;
    }
}


/*************************************************************************
This function returns vector of parameter values correspoding to points.

Same as PSpline2ParameterValues(), but for 3D.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3parametervalues(pspline3interpolant* p,
     ae_int_t* n,
     /* Real    */ ae_vector* t,
     ae_state *_state)
{

    *n = 0;
    ae_vector_clear(t);

    ae_assert(p->n>=2, "PSpline3ParameterValues: internal error!", _state);
    *n = p->n;
    ae_vector_set_length(t, *n, _state);
    ae_v_move(&t->ptr.p_double[0], 1, &p->p.ptr.p_double[0], 1, ae_v_len(0,*n-1));
    t->ptr.p_double[0] = 0;
    if( !p->periodic )
    {
        t->ptr.p_double[*n-1] = 1;
    }
}


/*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2calc(pspline2interpolant* p,
     double t,
     double* x,
     double* y,
     ae_state *_state)
{

    *x = 0;
    *y = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    *x = spline1dcalc(&p->x, t, _state);
    *y = spline1dcalc(&p->y, t, _state);
}


/*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position
    Z   -   Z-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3calc(pspline3interpolant* p,
     double t,
     double* x,
     double* y,
     double* z,
     ae_state *_state)
{

    *x = 0;
    *y = 0;
    *z = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    *x = spline1dcalc(&p->x, t, _state);
    *y = spline1dcalc(&p->y, t, _state);
    *z = spline1dcalc(&p->z, t, _state);
}


/*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    
NOTE:
    X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2tangent(pspline2interpolant* p,
     double t,
     double* x,
     double* y,
     ae_state *_state)
{
    double v;
    double v0;
    double v1;

    *x = 0;
    *y = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    pspline2diff(p, t, &v0, x, &v1, y, _state);
    if( ae_fp_neq(*x,0)||ae_fp_neq(*y,0) )
    {
        
        /*
         * this code is a bit more complex than X^2+Y^2 to avoid
         * overflow for large values of X and Y.
         */
        v = safepythag2(*x, *y, _state);
        *x = *x/v;
        *y = *y/v;
    }
}


/*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    Z    -   Z-component of tangent vector (normalized)

NOTE:
    X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3tangent(pspline3interpolant* p,
     double t,
     double* x,
     double* y,
     double* z,
     ae_state *_state)
{
    double v;
    double v0;
    double v1;
    double v2;

    *x = 0;
    *y = 0;
    *z = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    pspline3diff(p, t, &v0, x, &v1, y, &v2, z, _state);
    if( (ae_fp_neq(*x,0)||ae_fp_neq(*y,0))||ae_fp_neq(*z,0) )
    {
        v = safepythag3(*x, *y, *z, _state);
        *x = *x/v;
        *y = *y/v;
        *z = *z/v;
    }
}


/*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2diff(pspline2interpolant* p,
     double t,
     double* x,
     double* dx,
     double* y,
     double* dy,
     ae_state *_state)
{
    double d2s;

    *x = 0;
    *dx = 0;
    *y = 0;
    *dy = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    spline1ddiff(&p->x, t, x, dx, &d2s, _state);
    spline1ddiff(&p->y, t, y, dy, &d2s, _state);
}


/*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT,dZ/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative
    Z   -   Z-value
    DZ  -   Z-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3diff(pspline3interpolant* p,
     double t,
     double* x,
     double* dx,
     double* y,
     double* dy,
     double* z,
     double* dz,
     ae_state *_state)
{
    double d2s;

    *x = 0;
    *dx = 0;
    *y = 0;
    *dy = 0;
    *z = 0;
    *dz = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    spline1ddiff(&p->x, t, x, dx, &d2s, _state);
    spline1ddiff(&p->y, t, y, dy, &d2s, _state);
    spline1ddiff(&p->z, t, z, dz, &d2s, _state);
}


/*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2diff2(pspline2interpolant* p,
     double t,
     double* x,
     double* dx,
     double* d2x,
     double* y,
     double* dy,
     double* d2y,
     ae_state *_state)
{

    *x = 0;
    *dx = 0;
    *d2x = 0;
    *y = 0;
    *dy = 0;
    *d2y = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    spline1ddiff(&p->x, t, x, dx, d2x, _state);
    spline1ddiff(&p->y, t, y, dy, d2y, _state);
}


/*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative
    Z   -   Z-value
    DZ  -   derivative
    D2Z -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3diff2(pspline3interpolant* p,
     double t,
     double* x,
     double* dx,
     double* d2x,
     double* y,
     double* dy,
     double* d2y,
     double* z,
     double* dz,
     double* d2z,
     ae_state *_state)
{

    *x = 0;
    *dx = 0;
    *d2x = 0;
    *y = 0;
    *dy = 0;
    *d2y = 0;
    *z = 0;
    *dz = 0;
    *d2z = 0;

    if( p->periodic )
    {
        t = t-ae_ifloor(t, _state);
    }
    spline1ddiff(&p->x, t, x, dx, d2x, _state);
    spline1ddiff(&p->y, t, y, dy, d2y, _state);
    spline1ddiff(&p->z, t, z, dz, d2z, _state);
}


/*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************/
double pspline2arclength(pspline2interpolant* p,
     double a,
     double b,
     ae_state *_state)
{
    ae_frame _frame_block;
    autogkstate state;
    autogkreport rep;
    double sx;
    double dsx;
    double d2sx;
    double sy;
    double dsy;
    double d2sy;
    double result;

    ae_frame_make(_state, &_frame_block);
    _autogkstate_init(&state, _state, ae_true);
    _autogkreport_init(&rep, _state, ae_true);

    autogksmooth(a, b, &state, _state);
    while(autogkiteration(&state, _state))
    {
        spline1ddiff(&p->x, state.x, &sx, &dsx, &d2sx, _state);
        spline1ddiff(&p->y, state.x, &sy, &dsy, &d2sy, _state);
        state.f = safepythag2(dsx, dsy, _state);
    }
    autogkresults(&state, &result, &rep, _state);
    ae_assert(rep.terminationtype>0, "PSpline2ArcLength: internal error!", _state);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************/
double pspline3arclength(pspline3interpolant* p,
     double a,
     double b,
     ae_state *_state)
{
    ae_frame _frame_block;
    autogkstate state;
    autogkreport rep;
    double sx;
    double dsx;
    double d2sx;
    double sy;
    double dsy;
    double d2sy;
    double sz;
    double dsz;
    double d2sz;
    double result;

    ae_frame_make(_state, &_frame_block);
    _autogkstate_init(&state, _state, ae_true);
    _autogkreport_init(&rep, _state, ae_true);

    autogksmooth(a, b, &state, _state);
    while(autogkiteration(&state, _state))
    {
        spline1ddiff(&p->x, state.x, &sx, &dsx, &d2sx, _state);
        spline1ddiff(&p->y, state.x, &sy, &dsy, &d2sy, _state);
        spline1ddiff(&p->z, state.x, &sz, &dsz, &d2sz, _state);
        state.f = safepythag3(dsx, dsy, dsz, _state);
    }
    autogkresults(&state, &result, &rep, _state);
    ae_assert(rep.terminationtype>0, "PSpline3ArcLength: internal error!", _state);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Builds non-periodic parameterization for 2-dimensional spline
*************************************************************************/
static void pspline_pspline2par(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t pt,
     /* Real    */ ae_vector* p,
     ae_state *_state)
{
    double v;
    ae_int_t i;

    ae_vector_clear(p);

    ae_assert(pt>=0&&pt<=2, "PSpline2Par: internal error!", _state);
    
    /*
     * Build parameterization:
     * * fill by non-normalized values
     * * normalize them so we have P[0]=0, P[N-1]=1.
     */
    ae_vector_set_length(p, n, _state);
    if( pt==0 )
    {
        for(i=0; i<=n-1; i++)
        {
            p->ptr.p_double[i] = i;
        }
    }
    if( pt==1 )
    {
        p->ptr.p_double[0] = 0;
        for(i=1; i<=n-1; i++)
        {
            p->ptr.p_double[i] = p->ptr.p_double[i-1]+safepythag2(xy->ptr.pp_double[i][0]-xy->ptr.pp_double[i-1][0], xy->ptr.pp_double[i][1]-xy->ptr.pp_double[i-1][1], _state);
        }
    }
    if( pt==2 )
    {
        p->ptr.p_double[0] = 0;
        for(i=1; i<=n-1; i++)
        {
            p->ptr.p_double[i] = p->ptr.p_double[i-1]+ae_sqrt(safepythag2(xy->ptr.pp_double[i][0]-xy->ptr.pp_double[i-1][0], xy->ptr.pp_double[i][1]-xy->ptr.pp_double[i-1][1], _state), _state);
        }
    }
    v = 1/p->ptr.p_double[n-1];
    ae_v_muld(&p->ptr.p_double[0], 1, ae_v_len(0,n-1), v);
}


/*************************************************************************
Builds non-periodic parameterization for 3-dimensional spline
*************************************************************************/
static void pspline_pspline3par(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t pt,
     /* Real    */ ae_vector* p,
     ae_state *_state)
{
    double v;
    ae_int_t i;

    ae_vector_clear(p);

    ae_assert(pt>=0&&pt<=2, "PSpline3Par: internal error!", _state);
    
    /*
     * Build parameterization:
     * * fill by non-normalized values
     * * normalize them so we have P[0]=0, P[N-1]=1.
     */
    ae_vector_set_length(p, n, _state);
    if( pt==0 )
    {
        for(i=0; i<=n-1; i++)
        {
            p->ptr.p_double[i] = i;
        }
    }
    if( pt==1 )
    {
        p->ptr.p_double[0] = 0;
        for(i=1; i<=n-1; i++)
        {
            p->ptr.p_double[i] = p->ptr.p_double[i-1]+safepythag3(xy->ptr.pp_double[i][0]-xy->ptr.pp_double[i-1][0], xy->ptr.pp_double[i][1]-xy->ptr.pp_double[i-1][1], xy->ptr.pp_double[i][2]-xy->ptr.pp_double[i-1][2], _state);
        }
    }
    if( pt==2 )
    {
        p->ptr.p_double[0] = 0;
        for(i=1; i<=n-1; i++)
        {
            p->ptr.p_double[i] = p->ptr.p_double[i-1]+ae_sqrt(safepythag3(xy->ptr.pp_double[i][0]-xy->ptr.pp_double[i-1][0], xy->ptr.pp_double[i][1]-xy->ptr.pp_double[i-1][1], xy->ptr.pp_double[i][2]-xy->ptr.pp_double[i-1][2], _state), _state);
        }
    }
    v = 1/p->ptr.p_double[n-1];
    ae_v_muld(&p->ptr.p_double[0], 1, ae_v_len(0,n-1), v);
}


ae_bool _pspline2interpolant_init(pspline2interpolant* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->p, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init(&p->x, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init(&p->y, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _pspline2interpolant_init_copy(pspline2interpolant* dst, pspline2interpolant* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->periodic = src->periodic;
    if( !ae_vector_init_copy(&dst->p, &src->p, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _pspline2interpolant_clear(pspline2interpolant* p)
{
    ae_vector_clear(&p->p);
    _spline1dinterpolant_clear(&p->x);
    _spline1dinterpolant_clear(&p->y);
}


ae_bool _pspline3interpolant_init(pspline3interpolant* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->p, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init(&p->x, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init(&p->y, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init(&p->z, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _pspline3interpolant_init_copy(pspline3interpolant* dst, pspline3interpolant* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->periodic = src->periodic;
    if( !ae_vector_init_copy(&dst->p, &src->p, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    if( !_spline1dinterpolant_init_copy(&dst->z, &src->z, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _pspline3interpolant_clear(pspline3interpolant* p)
{
    ae_vector_clear(&p->p);
    _spline1dinterpolant_clear(&p->x);
    _spline1dinterpolant_clear(&p->y);
    _spline1dinterpolant_clear(&p->z);
}


/*$ End $*/
