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
#include "spline2d.h"


/*$ Declarations $*/
static void spline2d_bicubiccalcderivatives(/* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* dx,
     /* Real    */ ae_matrix* dy,
     /* Real    */ ae_matrix* dxy,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
This subroutine builds bilinear spline coefficients table.

Input parameters:
    X   -   spline abscissas, array[0..N-1]
    Y   -   spline ordinates, array[0..M-1]
    F   -   function values, array[0..M-1,0..N-1]
    M,N -   grid size, M>=2, N>=2

Output parameters:
    C   -   spline interpolant

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dbuildbilinear(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* f,
     ae_int_t m,
     ae_int_t n,
     spline2dinterpolant* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_matrix _f;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tblsize;
    ae_int_t shift;
    double t;
    ae_matrix dx;
    ae_matrix dy;
    ae_matrix dxy;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_matrix_init_copy(&_f, f, _state, ae_true);
    f = &_f;
    _spline2dinterpolant_clear(c);
    ae_matrix_init(&dx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&dy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&dxy, 0, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=2&&m>=2, "Spline2DBuildBilinear: N<2 or M<2!", _state);
    
    /*
     * Sort points
     */
    for(j=0; j<=n-1; j++)
    {
        k = j;
        for(i=j+1; i<=n-1; i++)
        {
            if( ae_fp_less(x->ptr.p_double[i],x->ptr.p_double[k]) )
            {
                k = i;
            }
        }
        if( k!=j )
        {
            for(i=0; i<=m-1; i++)
            {
                t = f->ptr.pp_double[i][j];
                f->ptr.pp_double[i][j] = f->ptr.pp_double[i][k];
                f->ptr.pp_double[i][k] = t;
            }
            t = x->ptr.p_double[j];
            x->ptr.p_double[j] = x->ptr.p_double[k];
            x->ptr.p_double[k] = t;
        }
    }
    for(i=0; i<=m-1; i++)
    {
        k = i;
        for(j=i+1; j<=m-1; j++)
        {
            if( ae_fp_less(y->ptr.p_double[j],y->ptr.p_double[k]) )
            {
                k = j;
            }
        }
        if( k!=i )
        {
            for(j=0; j<=n-1; j++)
            {
                t = f->ptr.pp_double[i][j];
                f->ptr.pp_double[i][j] = f->ptr.pp_double[k][j];
                f->ptr.pp_double[k][j] = t;
            }
            t = y->ptr.p_double[i];
            y->ptr.p_double[i] = y->ptr.p_double[k];
            y->ptr.p_double[k] = t;
        }
    }
    
    /*
     * Fill C:
     *  C[0]            -   length(C)
     *  C[1]            -   type(C):
     *                      -1 = bilinear interpolant
     *                      -3 = general cubic spline
     *                           (see BuildBicubicSpline)
     *  C[2]:
     *      N (x count)
     *  C[3]:
     *      M (y count)
     *  C[4]...C[4+N-1]:
     *      x[i], i = 0...N-1
     *  C[4+N]...C[4+N+M-1]:
     *      y[i], i = 0...M-1
     *  C[4+N+M]...C[4+N+M+(N*M-1)]:
     *      f(i,j) table. f(0,0), f(0, 1), f(0,2) and so on...
     */
    c->k = 1;
    tblsize = 4+n+m+n*m;
    ae_vector_set_length(&c->c, tblsize-1+1, _state);
    c->c.ptr.p_double[0] = tblsize;
    c->c.ptr.p_double[1] = -1;
    c->c.ptr.p_double[2] = n;
    c->c.ptr.p_double[3] = m;
    for(i=0; i<=n-1; i++)
    {
        c->c.ptr.p_double[4+i] = x->ptr.p_double[i];
    }
    for(i=0; i<=m-1; i++)
    {
        c->c.ptr.p_double[4+n+i] = y->ptr.p_double[i];
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            shift = i*n+j;
            c->c.ptr.p_double[4+n+m+shift] = f->ptr.pp_double[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine builds bicubic spline coefficients table.

Input parameters:
    X   -   spline abscissas, array[0..N-1]
    Y   -   spline ordinates, array[0..M-1]
    F   -   function values, array[0..M-1,0..N-1]
    M,N -   grid size, M>=2, N>=2

Output parameters:
    C   -   spline interpolant

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dbuildbicubic(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_matrix* f,
     ae_int_t m,
     ae_int_t n,
     spline2dinterpolant* c,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_matrix _f;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tblsize;
    ae_int_t shift;
    double t;
    ae_matrix dx;
    ae_matrix dy;
    ae_matrix dxy;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state, ae_true);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state, ae_true);
    y = &_y;
    ae_matrix_init_copy(&_f, f, _state, ae_true);
    f = &_f;
    _spline2dinterpolant_clear(c);
    ae_matrix_init(&dx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&dy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&dxy, 0, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=2&&m>=2, "BuildBicubicSpline: N<2 or M<2!", _state);
    
    /*
     * Sort points
     */
    for(j=0; j<=n-1; j++)
    {
        k = j;
        for(i=j+1; i<=n-1; i++)
        {
            if( ae_fp_less(x->ptr.p_double[i],x->ptr.p_double[k]) )
            {
                k = i;
            }
        }
        if( k!=j )
        {
            for(i=0; i<=m-1; i++)
            {
                t = f->ptr.pp_double[i][j];
                f->ptr.pp_double[i][j] = f->ptr.pp_double[i][k];
                f->ptr.pp_double[i][k] = t;
            }
            t = x->ptr.p_double[j];
            x->ptr.p_double[j] = x->ptr.p_double[k];
            x->ptr.p_double[k] = t;
        }
    }
    for(i=0; i<=m-1; i++)
    {
        k = i;
        for(j=i+1; j<=m-1; j++)
        {
            if( ae_fp_less(y->ptr.p_double[j],y->ptr.p_double[k]) )
            {
                k = j;
            }
        }
        if( k!=i )
        {
            for(j=0; j<=n-1; j++)
            {
                t = f->ptr.pp_double[i][j];
                f->ptr.pp_double[i][j] = f->ptr.pp_double[k][j];
                f->ptr.pp_double[k][j] = t;
            }
            t = y->ptr.p_double[i];
            y->ptr.p_double[i] = y->ptr.p_double[k];
            y->ptr.p_double[k] = t;
        }
    }
    
    /*
     * Fill C:
     *  C[0]            -   length(C)
     *  C[1]            -   type(C):
     *                      -1 = bilinear interpolant
     *                           (see BuildBilinearInterpolant)
     *                      -3 = general cubic spline
     *  C[2]:
     *      N (x count)
     *  C[3]:
     *      M (y count)
     *  C[4]...C[4+N-1]:
     *      x[i], i = 0...N-1
     *  C[4+N]...C[4+N+M-1]:
     *      y[i], i = 0...M-1
     *  C[4+N+M]...C[4+N+M+(N*M-1)]:
     *      f(i,j) table. f(0,0), f(0, 1), f(0,2) and so on...
     *  C[4+N+M+N*M]...C[4+N+M+(2*N*M-1)]:
     *      df(i,j)/dx table.
     *  C[4+N+M+2*N*M]...C[4+N+M+(3*N*M-1)]:
     *      df(i,j)/dy table.
     *  C[4+N+M+3*N*M]...C[4+N+M+(4*N*M-1)]:
     *      d2f(i,j)/dxdy table.
     */
    c->k = 3;
    tblsize = 4+n+m+4*n*m;
    ae_vector_set_length(&c->c, tblsize-1+1, _state);
    c->c.ptr.p_double[0] = tblsize;
    c->c.ptr.p_double[1] = -3;
    c->c.ptr.p_double[2] = n;
    c->c.ptr.p_double[3] = m;
    for(i=0; i<=n-1; i++)
    {
        c->c.ptr.p_double[4+i] = x->ptr.p_double[i];
    }
    for(i=0; i<=m-1; i++)
    {
        c->c.ptr.p_double[4+n+i] = y->ptr.p_double[i];
    }
    spline2d_bicubiccalcderivatives(f, x, y, m, n, &dx, &dy, &dxy, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            shift = i*n+j;
            c->c.ptr.p_double[4+n+m+shift] = f->ptr.pp_double[i][j];
            c->c.ptr.p_double[4+n+m+n*m+shift] = dx.ptr.pp_double[i][j];
            c->c.ptr.p_double[4+n+m+2*n*m+shift] = dy.ptr.pp_double[i][j];
            c->c.ptr.p_double[4+n+m+3*n*m+shift] = dxy.ptr.pp_double[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine calculates the value of the bilinear or bicubic spline  at
the given point X.

Input parameters:
    C   -   coefficients table.
            Built by BuildBilinearSpline or BuildBicubicSpline.
    X, Y-   point

Result:
    S(x,y)

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************/
double spline2dcalc(spline2dinterpolant* c,
     double x,
     double y,
     ae_state *_state)
{
    double v;
    double vx;
    double vy;
    double vxy;
    double result;


    spline2ddiff(c, x, y, &v, &vx, &vy, &vxy, _state);
    result = v;
    return result;
}


/*************************************************************************
This subroutine calculates the value of the bilinear or bicubic spline  at
the given point X and its derivatives.

Input parameters:
    C   -   spline interpolant.
    X, Y-   point

Output parameters:
    F   -   S(x,y)
    FX  -   dS(x,y)/dX
    FY  -   dS(x,y)/dY
    FXY -   d2S(x,y)/dXdY

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************/
void spline2ddiff(spline2dinterpolant* c,
     double x,
     double y,
     double* f,
     double* fx,
     double* fy,
     double* fxy,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    double t;
    double dt;
    double u;
    double du;
    ae_int_t ix;
    ae_int_t iy;
    ae_int_t l;
    ae_int_t r;
    ae_int_t h;
    ae_int_t shift1;
    ae_int_t s1;
    ae_int_t s2;
    ae_int_t s3;
    ae_int_t s4;
    ae_int_t sf;
    ae_int_t sfx;
    ae_int_t sfy;
    ae_int_t sfxy;
    double y1;
    double y2;
    double y3;
    double y4;
    double v;
    double t0;
    double t1;
    double t2;
    double t3;
    double u0;
    double u1;
    double u2;
    double u3;

    *f = 0;
    *fx = 0;
    *fy = 0;
    *fxy = 0;

    ae_assert(ae_round(c->c.ptr.p_double[1], _state)==-1||ae_round(c->c.ptr.p_double[1], _state)==-3, "Spline2DDiff: incorrect C!", _state);
    n = ae_round(c->c.ptr.p_double[2], _state);
    m = ae_round(c->c.ptr.p_double[3], _state);
    
    /*
     * Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
     */
    l = 4;
    r = 4+n-2+1;
    while(l!=r-1)
    {
        h = (l+r)/2;
        if( ae_fp_greater_eq(c->c.ptr.p_double[h],x) )
        {
            r = h;
        }
        else
        {
            l = h;
        }
    }
    t = (x-c->c.ptr.p_double[l])/(c->c.ptr.p_double[l+1]-c->c.ptr.p_double[l]);
    dt = 1.0/(c->c.ptr.p_double[l+1]-c->c.ptr.p_double[l]);
    ix = l-4;
    
    /*
     * Binary search in the [ y[0], ..., y[m-2] ] (y[m-1] is not included)
     */
    l = 4+n;
    r = 4+n+(m-2)+1;
    while(l!=r-1)
    {
        h = (l+r)/2;
        if( ae_fp_greater_eq(c->c.ptr.p_double[h],y) )
        {
            r = h;
        }
        else
        {
            l = h;
        }
    }
    u = (y-c->c.ptr.p_double[l])/(c->c.ptr.p_double[l+1]-c->c.ptr.p_double[l]);
    du = 1.0/(c->c.ptr.p_double[l+1]-c->c.ptr.p_double[l]);
    iy = l-(4+n);
    
    /*
     * Prepare F, dF/dX, dF/dY, d2F/dXdY
     */
    *f = 0;
    *fx = 0;
    *fy = 0;
    *fxy = 0;
    
    /*
     * Bilinear interpolation
     */
    if( ae_round(c->c.ptr.p_double[1], _state)==-1 )
    {
        shift1 = 4+n+m;
        y1 = c->c.ptr.p_double[shift1+n*iy+ix];
        y2 = c->c.ptr.p_double[shift1+n*iy+(ix+1)];
        y3 = c->c.ptr.p_double[shift1+n*(iy+1)+(ix+1)];
        y4 = c->c.ptr.p_double[shift1+n*(iy+1)+ix];
        *f = (1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
        *fx = (-(1-u)*y1+(1-u)*y2+u*y3-u*y4)*dt;
        *fy = (-(1-t)*y1-t*y2+t*y3+(1-t)*y4)*du;
        *fxy = (y1-y2+y3-y4)*du*dt;
        return;
    }
    
    /*
     * Bicubic interpolation
     */
    if( ae_round(c->c.ptr.p_double[1], _state)==-3 )
    {
        
        /*
         * Prepare info
         */
        t0 = 1;
        t1 = t;
        t2 = ae_sqr(t, _state);
        t3 = t*t2;
        u0 = 1;
        u1 = u;
        u2 = ae_sqr(u, _state);
        u3 = u*u2;
        sf = 4+n+m;
        sfx = 4+n+m+n*m;
        sfy = 4+n+m+2*n*m;
        sfxy = 4+n+m+3*n*m;
        s1 = n*iy+ix;
        s2 = n*iy+(ix+1);
        s3 = n*(iy+1)+(ix+1);
        s4 = n*(iy+1)+ix;
        
        /*
         * Calculate
         */
        v = 1*c->c.ptr.p_double[sf+s1];
        *f = *f+v*t0*u0;
        v = 1*c->c.ptr.p_double[sfy+s1]/du;
        *f = *f+v*t0*u1;
        *fy = *fy+1*v*t0*u0*du;
        v = -3*c->c.ptr.p_double[sf+s1]+3*c->c.ptr.p_double[sf+s4]-2*c->c.ptr.p_double[sfy+s1]/du-1*c->c.ptr.p_double[sfy+s4]/du;
        *f = *f+v*t0*u2;
        *fy = *fy+2*v*t0*u1*du;
        v = 2*c->c.ptr.p_double[sf+s1]-2*c->c.ptr.p_double[sf+s4]+1*c->c.ptr.p_double[sfy+s1]/du+1*c->c.ptr.p_double[sfy+s4]/du;
        *f = *f+v*t0*u3;
        *fy = *fy+3*v*t0*u2*du;
        v = 1*c->c.ptr.p_double[sfx+s1]/dt;
        *f = *f+v*t1*u0;
        *fx = *fx+1*v*t0*u0*dt;
        v = 1*c->c.ptr.p_double[sfxy+s1]/(dt*du);
        *f = *f+v*t1*u1;
        *fx = *fx+1*v*t0*u1*dt;
        *fy = *fy+1*v*t1*u0*du;
        *fxy = *fxy+1*v*t0*u0*dt*du;
        v = -3*c->c.ptr.p_double[sfx+s1]/dt+3*c->c.ptr.p_double[sfx+s4]/dt-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t1*u2;
        *fx = *fx+1*v*t0*u2*dt;
        *fy = *fy+2*v*t1*u1*du;
        *fxy = *fxy+2*v*t0*u1*dt*du;
        v = 2*c->c.ptr.p_double[sfx+s1]/dt-2*c->c.ptr.p_double[sfx+s4]/dt+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t1*u3;
        *fx = *fx+1*v*t0*u3*dt;
        *fy = *fy+3*v*t1*u2*du;
        *fxy = *fxy+3*v*t0*u2*dt*du;
        v = -3*c->c.ptr.p_double[sf+s1]+3*c->c.ptr.p_double[sf+s2]-2*c->c.ptr.p_double[sfx+s1]/dt-1*c->c.ptr.p_double[sfx+s2]/dt;
        *f = *f+v*t2*u0;
        *fx = *fx+2*v*t1*u0*dt;
        v = -3*c->c.ptr.p_double[sfy+s1]/du+3*c->c.ptr.p_double[sfy+s2]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s2]/(dt*du);
        *f = *f+v*t2*u1;
        *fx = *fx+2*v*t1*u1*dt;
        *fy = *fy+1*v*t2*u0*du;
        *fxy = *fxy+2*v*t1*u0*dt*du;
        v = 9*c->c.ptr.p_double[sf+s1]-9*c->c.ptr.p_double[sf+s2]+9*c->c.ptr.p_double[sf+s3]-9*c->c.ptr.p_double[sf+s4]+6*c->c.ptr.p_double[sfx+s1]/dt+3*c->c.ptr.p_double[sfx+s2]/dt-3*c->c.ptr.p_double[sfx+s3]/dt-6*c->c.ptr.p_double[sfx+s4]/dt+6*c->c.ptr.p_double[sfy+s1]/du-6*c->c.ptr.p_double[sfy+s2]/du-3*c->c.ptr.p_double[sfy+s3]/du+3*c->c.ptr.p_double[sfy+s4]/du+4*c->c.ptr.p_double[sfxy+s1]/(dt*du)+2*c->c.ptr.p_double[sfxy+s2]/(dt*du)+1*c->c.ptr.p_double[sfxy+s3]/(dt*du)+2*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t2*u2;
        *fx = *fx+2*v*t1*u2*dt;
        *fy = *fy+2*v*t2*u1*du;
        *fxy = *fxy+4*v*t1*u1*dt*du;
        v = -6*c->c.ptr.p_double[sf+s1]+6*c->c.ptr.p_double[sf+s2]-6*c->c.ptr.p_double[sf+s3]+6*c->c.ptr.p_double[sf+s4]-4*c->c.ptr.p_double[sfx+s1]/dt-2*c->c.ptr.p_double[sfx+s2]/dt+2*c->c.ptr.p_double[sfx+s3]/dt+4*c->c.ptr.p_double[sfx+s4]/dt-3*c->c.ptr.p_double[sfy+s1]/du+3*c->c.ptr.p_double[sfy+s2]/du+3*c->c.ptr.p_double[sfy+s3]/du-3*c->c.ptr.p_double[sfy+s4]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s2]/(dt*du)-1*c->c.ptr.p_double[sfxy+s3]/(dt*du)-2*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t2*u3;
        *fx = *fx+2*v*t1*u3*dt;
        *fy = *fy+3*v*t2*u2*du;
        *fxy = *fxy+6*v*t1*u2*dt*du;
        v = 2*c->c.ptr.p_double[sf+s1]-2*c->c.ptr.p_double[sf+s2]+1*c->c.ptr.p_double[sfx+s1]/dt+1*c->c.ptr.p_double[sfx+s2]/dt;
        *f = *f+v*t3*u0;
        *fx = *fx+3*v*t2*u0*dt;
        v = 2*c->c.ptr.p_double[sfy+s1]/du-2*c->c.ptr.p_double[sfy+s2]/du+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s2]/(dt*du);
        *f = *f+v*t3*u1;
        *fx = *fx+3*v*t2*u1*dt;
        *fy = *fy+1*v*t3*u0*du;
        *fxy = *fxy+3*v*t2*u0*dt*du;
        v = -6*c->c.ptr.p_double[sf+s1]+6*c->c.ptr.p_double[sf+s2]-6*c->c.ptr.p_double[sf+s3]+6*c->c.ptr.p_double[sf+s4]-3*c->c.ptr.p_double[sfx+s1]/dt-3*c->c.ptr.p_double[sfx+s2]/dt+3*c->c.ptr.p_double[sfx+s3]/dt+3*c->c.ptr.p_double[sfx+s4]/dt-4*c->c.ptr.p_double[sfy+s1]/du+4*c->c.ptr.p_double[sfy+s2]/du+2*c->c.ptr.p_double[sfy+s3]/du-2*c->c.ptr.p_double[sfy+s4]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-2*c->c.ptr.p_double[sfxy+s2]/(dt*du)-1*c->c.ptr.p_double[sfxy+s3]/(dt*du)-1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t3*u2;
        *fx = *fx+3*v*t2*u2*dt;
        *fy = *fy+2*v*t3*u1*du;
        *fxy = *fxy+6*v*t2*u1*dt*du;
        v = 4*c->c.ptr.p_double[sf+s1]-4*c->c.ptr.p_double[sf+s2]+4*c->c.ptr.p_double[sf+s3]-4*c->c.ptr.p_double[sf+s4]+2*c->c.ptr.p_double[sfx+s1]/dt+2*c->c.ptr.p_double[sfx+s2]/dt-2*c->c.ptr.p_double[sfx+s3]/dt-2*c->c.ptr.p_double[sfx+s4]/dt+2*c->c.ptr.p_double[sfy+s1]/du-2*c->c.ptr.p_double[sfy+s2]/du-2*c->c.ptr.p_double[sfy+s3]/du+2*c->c.ptr.p_double[sfy+s4]/du+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s2]/(dt*du)+1*c->c.ptr.p_double[sfxy+s3]/(dt*du)+1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
        *f = *f+v*t3*u3;
        *fx = *fx+3*v*t2*u3*dt;
        *fy = *fy+3*v*t3*u2*du;
        *fxy = *fxy+9*v*t2*u2*dt*du;
        return;
    }
}


/*************************************************************************
This subroutine unpacks two-dimensional spline into the coefficients table

Input parameters:
    C   -   spline interpolant.

Result:
    M, N-   grid size (x-axis and y-axis)
    Tbl -   coefficients table, unpacked format,
            [0..(N-1)*(M-1)-1, 0..19].
            For I = 0...M-2, J=0..N-2:
                K =  I*(N-1)+J
                Tbl[K,0] = X[j]
                Tbl[K,1] = X[j+1]
                Tbl[K,2] = Y[i]
                Tbl[K,3] = Y[i+1]
                Tbl[K,4] = C00
                Tbl[K,5] = C01
                Tbl[K,6] = C02
                Tbl[K,7] = C03
                Tbl[K,8] = C10
                Tbl[K,9] = C11
                ...
                Tbl[K,19] = C33
            On each grid square spline is equals to:
                S(x) = SUM(c[i,j]*(x^i)*(y^j), i=0..3, j=0..3)
                t = x-x[j]
                u = y-y[i]

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dunpack(spline2dinterpolant* c,
     ae_int_t* m,
     ae_int_t* n,
     /* Real    */ ae_matrix* tbl,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t ci;
    ae_int_t cj;
    ae_int_t k;
    ae_int_t p;
    ae_int_t shift;
    ae_int_t s1;
    ae_int_t s2;
    ae_int_t s3;
    ae_int_t s4;
    ae_int_t sf;
    ae_int_t sfx;
    ae_int_t sfy;
    ae_int_t sfxy;
    double y1;
    double y2;
    double y3;
    double y4;
    double dt;
    double du;

    *m = 0;
    *n = 0;
    ae_matrix_clear(tbl);

    ae_assert(ae_round(c->c.ptr.p_double[1], _state)==-3||ae_round(c->c.ptr.p_double[1], _state)==-1, "SplineUnpack2D: incorrect C!", _state);
    *n = ae_round(c->c.ptr.p_double[2], _state);
    *m = ae_round(c->c.ptr.p_double[3], _state);
    ae_matrix_set_length(tbl, (*n-1)*(*m-1)-1+1, 19+1, _state);
    
    /*
     * Fill
     */
    for(i=0; i<=*m-2; i++)
    {
        for(j=0; j<=*n-2; j++)
        {
            p = i*(*n-1)+j;
            tbl->ptr.pp_double[p][0] = c->c.ptr.p_double[4+j];
            tbl->ptr.pp_double[p][1] = c->c.ptr.p_double[4+j+1];
            tbl->ptr.pp_double[p][2] = c->c.ptr.p_double[4+(*n)+i];
            tbl->ptr.pp_double[p][3] = c->c.ptr.p_double[4+(*n)+i+1];
            dt = 1/(tbl->ptr.pp_double[p][1]-tbl->ptr.pp_double[p][0]);
            du = 1/(tbl->ptr.pp_double[p][3]-tbl->ptr.pp_double[p][2]);
            
            /*
             * Bilinear interpolation
             */
            if( ae_round(c->c.ptr.p_double[1], _state)==-1 )
            {
                for(k=4; k<=19; k++)
                {
                    tbl->ptr.pp_double[p][k] = 0;
                }
                shift = 4+(*n)+(*m);
                y1 = c->c.ptr.p_double[shift+*n*i+j];
                y2 = c->c.ptr.p_double[shift+*n*i+(j+1)];
                y3 = c->c.ptr.p_double[shift+*n*(i+1)+(j+1)];
                y4 = c->c.ptr.p_double[shift+*n*(i+1)+j];
                tbl->ptr.pp_double[p][4] = y1;
                tbl->ptr.pp_double[p][4+1*4+0] = y2-y1;
                tbl->ptr.pp_double[p][4+0*4+1] = y4-y1;
                tbl->ptr.pp_double[p][4+1*4+1] = y3-y2-y4+y1;
            }
            
            /*
             * Bicubic interpolation
             */
            if( ae_round(c->c.ptr.p_double[1], _state)==-3 )
            {
                sf = 4+(*n)+(*m);
                sfx = 4+(*n)+(*m)+*n*(*m);
                sfy = 4+(*n)+(*m)+2*(*n)*(*m);
                sfxy = 4+(*n)+(*m)+3*(*n)*(*m);
                s1 = *n*i+j;
                s2 = *n*i+(j+1);
                s3 = *n*(i+1)+(j+1);
                s4 = *n*(i+1)+j;
                tbl->ptr.pp_double[p][4+0*4+0] = 1*c->c.ptr.p_double[sf+s1];
                tbl->ptr.pp_double[p][4+0*4+1] = 1*c->c.ptr.p_double[sfy+s1]/du;
                tbl->ptr.pp_double[p][4+0*4+2] = -3*c->c.ptr.p_double[sf+s1]+3*c->c.ptr.p_double[sf+s4]-2*c->c.ptr.p_double[sfy+s1]/du-1*c->c.ptr.p_double[sfy+s4]/du;
                tbl->ptr.pp_double[p][4+0*4+3] = 2*c->c.ptr.p_double[sf+s1]-2*c->c.ptr.p_double[sf+s4]+1*c->c.ptr.p_double[sfy+s1]/du+1*c->c.ptr.p_double[sfy+s4]/du;
                tbl->ptr.pp_double[p][4+1*4+0] = 1*c->c.ptr.p_double[sfx+s1]/dt;
                tbl->ptr.pp_double[p][4+1*4+1] = 1*c->c.ptr.p_double[sfxy+s1]/(dt*du);
                tbl->ptr.pp_double[p][4+1*4+2] = -3*c->c.ptr.p_double[sfx+s1]/dt+3*c->c.ptr.p_double[sfx+s4]/dt-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
                tbl->ptr.pp_double[p][4+1*4+3] = 2*c->c.ptr.p_double[sfx+s1]/dt-2*c->c.ptr.p_double[sfx+s4]/dt+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
                tbl->ptr.pp_double[p][4+2*4+0] = -3*c->c.ptr.p_double[sf+s1]+3*c->c.ptr.p_double[sf+s2]-2*c->c.ptr.p_double[sfx+s1]/dt-1*c->c.ptr.p_double[sfx+s2]/dt;
                tbl->ptr.pp_double[p][4+2*4+1] = -3*c->c.ptr.p_double[sfy+s1]/du+3*c->c.ptr.p_double[sfy+s2]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s2]/(dt*du);
                tbl->ptr.pp_double[p][4+2*4+2] = 9*c->c.ptr.p_double[sf+s1]-9*c->c.ptr.p_double[sf+s2]+9*c->c.ptr.p_double[sf+s3]-9*c->c.ptr.p_double[sf+s4]+6*c->c.ptr.p_double[sfx+s1]/dt+3*c->c.ptr.p_double[sfx+s2]/dt-3*c->c.ptr.p_double[sfx+s3]/dt-6*c->c.ptr.p_double[sfx+s4]/dt+6*c->c.ptr.p_double[sfy+s1]/du-6*c->c.ptr.p_double[sfy+s2]/du-3*c->c.ptr.p_double[sfy+s3]/du+3*c->c.ptr.p_double[sfy+s4]/du+4*c->c.ptr.p_double[sfxy+s1]/(dt*du)+2*c->c.ptr.p_double[sfxy+s2]/(dt*du)+1*c->c.ptr.p_double[sfxy+s3]/(dt*du)+2*c->c.ptr.p_double[sfxy+s4]/(dt*du);
                tbl->ptr.pp_double[p][4+2*4+3] = -6*c->c.ptr.p_double[sf+s1]+6*c->c.ptr.p_double[sf+s2]-6*c->c.ptr.p_double[sf+s3]+6*c->c.ptr.p_double[sf+s4]-4*c->c.ptr.p_double[sfx+s1]/dt-2*c->c.ptr.p_double[sfx+s2]/dt+2*c->c.ptr.p_double[sfx+s3]/dt+4*c->c.ptr.p_double[sfx+s4]/dt-3*c->c.ptr.p_double[sfy+s1]/du+3*c->c.ptr.p_double[sfy+s2]/du+3*c->c.ptr.p_double[sfy+s3]/du-3*c->c.ptr.p_double[sfy+s4]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-1*c->c.ptr.p_double[sfxy+s2]/(dt*du)-1*c->c.ptr.p_double[sfxy+s3]/(dt*du)-2*c->c.ptr.p_double[sfxy+s4]/(dt*du);
                tbl->ptr.pp_double[p][4+3*4+0] = 2*c->c.ptr.p_double[sf+s1]-2*c->c.ptr.p_double[sf+s2]+1*c->c.ptr.p_double[sfx+s1]/dt+1*c->c.ptr.p_double[sfx+s2]/dt;
                tbl->ptr.pp_double[p][4+3*4+1] = 2*c->c.ptr.p_double[sfy+s1]/du-2*c->c.ptr.p_double[sfy+s2]/du+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s2]/(dt*du);
                tbl->ptr.pp_double[p][4+3*4+2] = -6*c->c.ptr.p_double[sf+s1]+6*c->c.ptr.p_double[sf+s2]-6*c->c.ptr.p_double[sf+s3]+6*c->c.ptr.p_double[sf+s4]-3*c->c.ptr.p_double[sfx+s1]/dt-3*c->c.ptr.p_double[sfx+s2]/dt+3*c->c.ptr.p_double[sfx+s3]/dt+3*c->c.ptr.p_double[sfx+s4]/dt-4*c->c.ptr.p_double[sfy+s1]/du+4*c->c.ptr.p_double[sfy+s2]/du+2*c->c.ptr.p_double[sfy+s3]/du-2*c->c.ptr.p_double[sfy+s4]/du-2*c->c.ptr.p_double[sfxy+s1]/(dt*du)-2*c->c.ptr.p_double[sfxy+s2]/(dt*du)-1*c->c.ptr.p_double[sfxy+s3]/(dt*du)-1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
                tbl->ptr.pp_double[p][4+3*4+3] = 4*c->c.ptr.p_double[sf+s1]-4*c->c.ptr.p_double[sf+s2]+4*c->c.ptr.p_double[sf+s3]-4*c->c.ptr.p_double[sf+s4]+2*c->c.ptr.p_double[sfx+s1]/dt+2*c->c.ptr.p_double[sfx+s2]/dt-2*c->c.ptr.p_double[sfx+s3]/dt-2*c->c.ptr.p_double[sfx+s4]/dt+2*c->c.ptr.p_double[sfy+s1]/du-2*c->c.ptr.p_double[sfy+s2]/du-2*c->c.ptr.p_double[sfy+s3]/du+2*c->c.ptr.p_double[sfy+s4]/du+1*c->c.ptr.p_double[sfxy+s1]/(dt*du)+1*c->c.ptr.p_double[sfxy+s2]/(dt*du)+1*c->c.ptr.p_double[sfxy+s3]/(dt*du)+1*c->c.ptr.p_double[sfxy+s4]/(dt*du);
            }
            
            /*
             * Rescale Cij
             */
            for(ci=0; ci<=3; ci++)
            {
                for(cj=0; cj<=3; cj++)
                {
                    tbl->ptr.pp_double[p][4+ci*4+cj] = tbl->ptr.pp_double[p][4+ci*4+cj]*ae_pow(dt, ci, _state)*ae_pow(du, cj, _state);
                }
            }
        }
    }
}


/*************************************************************************
This subroutine performs linear transformation of the spline argument.

Input parameters:
    C       -   spline interpolant
    AX, BX  -   transformation coefficients: x = A*t + B
    AY, BY  -   transformation coefficients: y = A*u + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dlintransxy(spline2dinterpolant* c,
     double ax,
     double bx,
     double ay,
     double by,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    ae_int_t m;
    double v;
    ae_vector x;
    ae_vector y;
    ae_matrix f;
    ae_int_t typec;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&f, 0, 0, DT_REAL, _state, ae_true);

    typec = ae_round(c->c.ptr.p_double[1], _state);
    ae_assert(typec==-3||typec==-1, "Spline2DLinTransXY: incorrect C!", _state);
    n = ae_round(c->c.ptr.p_double[2], _state);
    m = ae_round(c->c.ptr.p_double[3], _state);
    ae_vector_set_length(&x, n-1+1, _state);
    ae_vector_set_length(&y, m-1+1, _state);
    ae_matrix_set_length(&f, m-1+1, n-1+1, _state);
    for(j=0; j<=n-1; j++)
    {
        x.ptr.p_double[j] = c->c.ptr.p_double[4+j];
    }
    for(i=0; i<=m-1; i++)
    {
        y.ptr.p_double[i] = c->c.ptr.p_double[4+n+i];
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            f.ptr.pp_double[i][j] = c->c.ptr.p_double[4+n+m+i*n+j];
        }
    }
    
    /*
     * Special case: AX=0 or AY=0
     */
    if( ae_fp_eq(ax,0) )
    {
        for(i=0; i<=m-1; i++)
        {
            v = spline2dcalc(c, bx, y.ptr.p_double[i], _state);
            for(j=0; j<=n-1; j++)
            {
                f.ptr.pp_double[i][j] = v;
            }
        }
        if( typec==-3 )
        {
            spline2dbuildbicubic(&x, &y, &f, m, n, c, _state);
        }
        if( typec==-1 )
        {
            spline2dbuildbilinear(&x, &y, &f, m, n, c, _state);
        }
        ax = 1;
        bx = 0;
    }
    if( ae_fp_eq(ay,0) )
    {
        for(j=0; j<=n-1; j++)
        {
            v = spline2dcalc(c, x.ptr.p_double[j], by, _state);
            for(i=0; i<=m-1; i++)
            {
                f.ptr.pp_double[i][j] = v;
            }
        }
        if( typec==-3 )
        {
            spline2dbuildbicubic(&x, &y, &f, m, n, c, _state);
        }
        if( typec==-1 )
        {
            spline2dbuildbilinear(&x, &y, &f, m, n, c, _state);
        }
        ay = 1;
        by = 0;
    }
    
    /*
     * General case: AX<>0, AY<>0
     * Unpack, scale and pack again.
     */
    for(j=0; j<=n-1; j++)
    {
        x.ptr.p_double[j] = (x.ptr.p_double[j]-bx)/ax;
    }
    for(i=0; i<=m-1; i++)
    {
        y.ptr.p_double[i] = (y.ptr.p_double[i]-by)/ay;
    }
    if( typec==-3 )
    {
        spline2dbuildbicubic(&x, &y, &f, m, n, c, _state);
    }
    if( typec==-1 )
    {
        spline2dbuildbilinear(&x, &y, &f, m, n, c, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine performs linear transformation of the spline.

Input parameters:
    C   -   spline interpolant.
    A, B-   transformation coefficients: S2(x,y) = A*S(x,y) + B
    
Output parameters:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dlintransf(spline2dinterpolant* c,
     double a,
     double b,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    ae_int_t m;
    ae_vector x;
    ae_vector y;
    ae_matrix f;
    ae_int_t typec;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&f, 0, 0, DT_REAL, _state, ae_true);

    typec = ae_round(c->c.ptr.p_double[1], _state);
    ae_assert(typec==-3||typec==-1, "Spline2DLinTransXY: incorrect C!", _state);
    n = ae_round(c->c.ptr.p_double[2], _state);
    m = ae_round(c->c.ptr.p_double[3], _state);
    ae_vector_set_length(&x, n-1+1, _state);
    ae_vector_set_length(&y, m-1+1, _state);
    ae_matrix_set_length(&f, m-1+1, n-1+1, _state);
    for(j=0; j<=n-1; j++)
    {
        x.ptr.p_double[j] = c->c.ptr.p_double[4+j];
    }
    for(i=0; i<=m-1; i++)
    {
        y.ptr.p_double[i] = c->c.ptr.p_double[4+n+i];
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            f.ptr.pp_double[i][j] = a*c->c.ptr.p_double[4+n+m+i*n+j]+b;
        }
    }
    if( typec==-3 )
    {
        spline2dbuildbicubic(&x, &y, &f, m, n, c, _state);
    }
    if( typec==-1 )
    {
        spline2dbuildbilinear(&x, &y, &f, m, n, c, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine makes the copy of the spline model.

Input parameters:
    C   -   spline interpolant

Output parameters:
    CC  -   spline copy

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline2dcopy(spline2dinterpolant* c,
     spline2dinterpolant* cc,
     ae_state *_state)
{
    ae_int_t n;

    _spline2dinterpolant_clear(cc);

    ae_assert(c->k==1||c->k==3, "Spline2DCopy: incorrect C!", _state);
    cc->k = c->k;
    n = ae_round(c->c.ptr.p_double[0], _state);
    ae_vector_set_length(&cc->c, n, _state);
    ae_v_move(&cc->c.ptr.p_double[0], 1, &c->c.ptr.p_double[0], 1, ae_v_len(0,n-1));
}


/*************************************************************************
Bicubic spline resampling

Input parameters:
    A           -   function values at the old grid,
                    array[0..OldHeight-1, 0..OldWidth-1]
    OldHeight   -   old grid height, OldHeight>1
    OldWidth    -   old grid width, OldWidth>1
    NewHeight   -   new grid height, NewHeight>1
    NewWidth    -   new grid width, NewWidth>1
    
Output parameters:
    B           -   function values at the new grid,
                    array[0..NewHeight-1, 0..NewWidth-1]

  -- ALGLIB routine --
     15 May, 2007
     Copyright by Bochkanov Sergey
*************************************************************************/
void spline2dresamplebicubic(/* Real    */ ae_matrix* a,
     ae_int_t oldheight,
     ae_int_t oldwidth,
     /* Real    */ ae_matrix* b,
     ae_int_t newheight,
     ae_int_t newwidth,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix buf;
    ae_vector x;
    ae_vector y;
    spline1dinterpolant c;
    ae_int_t i;
    ae_int_t j;
    ae_int_t mw;
    ae_int_t mh;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_clear(b);
    ae_matrix_init(&buf, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    _spline1dinterpolant_init(&c, _state, ae_true);

    ae_assert(oldwidth>1&&oldheight>1, "Spline2DResampleBicubic: width/height less than 1", _state);
    ae_assert(newwidth>1&&newheight>1, "Spline2DResampleBicubic: width/height less than 1", _state);
    
    /*
     * Prepare
     */
    mw = ae_maxint(oldwidth, newwidth, _state);
    mh = ae_maxint(oldheight, newheight, _state);
    ae_matrix_set_length(b, newheight-1+1, newwidth-1+1, _state);
    ae_matrix_set_length(&buf, oldheight-1+1, newwidth-1+1, _state);
    ae_vector_set_length(&x, ae_maxint(mw, mh, _state)-1+1, _state);
    ae_vector_set_length(&y, ae_maxint(mw, mh, _state)-1+1, _state);
    
    /*
     * Horizontal interpolation
     */
    for(i=0; i<=oldheight-1; i++)
    {
        
        /*
         * Fill X, Y
         */
        for(j=0; j<=oldwidth-1; j++)
        {
            x.ptr.p_double[j] = (double)j/(double)(oldwidth-1);
            y.ptr.p_double[j] = a->ptr.pp_double[i][j];
        }
        
        /*
         * Interpolate and place result into temporary matrix
         */
        spline1dbuildcubic(&x, &y, oldwidth, 0, 0.0, 0, 0.0, &c, _state);
        for(j=0; j<=newwidth-1; j++)
        {
            buf.ptr.pp_double[i][j] = spline1dcalc(&c, (double)j/(double)(newwidth-1), _state);
        }
    }
    
    /*
     * Vertical interpolation
     */
    for(j=0; j<=newwidth-1; j++)
    {
        
        /*
         * Fill X, Y
         */
        for(i=0; i<=oldheight-1; i++)
        {
            x.ptr.p_double[i] = (double)i/(double)(oldheight-1);
            y.ptr.p_double[i] = buf.ptr.pp_double[i][j];
        }
        
        /*
         * Interpolate and place result into B
         */
        spline1dbuildcubic(&x, &y, oldheight, 0, 0.0, 0, 0.0, &c, _state);
        for(i=0; i<=newheight-1; i++)
        {
            b->ptr.pp_double[i][j] = spline1dcalc(&c, (double)i/(double)(newheight-1), _state);
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Bilinear spline resampling

Input parameters:
    A           -   function values at the old grid,
                    array[0..OldHeight-1, 0..OldWidth-1]
    OldHeight   -   old grid height, OldHeight>1
    OldWidth    -   old grid width, OldWidth>1
    NewHeight   -   new grid height, NewHeight>1
    NewWidth    -   new grid width, NewWidth>1

Output parameters:
    B           -   function values at the new grid,
                    array[0..NewHeight-1, 0..NewWidth-1]

  -- ALGLIB routine --
     09.07.2007
     Copyright by Bochkanov Sergey
*************************************************************************/
void spline2dresamplebilinear(/* Real    */ ae_matrix* a,
     ae_int_t oldheight,
     ae_int_t oldwidth,
     /* Real    */ ae_matrix* b,
     ae_int_t newheight,
     ae_int_t newwidth,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t l;
    ae_int_t c;
    double t;
    double u;

    ae_matrix_clear(b);

    ae_matrix_set_length(b, newheight-1+1, newwidth-1+1, _state);
    for(i=0; i<=newheight-1; i++)
    {
        for(j=0; j<=newwidth-1; j++)
        {
            l = i*(oldheight-1)/(newheight-1);
            if( l==oldheight-1 )
            {
                l = oldheight-2;
            }
            u = (double)i/(double)(newheight-1)*(oldheight-1)-l;
            c = j*(oldwidth-1)/(newwidth-1);
            if( c==oldwidth-1 )
            {
                c = oldwidth-2;
            }
            t = (double)(j*(oldwidth-1))/(double)(newwidth-1)-c;
            b->ptr.pp_double[i][j] = (1-t)*(1-u)*a->ptr.pp_double[l][c]+t*(1-u)*a->ptr.pp_double[l][c+1]+t*u*a->ptr.pp_double[l+1][c+1]+(1-t)*u*a->ptr.pp_double[l+1][c];
        }
    }
}


/*************************************************************************
Internal subroutine.
Calculation of the first derivatives and the cross-derivative.
*************************************************************************/
static void spline2d_bicubiccalcderivatives(/* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* dx,
     /* Real    */ ae_matrix* dy,
     /* Real    */ ae_matrix* dxy,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_vector xt;
    ae_vector ft;
    double s;
    double ds;
    double d2s;
    spline1dinterpolant c;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_clear(dx);
    ae_matrix_clear(dy);
    ae_matrix_clear(dxy);
    ae_vector_init(&xt, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ft, 0, DT_REAL, _state, ae_true);
    _spline1dinterpolant_init(&c, _state, ae_true);

    ae_matrix_set_length(dx, m, n, _state);
    ae_matrix_set_length(dy, m, n, _state);
    ae_matrix_set_length(dxy, m, n, _state);
    
    /*
     * dF/dX
     */
    ae_vector_set_length(&xt, n, _state);
    ae_vector_set_length(&ft, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            xt.ptr.p_double[j] = x->ptr.p_double[j];
            ft.ptr.p_double[j] = a->ptr.pp_double[i][j];
        }
        spline1dbuildcubic(&xt, &ft, n, 0, 0.0, 0, 0.0, &c, _state);
        for(j=0; j<=n-1; j++)
        {
            spline1ddiff(&c, x->ptr.p_double[j], &s, &ds, &d2s, _state);
            dx->ptr.pp_double[i][j] = ds;
        }
    }
    
    /*
     * dF/dY
     */
    ae_vector_set_length(&xt, m, _state);
    ae_vector_set_length(&ft, m, _state);
    for(j=0; j<=n-1; j++)
    {
        for(i=0; i<=m-1; i++)
        {
            xt.ptr.p_double[i] = y->ptr.p_double[i];
            ft.ptr.p_double[i] = a->ptr.pp_double[i][j];
        }
        spline1dbuildcubic(&xt, &ft, m, 0, 0.0, 0, 0.0, &c, _state);
        for(i=0; i<=m-1; i++)
        {
            spline1ddiff(&c, y->ptr.p_double[i], &s, &ds, &d2s, _state);
            dy->ptr.pp_double[i][j] = ds;
        }
    }
    
    /*
     * d2F/dXdY
     */
    ae_vector_set_length(&xt, n, _state);
    ae_vector_set_length(&ft, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            xt.ptr.p_double[j] = x->ptr.p_double[j];
            ft.ptr.p_double[j] = dy->ptr.pp_double[i][j];
        }
        spline1dbuildcubic(&xt, &ft, n, 0, 0.0, 0, 0.0, &c, _state);
        for(j=0; j<=n-1; j++)
        {
            spline1ddiff(&c, x->ptr.p_double[j], &s, &ds, &d2s, _state);
            dxy->ptr.pp_double[i][j] = ds;
        }
    }
    ae_frame_leave(_state);
}


ae_bool _spline2dinterpolant_init(spline2dinterpolant* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->c, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _spline2dinterpolant_init_copy(spline2dinterpolant* dst, spline2dinterpolant* src, ae_state *_state, ae_bool make_automatic)
{
    dst->k = src->k;
    if( !ae_vector_init_copy(&dst->c, &src->c, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _spline2dinterpolant_clear(spline2dinterpolant* p)
{
    ae_vector_clear(&p->c);
}


/*$ End $*/
