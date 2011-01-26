/*************************************************************************
Copyright 2010 by Sergey Bochkanov (ALGLIB project).

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
#include "basicstatops.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Internal ranking subroutine
*************************************************************************/
void rankx(/* Real    */ ae_vector* x,
     ae_int_t n,
     apbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t t;
    double tmp;
    ae_int_t tmpi;


    
    /*
     * Prepare
     */
    if( n<1 )
    {
        return;
    }
    if( n==1 )
    {
        x->ptr.p_double[0] = 1;
        return;
    }
    if( buf->ra1.cnt<n )
    {
        ae_vector_set_length(&buf->ra1, n, _state);
    }
    if( buf->ia1.cnt<n )
    {
        ae_vector_set_length(&buf->ia1, n, _state);
    }
    for(i=0; i<=n-1; i++)
    {
        buf->ra1.ptr.p_double[i] = x->ptr.p_double[i];
        buf->ia1.ptr.p_int[i] = i;
    }
    
    /*
     * sort {R, C}
     */
    if( n!=1 )
    {
        i = 2;
        do
        {
            t = i;
            while(t!=1)
            {
                k = t/2;
                if( ae_fp_greater_eq(buf->ra1.ptr.p_double[k-1],buf->ra1.ptr.p_double[t-1]) )
                {
                    t = 1;
                }
                else
                {
                    tmp = buf->ra1.ptr.p_double[k-1];
                    buf->ra1.ptr.p_double[k-1] = buf->ra1.ptr.p_double[t-1];
                    buf->ra1.ptr.p_double[t-1] = tmp;
                    tmpi = buf->ia1.ptr.p_int[k-1];
                    buf->ia1.ptr.p_int[k-1] = buf->ia1.ptr.p_int[t-1];
                    buf->ia1.ptr.p_int[t-1] = tmpi;
                    t = k;
                }
            }
            i = i+1;
        }
        while(i<=n);
        i = n-1;
        do
        {
            tmp = buf->ra1.ptr.p_double[i];
            buf->ra1.ptr.p_double[i] = buf->ra1.ptr.p_double[0];
            buf->ra1.ptr.p_double[0] = tmp;
            tmpi = buf->ia1.ptr.p_int[i];
            buf->ia1.ptr.p_int[i] = buf->ia1.ptr.p_int[0];
            buf->ia1.ptr.p_int[0] = tmpi;
            t = 1;
            while(t!=0)
            {
                k = 2*t;
                if( k>i )
                {
                    t = 0;
                }
                else
                {
                    if( k<i )
                    {
                        if( ae_fp_greater(buf->ra1.ptr.p_double[k],buf->ra1.ptr.p_double[k-1]) )
                        {
                            k = k+1;
                        }
                    }
                    if( ae_fp_greater_eq(buf->ra1.ptr.p_double[t-1],buf->ra1.ptr.p_double[k-1]) )
                    {
                        t = 0;
                    }
                    else
                    {
                        tmp = buf->ra1.ptr.p_double[k-1];
                        buf->ra1.ptr.p_double[k-1] = buf->ra1.ptr.p_double[t-1];
                        buf->ra1.ptr.p_double[t-1] = tmp;
                        tmpi = buf->ia1.ptr.p_int[k-1];
                        buf->ia1.ptr.p_int[k-1] = buf->ia1.ptr.p_int[t-1];
                        buf->ia1.ptr.p_int[t-1] = tmpi;
                        t = k;
                    }
                }
            }
            i = i-1;
        }
        while(i>=1);
    }
    
    /*
     * compute tied ranks
     */
    i = 0;
    while(i<=n-1)
    {
        j = i+1;
        while(j<=n-1)
        {
            if( ae_fp_neq(buf->ra1.ptr.p_double[j],buf->ra1.ptr.p_double[i]) )
            {
                break;
            }
            j = j+1;
        }
        for(k=i; k<=j-1; k++)
        {
            buf->ra1.ptr.p_double[k] = 1+(double)(i+j-1)/(double)2;
        }
        i = j;
    }
    
    /*
     * back to x
     */
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[buf->ia1.ptr.p_int[i]] = buf->ra1.ptr.p_double[i];
    }
}


/*$ End $*/
