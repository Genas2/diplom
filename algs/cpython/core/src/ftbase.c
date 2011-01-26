/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
#include "ftbase.h"


/*$ Declarations $*/
static ae_int_t ftbase_ftbaseplanentrysize = 8;
static ae_int_t ftbase_ftbasecffttask = 0;
static ae_int_t ftbase_ftbaserfhttask = 1;
static ae_int_t ftbase_ftbaserffttask = 2;
static ae_int_t ftbase_fftcooleytukeyplan = 0;
static ae_int_t ftbase_fftbluesteinplan = 1;
static ae_int_t ftbase_fftcodeletplan = 2;
static ae_int_t ftbase_fhtcooleytukeyplan = 3;
static ae_int_t ftbase_fhtcodeletplan = 4;
static ae_int_t ftbase_fftrealcooleytukeyplan = 5;
static ae_int_t ftbase_fftemptyplan = 6;
static ae_int_t ftbase_fhtn2plan = 999;
static ae_int_t ftbase_ftbaseupdatetw = 4;
static ae_int_t ftbase_ftbasecodeletrecommended = 5;
static double ftbase_ftbaseinefficiencyfactor = 1.3;
static ae_int_t ftbase_ftbasemaxsmoothfactor = 5;
static void ftbase_ftbasegenerateplanrec(ae_int_t n,
     ae_int_t tasktype,
     ftplan* plan,
     ae_int_t* plansize,
     ae_int_t* precomputedsize,
     ae_int_t* planarraysize,
     ae_int_t* tmpmemsize,
     ae_int_t* stackmemsize,
     ae_int_t stackptr,
     ae_state *_state);
static void ftbase_ftbaseprecomputeplanrec(ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state);
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state);
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_internalreallintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state);
static void ftbase_fftarrayresize(/* Integer */ ae_vector* a,
     ae_int_t* asize,
     ae_int_t newasize,
     ae_state *_state);
static void ftbase_reffht(/* Real    */ ae_vector* a,
     ae_int_t n,
     ae_int_t offs,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
This subroutine generates FFT plan - a decomposition of a N-length FFT to
the more simpler operations. Plan consists of the root entry and the child
entries.

Subroutine parameters:
    N               task size
    
Output parameters:
    Plan            plan

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasegeneratecomplexfftplan(ae_int_t n,
     ftplan* plan,
     ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbasecffttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!", _state);
}


/*************************************************************************
Generates real FFT plan
*************************************************************************/
void ftbasegeneraterealfftplan(ae_int_t n, ftplan* plan, ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbaserffttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!", _state);
}


/*************************************************************************
Generates real FHT plan
*************************************************************************/
void ftbasegeneraterealfhtplan(ae_int_t n, ftplan* plan, ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbaserfhttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!", _state);
}


/*************************************************************************
This subroutine executes FFT/FHT plan.

If Plan is a:
* complex FFT plan  -   sizeof(A)=2*N,
                        A contains interleaved real/imaginary values
* real FFT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros
* real FHT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplan(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n,
     ftplan* plan,
     ae_state *_state)
{
    ae_int_t stackptr;


    stackptr = 0;
    ftbaseexecuteplanrec(a, aoffset, plan, 0, stackptr, _state);
}


/*************************************************************************
Recurrent subroutine for the FTBaseExecutePlan

Parameters:
    A           FFT'ed array
    AOffset     offset of the FFT'ed part (distance is measured in doubles)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplanrec(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t offs;
    ae_int_t offs1;
    ae_int_t offs2;
    ae_int_t offsa;
    ae_int_t offsb;
    ae_int_t offsp;
    double hk;
    double hnk;
    double x;
    double y;
    double bx;
    double by;
    ae_vector emptyarray;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double a2x;
    double a2y;
    double a3x;
    double a3y;
    double v0;
    double v1;
    double v2;
    double v3;
    double t1x;
    double t1y;
    double t2x;
    double t2y;
    double t3x;
    double t3y;
    double t4x;
    double t4y;
    double t5x;
    double t5y;
    double m1x;
    double m1y;
    double m2x;
    double m2y;
    double m3x;
    double m3y;
    double m4x;
    double m4y;
    double m5x;
    double m5y;
    double s1x;
    double s1y;
    double s2x;
    double s2y;
    double s3x;
    double s3y;
    double s4x;
    double s4y;
    double s5x;
    double s5y;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&emptyarray, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftemptyplan )
    {
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey plan
         * * transposition
         * * row-wise FFT
         * * twiddle factors:
         *   - TwBase is a basis twiddle factor for I=1, J=1
         *   - TwRow is a twiddle factor for a second element in a row (J=1)
         *   - Tw is a twiddle factor for a current element
         * * transposition again
         * * row-wise FFT again
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        ftbase_ffttwcalc(a, aoffset, n1, n2, _state);
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2*2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftrealcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey plan
         * * transposition
         * * row-wise FFT
         * * twiddle factors:
         *   - TwBase is a basis twiddle factor for I=1, J=1
         *   - TwRow is a twiddle factor for a second element in a row (J=1)
         *   - Tw is a twiddle factor for a current element
         * * transposition again
         * * row-wise FFT again
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1/2-1; i++)
        {
            
            /*
             * pack two adjacent smaller real FFT's together,
             * make one complex FFT,
             * unpack result
             */
            offs = aoffset+2*i*n2*2;
            for(k=0; k<=n2-1; k++)
            {
                a->ptr.p_double[offs+2*k+1] = a->ptr.p_double[offs+2*n2+2*k+0];
            }
            ftbaseexecuteplanrec(a, offs, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
            plan->tmpbuf.ptr.p_double[0] = a->ptr.p_double[offs+0];
            plan->tmpbuf.ptr.p_double[1] = 0;
            plan->tmpbuf.ptr.p_double[2*n2+0] = a->ptr.p_double[offs+1];
            plan->tmpbuf.ptr.p_double[2*n2+1] = 0;
            for(k=1; k<=n2-1; k++)
            {
                offs1 = 2*k;
                offs2 = 2*n2+2*k;
                hk = a->ptr.p_double[offs+2*k+0];
                hnk = a->ptr.p_double[offs+2*(n2-k)+0];
                plan->tmpbuf.ptr.p_double[offs1+0] = 0.5*(hk+hnk);
                plan->tmpbuf.ptr.p_double[offs2+1] = -0.5*(hk-hnk);
                hk = a->ptr.p_double[offs+2*k+1];
                hnk = a->ptr.p_double[offs+2*(n2-k)+1];
                plan->tmpbuf.ptr.p_double[offs2+0] = 0.5*(hk+hnk);
                plan->tmpbuf.ptr.p_double[offs1+1] = 0.5*(hk-hnk);
            }
            ae_v_move(&a->ptr.p_double[offs], 1, &plan->tmpbuf.ptr.p_double[0], 1, ae_v_len(offs,offs+2*n2*2-1));
        }
        if( n1%2!=0 )
        {
            ftbaseexecuteplanrec(a, aoffset+(n1-1)*n2*2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_ffttwcalc(a, aoffset, n2, n1, _state);
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey FHT plan:
         * * transpose                    \
         * * smaller FHT's                |
         * * pre-process                  |
         * * multiply by twiddle factors  | corresponds to multiplication by H1
         * * post-process                 |
         * * transpose again              /
         * * multiply by H2 (smaller FHT's)
         * * final transposition
         *
         * For more details see Vitezslav Vesely, "Fast algorithms
         * of Fourier and Hartley transform and their implementation in MATLAB",
         * page 31.
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        ftbase_internalreallintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        for(i=0; i<=n2-1; i++)
        {
            for(j=0; j<=n1-1; j++)
            {
                offsa = aoffset+i*n1;
                hk = a->ptr.p_double[offsa+j];
                hnk = a->ptr.p_double[offsa+(n1-j)%n1];
                offs = 2*(i*n1+j);
                plan->tmpbuf.ptr.p_double[offs+0] = -0.5*(hnk-hk);
                plan->tmpbuf.ptr.p_double[offs+1] = 0.5*(hk+hnk);
            }
        }
        ftbase_ffttwcalc(&plan->tmpbuf, 0, n1, n2, _state);
        for(j=0; j<=n1-1; j++)
        {
            a->ptr.p_double[aoffset+j] = plan->tmpbuf.ptr.p_double[2*j+0]+plan->tmpbuf.ptr.p_double[2*j+1];
        }
        if( n2%2==0 )
        {
            offs = 2*(n2/2)*n1;
            offsa = aoffset+n2/2*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+0]+plan->tmpbuf.ptr.p_double[offs+2*j+1];
            }
        }
        for(i=1; i<=(n2+1)/2-1; i++)
        {
            offs = 2*i*n1;
            offs2 = 2*(n2-i)*n1;
            offsa = aoffset+i*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+1]+plan->tmpbuf.ptr.p_double[offs2+2*j+0];
            }
            offsa = aoffset+(n2-i)*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+0]+plan->tmpbuf.ptr.p_double[offs2+2*j+1];
            }
        }
        ftbase_internalreallintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_internalreallintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtn2plan )
    {
        
        /*
         * Cooley-Tukey FHT plan
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        ftbase_reffht(a, n, aoffset, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==2 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            v0 = a0x+a1x;
            v1 = a0y+a1y;
            v2 = a0x-a1x;
            v3 = a0y-a1y;
            a->ptr.p_double[aoffset+0] = v0;
            a->ptr.p_double[aoffset+1] = v1;
            a->ptr.p_double[aoffset+2] = v2;
            a->ptr.p_double[aoffset+3] = v3;
            ae_frame_leave(_state);
            return;
        }
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            t1x = a1x+a2x;
            t1y = a1y+a2y;
            a0x = a0x+t1x;
            a0y = a0y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a1y-a2y);
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            s1y = a0y+m1y;
            a1x = s1x+m2x;
            a1y = s1y+m2y;
            a2x = s1x-m2x;
            a2y = s1y-m2y;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = a0y;
            a->ptr.p_double[aoffset+2] = a1x;
            a->ptr.p_double[aoffset+3] = a1y;
            a->ptr.p_double[aoffset+4] = a2x;
            a->ptr.p_double[aoffset+5] = a2y;
            ae_frame_leave(_state);
            return;
        }
        if( n==4 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            a3x = a->ptr.p_double[aoffset+6];
            a3y = a->ptr.p_double[aoffset+7];
            t1x = a0x+a2x;
            t1y = a0y+a2y;
            t2x = a1x+a3x;
            t2y = a1y+a3y;
            m2x = a0x-a2x;
            m2y = a0y-a2y;
            m3x = a1y-a3y;
            m3y = a3x-a1x;
            a->ptr.p_double[aoffset+0] = t1x+t2x;
            a->ptr.p_double[aoffset+1] = t1y+t2y;
            a->ptr.p_double[aoffset+4] = t1x-t2x;
            a->ptr.p_double[aoffset+5] = t1y-t2y;
            a->ptr.p_double[aoffset+2] = m2x+m3x;
            a->ptr.p_double[aoffset+3] = m2y+m3y;
            a->ptr.p_double[aoffset+6] = m2x-m3x;
            a->ptr.p_double[aoffset+7] = m2y-m3y;
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            c3 = plan->precomputed.ptr.p_double[offs+2];
            c4 = plan->precomputed.ptr.p_double[offs+3];
            c5 = plan->precomputed.ptr.p_double[offs+4];
            t1x = a->ptr.p_double[aoffset+2]+a->ptr.p_double[aoffset+8];
            t1y = a->ptr.p_double[aoffset+3]+a->ptr.p_double[aoffset+9];
            t2x = a->ptr.p_double[aoffset+4]+a->ptr.p_double[aoffset+6];
            t2y = a->ptr.p_double[aoffset+5]+a->ptr.p_double[aoffset+7];
            t3x = a->ptr.p_double[aoffset+2]-a->ptr.p_double[aoffset+8];
            t3y = a->ptr.p_double[aoffset+3]-a->ptr.p_double[aoffset+9];
            t4x = a->ptr.p_double[aoffset+6]-a->ptr.p_double[aoffset+4];
            t4y = a->ptr.p_double[aoffset+7]-a->ptr.p_double[aoffset+5];
            t5x = t1x+t2x;
            t5y = t1y+t2y;
            a->ptr.p_double[aoffset+0] = a->ptr.p_double[aoffset+0]+t5x;
            a->ptr.p_double[aoffset+1] = a->ptr.p_double[aoffset+1]+t5y;
            m1x = c1*t5x;
            m1y = c1*t5y;
            m2x = c2*(t1x-t2x);
            m2y = c2*(t1y-t2y);
            m3x = -c3*(t3y+t4y);
            m3y = c3*(t3x+t4x);
            m4x = -c4*t4y;
            m4y = c4*t4x;
            m5x = -c5*t3y;
            m5y = c5*t3x;
            s3x = m3x-m4x;
            s3y = m3y-m4y;
            s5x = m3x+m5x;
            s5y = m3y+m5y;
            s1x = a->ptr.p_double[aoffset+0]+m1x;
            s1y = a->ptr.p_double[aoffset+1]+m1y;
            s2x = s1x+m2x;
            s2y = s1y+m2y;
            s4x = s1x-m2x;
            s4y = s1y-m2y;
            a->ptr.p_double[aoffset+2] = s2x+s3x;
            a->ptr.p_double[aoffset+3] = s2y+s3y;
            a->ptr.p_double[aoffset+4] = s4x+s5x;
            a->ptr.p_double[aoffset+5] = s4y+s5y;
            a->ptr.p_double[aoffset+6] = s4x-s5x;
            a->ptr.p_double[aoffset+7] = s4y-s5y;
            a->ptr.p_double[aoffset+8] = s2x-s3x;
            a->ptr.p_double[aoffset+9] = s2y-s3y;
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==2 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a->ptr.p_double[aoffset+0] = a0x+a1x;
            a->ptr.p_double[aoffset+1] = a0x-a1x;
            ae_frame_leave(_state);
            return;
        }
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a2x = a->ptr.p_double[aoffset+2];
            t1x = a1x+a2x;
            a0x = a0x+t1x;
            m1x = c1*t1x;
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = s1x-m2y;
            a->ptr.p_double[aoffset+2] = s1x+m2y;
            ae_frame_leave(_state);
            return;
        }
        if( n==4 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a2x = a->ptr.p_double[aoffset+2];
            a3x = a->ptr.p_double[aoffset+3];
            t1x = a0x+a2x;
            t2x = a1x+a3x;
            m2x = a0x-a2x;
            m3y = a3x-a1x;
            a->ptr.p_double[aoffset+0] = t1x+t2x;
            a->ptr.p_double[aoffset+1] = m2x-m3y;
            a->ptr.p_double[aoffset+2] = t1x-t2x;
            a->ptr.p_double[aoffset+3] = m2x+m3y;
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            c3 = plan->precomputed.ptr.p_double[offs+2];
            c4 = plan->precomputed.ptr.p_double[offs+3];
            c5 = plan->precomputed.ptr.p_double[offs+4];
            t1x = a->ptr.p_double[aoffset+1]+a->ptr.p_double[aoffset+4];
            t2x = a->ptr.p_double[aoffset+2]+a->ptr.p_double[aoffset+3];
            t3x = a->ptr.p_double[aoffset+1]-a->ptr.p_double[aoffset+4];
            t4x = a->ptr.p_double[aoffset+3]-a->ptr.p_double[aoffset+2];
            t5x = t1x+t2x;
            v0 = a->ptr.p_double[aoffset+0]+t5x;
            a->ptr.p_double[aoffset+0] = v0;
            m2x = c2*(t1x-t2x);
            m3y = c3*(t3x+t4x);
            s3y = m3y-c4*t4x;
            s5y = m3y+c5*t3x;
            s1x = v0+c1*t5x;
            s2x = s1x+m2x;
            s4x = s1x-m2x;
            a->ptr.p_double[aoffset+1] = s2x-s3y;
            a->ptr.p_double[aoffset+2] = s4x-s5y;
            a->ptr.p_double[aoffset+3] = s4x+s5y;
            a->ptr.p_double[aoffset+4] = s2x+s3y;
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftbluesteinplan )
    {
        
        /*
         * Bluestein plan:
         * 1. multiply by precomputed coefficients
         * 2. make convolution: forward FFT, multiplication by precomputed FFT
         *    and backward FFT. backward FFT is represented as
         *
         *        invfft(x) = fft(x')'/M
         *
         *    for performance reasons reduction of inverse FFT to
         *    forward FFT is merged with multiplication of FFT components
         *    and last stage of Bluestein's transformation.
         * 3. post-multiplication by Bluestein factors
         */
        n = plan->plan.ptr.p_int[entryoffset+1];
        m = plan->plan.ptr.p_int[entryoffset+4];
        offs = plan->plan.ptr.p_int[entryoffset+7];
        for(i=stackptr+2*n; i<=stackptr+2*m-1; i++)
        {
            plan->stackbuf.ptr.p_double[i] = 0;
        }
        offsp = offs+2*m;
        offsa = aoffset;
        offsb = stackptr;
        for(i=0; i<=n-1; i++)
        {
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            x = a->ptr.p_double[offsa+0];
            y = a->ptr.p_double[offsa+1];
            plan->stackbuf.ptr.p_double[offsb+0] = x*bx-y*(-by);
            plan->stackbuf.ptr.p_double[offsb+1] = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        ftbaseexecuteplanrec(&plan->stackbuf, stackptr, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr+2*2*m, _state);
        offsb = stackptr;
        offsp = offs;
        for(i=0; i<=m-1; i++)
        {
            x = plan->stackbuf.ptr.p_double[offsb+0];
            y = plan->stackbuf.ptr.p_double[offsb+1];
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            plan->stackbuf.ptr.p_double[offsb+0] = x*bx-y*by;
            plan->stackbuf.ptr.p_double[offsb+1] = -(x*by+y*bx);
            offsb = offsb+2;
            offsp = offsp+2;
        }
        ftbaseexecuteplanrec(&plan->stackbuf, stackptr, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr+2*2*m, _state);
        offsb = stackptr;
        offsp = offs+2*m;
        offsa = aoffset;
        for(i=0; i<=n-1; i++)
        {
            x = plan->stackbuf.ptr.p_double[offsb+0]/m;
            y = -plan->stackbuf.ptr.p_double[offsb+1]/m;
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            a->ptr.p_double[offsa+0] = x*bx-y*(-by);
            a->ptr.p_double[offsa+1] = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Returns good factorization N=N1*N2.

Usually N1<=N2 (but not always - small N's may be exception).
if N1<>1 then N2<>1.

Factorization is chosen depending on task type and codelets we have.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasefactorize(ae_int_t n,
     ae_int_t tasktype,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *_state)
{
    ae_int_t j;

    *n1 = 0;
    *n2 = 0;

    *n1 = 0;
    *n2 = 0;
    
    /*
     * try to find good codelet
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended; j>=2; j--)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * try to factorize N
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended+1; j<=n-1; j++)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * looks like N is prime :(
     */
    if( *n1*(*n2)!=n )
    {
        *n1 = 1;
        *n2 = n;
    }
    
    /*
     * normalize
     */
    if( *n2==1&&*n1!=1 )
    {
        *n2 = *n1;
        *n1 = 1;
    }
}


/*************************************************************************
Is number smooth?

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool ftbaseissmooth(ae_int_t n, ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    for(i=2; i<=ftbase_ftbasemaxsmoothfactor; i++)
    {
        while(n%i==0)
        {
            n = n/i;
        }
    }
    result = n==1;
    return result;
}


/*************************************************************************
Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 1, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
greater than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 2, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns estimate of FLOP count for the FFT.

It is only an estimate based on operations count for the PERFECT FFT
and relative inefficiency of the algorithm actually used.

N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
double ftbasegetflopestimate(ae_int_t n, ae_state *_state)
{
    double result;


    result = ftbase_ftbaseinefficiencyfactor*(4*n*ae_log(n, _state)/ae_log(2, _state)-6*n+8);
    return result;
}


/*************************************************************************
Recurrent subroutine for the FFTGeneratePlan:

PARAMETERS:
    N                   plan size
    IsReal              whether input is real or not.
                        subroutine MUST NOT ignore this flag because real
                        inputs comes with non-initialized imaginary parts,
                        so ignoring this flag will result in corrupted output
    HalfOut             whether full output or only half of it from 0 to
                        floor(N/2) is needed. This flag may be ignored if
                        doing so will simplify calculations
    Plan                plan array
    PlanSize            size of used part (in integers)
    PrecomputedSize     size of precomputed array allocated yet
    PlanArraySize       plan array size (actual)
    TmpMemSize          temporary memory required size
    BluesteinMemSize    temporary memory required size

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbasegenerateplanrec(ae_int_t n,
     ae_int_t tasktype,
     ftplan* plan,
     ae_int_t* plansize,
     ae_int_t* precomputedsize,
     ae_int_t* planarraysize,
     ae_int_t* tmpmemsize,
     ae_int_t* stackmemsize,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_int_t k;
    ae_int_t m;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t esize;
    ae_int_t entryoffset;


    
    /*
     * prepare
     */
    if( *plansize+ftbase_ftbaseplanentrysize>(*planarraysize) )
    {
        ftbase_fftarrayresize(&plan->plan, planarraysize, 8*(*planarraysize), _state);
    }
    entryoffset = *plansize;
    esize = ftbase_ftbaseplanentrysize;
    *plansize = *plansize+esize;
    
    /*
     * if N=1, generate empty plan and exit
     */
    if( n==1 )
    {
        plan->plan.ptr.p_int[entryoffset+0] = esize;
        plan->plan.ptr.p_int[entryoffset+1] = -1;
        plan->plan.ptr.p_int[entryoffset+2] = -1;
        plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftemptyplan;
        plan->plan.ptr.p_int[entryoffset+4] = -1;
        plan->plan.ptr.p_int[entryoffset+5] = -1;
        plan->plan.ptr.p_int[entryoffset+6] = -1;
        plan->plan.ptr.p_int[entryoffset+7] = -1;
        return;
    }
    
    /*
     * generate plans
     */
    ftbasefactorize(n, tasktype, &n1, &n2, _state);
    if( tasktype==ftbase_ftbasecffttask||tasktype==ftbase_ftbaserffttask )
    {
        
        /*
         * complex FFT plans
         */
        if( n1!=1 )
        {
            
            /*
             * Cooley-Tukey plan (real or complex)
             *
             * Note that child plans are COMPLEX
             * (whether plan itself is complex or not).
             */
            *tmpmemsize = ae_maxint(*tmpmemsize, 2*n1*n2, _state);
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            if( tasktype==ftbase_ftbasecffttask )
            {
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftcooleytukeyplan;
            }
            else
            {
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftrealcooleytukeyplan;
            }
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = *plansize;
            ftbase_ftbasegenerateplanrec(n1, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+6] = *plansize;
            ftbase_ftbasegenerateplanrec(n2, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            return;
        }
        else
        {
            if( ((n==2||n==3)||n==4)||n==5 )
            {
                
                /*
                 * hard-coded plan
                 */
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n1;
                plan->plan.ptr.p_int[entryoffset+2] = n2;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftcodeletplan;
                plan->plan.ptr.p_int[entryoffset+4] = 0;
                plan->plan.ptr.p_int[entryoffset+5] = -1;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                if( n==3 )
                {
                    *precomputedsize = *precomputedsize+2;
                }
                if( n==5 )
                {
                    *precomputedsize = *precomputedsize+5;
                }
                return;
            }
            else
            {
                
                /*
                 * Bluestein's plan
                 *
                 * Select such M that M>=2*N-1, M is composite, and M's
                 * factors are 2, 3, 5
                 */
                k = 2*n2-1;
                m = ftbasefindsmooth(k, _state);
                *tmpmemsize = ae_maxint(*tmpmemsize, 2*m, _state);
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n2;
                plan->plan.ptr.p_int[entryoffset+2] = -1;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftbluesteinplan;
                plan->plan.ptr.p_int[entryoffset+4] = m;
                plan->plan.ptr.p_int[entryoffset+5] = *plansize;
                stackptr = stackptr+2*2*m;
                *stackmemsize = ae_maxint(*stackmemsize, stackptr, _state);
                ftbase_ftbasegenerateplanrec(m, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
                stackptr = stackptr-2*2*m;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                *precomputedsize = *precomputedsize+2*m+2*n;
                return;
            }
        }
    }
    if( tasktype==ftbase_ftbaserfhttask )
    {
        
        /*
         * real FHT plans
         */
        if( n1!=1 )
        {
            
            /*
             * Cooley-Tukey plan
             *
             */
            *tmpmemsize = ae_maxint(*tmpmemsize, 2*n1*n2, _state);
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtcooleytukeyplan;
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = *plansize;
            ftbase_ftbasegenerateplanrec(n1, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+6] = *plansize;
            ftbase_ftbasegenerateplanrec(n2, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            return;
        }
        else
        {
            
            /*
             * N2 plan
             */
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtn2plan;
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = -1;
            plan->plan.ptr.p_int[entryoffset+6] = -1;
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            if( ((n==2||n==3)||n==4)||n==5 )
            {
                
                /*
                 * hard-coded plan
                 */
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n1;
                plan->plan.ptr.p_int[entryoffset+2] = n2;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtcodeletplan;
                plan->plan.ptr.p_int[entryoffset+4] = 0;
                plan->plan.ptr.p_int[entryoffset+5] = -1;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                if( n==3 )
                {
                    *precomputedsize = *precomputedsize+2;
                }
                if( n==5 )
                {
                    *precomputedsize = *precomputedsize+5;
                }
                return;
            }
            return;
        }
    }
}


/*************************************************************************
Recurrent subroutine for precomputing FFT plans

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbaseprecomputeplanrec(ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t offs;
    double v;
    ae_vector emptyarray;
    double bx;
    double by;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&emptyarray, 0, DT_REAL, _state, ae_true);

    if( (plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcooleytukeyplan||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftrealcooleytukeyplan)||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcooleytukeyplan )
    {
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcodeletplan||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            plan->precomputed.ptr.p_double[offs+0] = ae_cos(2*ae_pi/3, _state)-1;
            plan->precomputed.ptr.p_double[offs+1] = ae_sin(2*ae_pi/3, _state);
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            v = 2*ae_pi/5;
            plan->precomputed.ptr.p_double[offs+0] = (ae_cos(v, _state)+ae_cos(2*v, _state))/2-1;
            plan->precomputed.ptr.p_double[offs+1] = (ae_cos(v, _state)-ae_cos(2*v, _state))/2;
            plan->precomputed.ptr.p_double[offs+2] = -ae_sin(v, _state);
            plan->precomputed.ptr.p_double[offs+3] = -(ae_sin(v, _state)+ae_sin(2*v, _state));
            plan->precomputed.ptr.p_double[offs+4] = ae_sin(v, _state)-ae_sin(2*v, _state);
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftbluesteinplan )
    {
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        n = plan->plan.ptr.p_int[entryoffset+1];
        m = plan->plan.ptr.p_int[entryoffset+4];
        offs = plan->plan.ptr.p_int[entryoffset+7];
        for(i=0; i<=2*m-1; i++)
        {
            plan->precomputed.ptr.p_double[offs+i] = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            bx = ae_cos(ae_pi*ae_sqr(i, _state)/n, _state);
            by = ae_sin(ae_pi*ae_sqr(i, _state)/n, _state);
            plan->precomputed.ptr.p_double[offs+2*i+0] = bx;
            plan->precomputed.ptr.p_double[offs+2*i+1] = by;
            plan->precomputed.ptr.p_double[offs+2*m+2*i+0] = bx;
            plan->precomputed.ptr.p_double[offs+2*m+2*i+1] = by;
            if( i>0 )
            {
                plan->precomputed.ptr.p_double[offs+2*(m-i)+0] = bx;
                plan->precomputed.ptr.p_double[offs+2*(m-i)+1] = by;
            }
        }
        ftbaseexecuteplanrec(&plan->precomputed, offs, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Twiddle factors calculation

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    ae_int_t idx;
    ae_int_t offs;
    double x;
    double y;
    double twxm1;
    double twy;
    double twbasexm1;
    double twbasey;
    double twrowxm1;
    double twrowy;
    double tmpx;
    double tmpy;
    double v;


    n = n1*n2;
    v = -2*ae_pi/n;
    twbasexm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
    twbasey = ae_sin(v, _state);
    twrowxm1 = 0;
    twrowy = 0;
    for(i=0; i<=n2-1; i++)
    {
        twxm1 = 0;
        twy = 0;
        for(j=0; j<=n1-1; j++)
        {
            idx = i*n1+j;
            offs = aoffset+2*idx;
            x = a->ptr.p_double[offs+0];
            y = a->ptr.p_double[offs+1];
            tmpx = x*twxm1-y*twy;
            tmpy = x*twy+y*twxm1;
            a->ptr.p_double[offs+0] = x+tmpx;
            a->ptr.p_double[offs+1] = y+tmpy;
            
            /*
             * update Tw: Tw(new) = Tw(old)*TwRow
             */
            if( j<n1-1 )
            {
                if( j%ftbase_ftbaseupdatetw==0 )
                {
                    v = -2*ae_pi*i*(j+1)/n;
                    twxm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
                    twy = ae_sin(v, _state);
                }
                else
                {
                    tmpx = twrowxm1+twxm1*twrowxm1-twy*twrowy;
                    tmpy = twrowy+twxm1*twrowy+twy*twrowxm1;
                    twxm1 = twxm1+tmpx;
                    twy = twy+tmpy;
                }
            }
        }
        
        /*
         * update TwRow: TwRow(new) = TwRow(old)*TwBase
         */
        if( i<n2-1 )
        {
            if( j%ftbase_ftbaseupdatetw==0 )
            {
                v = -2*ae_pi*(i+1)/n;
                twrowxm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
                twrowy = ae_sin(v, _state);
            }
            else
            {
                tmpx = twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
                tmpy = twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
                twrowxm1 = twrowxm1+tmpx;
                twrowy = twrowy+tmpy;
            }
        }
    }
}


/*************************************************************************
Linear transpose: transpose complex matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{


    ftbase_ffticltrec(a, astart, n, buf, 0, m, m, n, _state);
    ae_v_move(&a->ptr.p_double[astart], 1, &buf->ptr.p_double[0], 1, ae_v_len(astart,astart+2*m*n-1));
}


/*************************************************************************
Linear transpose: transpose real matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_internalreallintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{


    ftbase_fftirltrec(a, astart, n, buf, 0, m, m, n, _state);
    ae_v_move(&a->ptr.p_double[astart], 1, &buf->ptr.p_double[0], 1, ae_v_len(astart,astart+m*n-1));
}


/*************************************************************************
Recurrent subroutine for a InternalComplexLinTranspose

Write A^T to B, where:
* A is m*n complex matrix stored in array A as pairs of real/image values,
  beginning from AStart position, with AStride stride
* B is n*m complex matrix stored in array B as pairs of real/image values,
  beginning from BStart position, with BStride stride
stride is measured in complex numbers, i.e. in real/image pairs.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        m2 = 2*bstride;
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+2*i;
            idx2 = astart+2*i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1+0] = a->ptr.p_double[idx2+0];
                b->ptr.p_double[idx1+1] = a->ptr.p_double[idx2+1];
                idx1 = idx1+m2;
                idx2 = idx2+2;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_ffticltrec(a, astart+2*n1, astride, b, bstart+2*n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_ffticltrec(a, astart+2*m1*astride, astride, b, bstart+2*m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
Recurrent subroutine for a InternalRealLinTranspose


  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+i;
            idx2 = astart+i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1] = a->ptr.p_double[idx2];
                idx1 = idx1+bstride;
                idx2 = idx2+1;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_fftirltrec(a, astart+n1, astride, b, bstart+n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_fftirltrec(a, astart+m1*astride, astride, b, bstart+m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
recurrent subroutine for FFTFindSmoothRec

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state)
{


    ae_assert(ftbase_ftbasemaxsmoothfactor<=5, "FTBaseFindSmoothRec: internal error!", _state);
    if( seed>=n )
    {
        *best = ae_minint(*best, seed, _state);
        return;
    }
    if( leastfactor<=2 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*2, 2, best, _state);
    }
    if( leastfactor<=3 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*3, 3, best, _state);
    }
    if( leastfactor<=5 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*5, 5, best, _state);
    }
}


/*************************************************************************
Internal subroutine: array resize

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_fftarrayresize(/* Integer */ ae_vector* a,
     ae_int_t* asize,
     ae_int_t newasize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_INT, _state, ae_true);

    ae_vector_set_length(&tmp, *asize, _state);
    for(i=0; i<=*asize-1; i++)
    {
        tmp.ptr.p_int[i] = a->ptr.p_int[i];
    }
    ae_vector_set_length(a, newasize, _state);
    for(i=0; i<=*asize-1; i++)
    {
        a->ptr.p_int[i] = tmp.ptr.p_int[i];
    }
    *asize = newasize;
    ae_frame_leave(_state);
}


/*************************************************************************
Reference FHT stub
*************************************************************************/
static void ftbase_reffht(/* Real    */ ae_vector* a,
     ae_int_t n,
     ae_int_t offs,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i;
    ae_int_t j;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);

    ae_assert(n>0, "RefFHTR1D: incorrect N!", _state);
    ae_vector_set_length(&buf, n, _state);
    for(i=0; i<=n-1; i++)
    {
        v = 0;
        for(j=0; j<=n-1; j++)
        {
            v = v+a->ptr.p_double[offs+j]*(ae_cos(2*ae_pi*i*j/n, _state)+ae_sin(2*ae_pi*i*j/n, _state));
        }
        buf.ptr.p_double[i] = v;
    }
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_double[offs+i] = buf.ptr.p_double[i];
    }
    ae_frame_leave(_state);
}


ae_bool _ftplan_init(ftplan* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->plan, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->precomputed, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->stackbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _ftplan_init_copy(ftplan* dst, ftplan* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->plan, &src->plan, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->precomputed, &src->precomputed, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpbuf, &src->tmpbuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->stackbuf, &src->stackbuf, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _ftplan_clear(ftplan* p)
{
    ae_vector_clear(&p->plan);
    ae_vector_clear(&p->precomputed);
    ae_vector_clear(&p->tmpbuf);
    ae_vector_clear(&p->stackbuf);
}


/*$ End $*/
