/*************************************************************************
Copyright 2008 by Sergey Bochkanov (ALGLIB project).

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
#include "tsort.h"


/*$ Declarations $*/
static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
This function sorts array of real keys by ascending.

Its results are:
* sorted array A
* permutation tables P1, P2

Algorithm outputs permutation tables using two formats:
* as usual permutation of [0..N-1]. If P1[i]=j, then sorted A[i]  contains
  value which was moved there from J-th position.
* as a sequence of pairwise permutations. Sorted A[] may  be  obtained  by
  swaping A[i] and A[P2[i]] for all i from 0 to N-1.
  
INPUT PARAMETERS:
    A       -   unsorted array
    N       -   array size

OUPUT PARAMETERS:
    A       -   sorted array
    P1, P2  -   permutation tables, array[N]
    
NOTES:
    this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.

  -- ALGLIB --
     Copyright 14.05.2008 by Bochkanov Sergey
*************************************************************************/
void tagsort(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector pv;
    ae_vector vp;
    ae_vector bufa;
    ae_vector bufb;
    ae_int_t lv;
    ae_int_t lp;
    ae_int_t rv;
    ae_int_t rp;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(p1);
    ae_vector_clear(p2);
    ae_vector_init(&pv, 0, DT_INT, _state, ae_true);
    ae_vector_init(&vp, 0, DT_INT, _state, ae_true);
    ae_vector_init(&bufa, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&bufb, 0, DT_INT, _state, ae_true);

    
    /*
     * Special cases
     */
    if( n<=0 )
    {
        ae_frame_leave(_state);
        return;
    }
    if( n==1 )
    {
        ae_vector_set_length(p1, 0+1, _state);
        ae_vector_set_length(p2, 0+1, _state);
        p1->ptr.p_int[0] = 0;
        p2->ptr.p_int[0] = 0;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case, N>1: prepare permutations table P1
     */
    ae_vector_set_length(p1, n-1+1, _state);
    for(i=0; i<=n-1; i++)
    {
        p1->ptr.p_int[i] = i;
    }
    
    /*
     * General case, N>1: sort, update P1
     */
    ae_vector_set_length(&bufa, n, _state);
    ae_vector_set_length(&bufb, n, _state);
    tagsortfasti(a, p1, &bufa, &bufb, n, _state);
    
    /*
     * General case, N>1: fill permutations table P2
     *
     * To fill P2 we maintain two arrays:
     * * PV, Position(Value). PV[i] contains position of I-th key at the moment
     * * VP, Value(Position). VP[i] contains key which has position I at the moment
     *
     * At each step we making permutation of two items:
     *   Left, which is given by position/value pair LP/LV
     *   and Right, which is given by RP/RV
     * and updating PV[] and VP[] correspondingly.
     */
    ae_vector_set_length(&pv, n-1+1, _state);
    ae_vector_set_length(&vp, n-1+1, _state);
    ae_vector_set_length(p2, n-1+1, _state);
    for(i=0; i<=n-1; i++)
    {
        pv.ptr.p_int[i] = i;
        vp.ptr.p_int[i] = i;
    }
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * calculate LP, LV, RP, RV
         */
        lp = i;
        lv = vp.ptr.p_int[lp];
        rv = p1->ptr.p_int[i];
        rp = pv.ptr.p_int[rv];
        
        /*
         * Fill P2
         */
        p2->ptr.p_int[i] = rp;
        
        /*
         * update PV and VP
         */
        vp.ptr.p_int[lp] = rv;
        vp.ptr.p_int[rp] = lv;
        pv.ptr.p_int[lv] = rp;
        pv.ptr.p_int[rv] = lp;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys and integer labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;
    ae_int_t tmpi;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpi = b->ptr.p_int[i];
            b->ptr.p_int[i] = b->ptr.p_int[j];
            b->ptr.p_int[j] = tmpi;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastirec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys and real labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfastr(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpr = b->ptr.p_double[i];
            b->ptr.p_double[i] = b->ptr.p_double[j];
            b->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastrrec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys without labels.

A is sorted, and that's all.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses buffer, BufA, which is N elements large. It may be
    preallocated (which will save some time) or not, in which case
    function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfast(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    tsort_tagsortfastrec(a, bufa, 0, n-1, _state);
}


/*************************************************************************
Heap operations: adds element to the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap (without new element).
                updated on output
    VA      -   value of the element being added
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappushi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     double va,
     ae_int_t vb,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t k;
    double v;


    if( *n<0 )
    {
        return;
    }
    
    /*
     * N=0 is a special case
     */
    if( *n==0 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        *n = *n+1;
        return;
    }
    
    /*
     * add current point to the heap
     * (add to the bottom, then move up)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = *n;
    *n = *n+1;
    while(j>0)
    {
        k = (j-1)/2;
        v = a->ptr.p_double[k];
        if( ae_fp_less(v,va) )
        {
            
            /*
             * swap with higher element
             */
            a->ptr.p_double[j] = v;
            b->ptr.p_int[j] = b->ptr.p_int[k];
            j = k;
        }
        else
        {
            
            /*
             * element in its place. terminate.
             */
            break;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: replaces top element with new element
(which is moved down)

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap
    VA      -   value of the element which replaces top element
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheapreplacetopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     double va,
     ae_int_t vb,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t k1;
    ae_int_t k2;
    double v;
    double v1;
    double v2;


    if( n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( n==1 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        return;
    }
    
    /*
     * move down through heap:
     * * J  -   current element
     * * K1 -   first child (always exists)
     * * K2 -   second child (may not exists)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = 0;
    k1 = 1;
    k2 = 2;
    while(k1<n)
    {
        if( k2>=n )
        {
            
            /*
             * only one child.
             *
             * swap and terminate (because this child
             * have no siblings due to heap structure)
             */
            v = a->ptr.p_double[k1];
            if( ae_fp_greater(v,va) )
            {
                a->ptr.p_double[j] = v;
                b->ptr.p_int[j] = b->ptr.p_int[k1];
                j = k1;
            }
            break;
        }
        else
        {
            
            /*
             * two childs
             */
            v1 = a->ptr.p_double[k1];
            v2 = a->ptr.p_double[k2];
            if( ae_fp_greater(v1,v2) )
            {
                if( ae_fp_less(va,v1) )
                {
                    a->ptr.p_double[j] = v1;
                    b->ptr.p_int[j] = b->ptr.p_int[k1];
                    j = k1;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if( ae_fp_less(va,v2) )
                {
                    a->ptr.p_double[j] = v2;
                    b->ptr.p_int[j] = b->ptr.p_int[k2];
                    j = k2;
                }
                else
                {
                    break;
                }
            }
            k1 = 2*j+1;
            k2 = 2*j+2;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: pops top element from the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap, N>=1

On output top element is moved to A[N-1], B[N-1], heap is reordered, N is
decreased by 1.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     ae_state *_state)
{
    double va;
    ae_int_t vb;


    if( *n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( *n==1 )
    {
        *n = 0;
        return;
    }
    
    /*
     * swap top element and last element,
     * then reorder heap
     */
    va = a->ptr.p_double[*n-1];
    vb = b->ptr.p_int[*n-1];
    a->ptr.p_double[*n-1] = a->ptr.p_double[0];
    b->ptr.p_int[*n-1] = b->ptr.p_int[0];
    *n = *n-1;
    tagheapreplacetopi(a, b, *n, va, vb, _state);
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediately if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpi = b->ptr.p_int[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_int[i+1] = b->ptr.p_int[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_int[k] = tmpi;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=2
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_int[k] = b->ptr.p_int[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_int[k] = b->ptr.p_int[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_int[k] = b->ptr.p_int[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastirec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastirec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastR: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    double tmpr2;
    ae_int_t tmpi;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpr2 = b->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_double[i+1] = b->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_double[k] = tmpr2;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_double[k] = b->ptr.p_double[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_double[k] = b->ptr.p_double[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_double[k] = b->ptr.p_double[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrec(a, bufa, i1, i1+cntless-1, _state);
    tsort_tagsortfastrec(a, bufa, i1+cntless+cnteq, i2, _state);
}


/*$ End $*/
