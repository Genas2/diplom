/*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _fbls_h
#define _fbls_h

#include "aenv.h"
#include "ialglib.h"
#include "ablasf.h"
#include "ablas.h"


/*$ Declarations $*/


/*************************************************************************
Structure which stores state of linear CG solver between subsequent calls
of FBLSCgIteration(). Initialized with FBLSCGCreate().

USAGE:
1. call to FBLSCGCreate()
2. F:=FBLSCgIteration(State)
3. if F is False, iterations are over
4. otherwise, fill State.AX with A*x, State.XAX with x'*A*x
5. goto 2

If you want to rerminate iterations, pass zero or negative value to XAX.

FIELDS:
    E1      -   2-norm of residual at the start
    E2      -   2-norm of residual at the end
    X       -   on return from FBLSCgIteration() it contains vector for
                matrix-vector product
    AX      -   must be filled with A*x if FBLSCgIteration() returned True
    XAX     -   must be filled with x'*A*x
    XK      -   contains result (if FBLSCgIteration() returned False)
    
Other fields are private and should not be used by outsiders.
*************************************************************************/
typedef struct
{
    double e1;
    double e2;
    ae_vector x;
    ae_vector ax;
    double xax;
    ae_int_t n;
    ae_vector rk;
    ae_vector rk1;
    ae_vector xk;
    ae_vector xk1;
    ae_vector pk;
    ae_vector pk1;
    ae_vector b;
    rcommstate rstate;
    ae_vector tmp2;
} fblslincgstate;


/*$ Body $*/


/*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

INPUT PARAMETERS:
    CHA     -   Cholesky decomposition of A
    SqrtScaleA- square root of scale factor ScaleA
    N       -   matrix size
    IsUpper -   storage type
    XB      -   right part
    Tmp     -   buffer; function automatically allocates it, if it is  too
                small.  It  can  be  reused  if function is called several
                times.
                
OUTPUT PARAMETERS:
    XB      -   solution

NOTES: no assertion or tests are done during algorithm operation

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void fblscholeskysolve(/* Real    */ ae_matrix* cha,
     double sqrtscalea,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* xb,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);


/*************************************************************************
Fast basic linear solver: linear SPD CG

Solves (A^T*A + alpha*I)*x = b where:
* A is MxN matrix
* alpha>0 is a scalar
* I is NxN identity matrix
* b is Nx1 vector
* X is Nx1 unknown vector.

N iterations of linear conjugate gradient are used to solve problem.

INPUT PARAMETERS:
    A   -   array[M,N], matrix
    M   -   number of rows
    N   -   number of unknowns
    B   -   array[N], right part
    X   -   initial approxumation, array[N]
    Buf -   buffer; function automatically allocates it, if it is too
            small. It can be reused if function is called several times
            with same M and N.
            
OUTPUT PARAMETERS:
    X   -   improved solution
    
NOTES:
*   solver checks quality of improved solution. If (because of problem
    condition number, numerical noise, etc.) new solution is WORSE than
    original approximation, then original approximation is returned.
*   solver assumes that both A, B, Alpha are well scaled (i.e. they are
    less than sqrt(overflow) and greater than sqrt(underflow)).
    
  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void fblssolvecgx(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     double alpha,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* buf,
     ae_state *_state);


/*************************************************************************
Construction of linear conjugate gradient solver.

State parameter passed using "var" semantics (i.e. previous state  is  NOT
erased). When it is already initialized, we can reause prevously allocated
memory.

INPUT PARAMETERS:
    X       -   initial solution
    B       -   right part
    N       -   system size
    State   -   structure; may be preallocated, if we want to reuse memory

OUTPUT PARAMETERS:
    State   -   structure which is used by FBLSCGIteration() to store
                algorithm state between subsequent calls.

NOTE: no error checking is done; caller must check all parameters, prevent
      overflows, and so on.

  -- ALGLIB --
     Copyright 22.10.2009 by Bochkanov Sergey
*************************************************************************/
void fblscgcreate(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     fblslincgstate* state,
     ae_state *_state);


/*************************************************************************
Linear CG solver, function relying on reverse communication to calculate
matrix-vector products.

See comments for FBLSLinCGState structure for more info.

  -- ALGLIB --
     Copyright 22.10.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool fblscgiteration(fblslincgstate* state, ae_state *_state);
ae_bool _fblslincgstate_init(fblslincgstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _fblslincgstate_init_copy(fblslincgstate* dst, fblslincgstate* src, ae_state *_state, ae_bool make_automatic);
void _fblslincgstate_clear(fblslincgstate* p);


/*$ End $*/
#endif

