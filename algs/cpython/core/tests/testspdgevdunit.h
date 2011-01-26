
#ifndef _testspdgevdunit_h
#define _testspdgevdunit_h

#include "aenv.h"
#include "ialglib.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "sblas.h"
#include "blas.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "apserv.h"
#include "matinv.h"
#include "hblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "spdgevd.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
ae_bool testspdgevd(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

