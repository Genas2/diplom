
#ifndef _testdensesolverunit_h
#define _testdensesolverunit_h

#include "aenv.h"
#include "ialglib.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Test
*************************************************************************/
ae_bool testdensesolver(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

