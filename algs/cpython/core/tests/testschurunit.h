
#ifndef _testschurunit_h
#define _testschurunit_h

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
#include "hsschur.h"
#include "schur.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Testing Schur decomposition subroutine
*************************************************************************/
ae_bool testschur(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

