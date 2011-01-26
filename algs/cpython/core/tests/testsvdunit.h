
#ifndef _testsvdunit_h
#define _testsvdunit_h

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


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Testing SVD decomposition subroutine
*************************************************************************/
ae_bool testsvd(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

