
#ifndef _testidwintunit_h
#define _testidwintunit_h

#include "aenv.h"
#include "ialglib.h"
#include "tsort.h"
#include "apserv.h"
#include "nearestneighbor.h"
#include "reflections.h"
#include "hblas.h"
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
#include "idwint.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Testing IDW interpolation
*************************************************************************/
ae_bool testidwint(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

