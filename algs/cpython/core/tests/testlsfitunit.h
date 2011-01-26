
#ifndef _testlsfitunit_h
#define _testlsfitunit_h

#include "aenv.h"
#include "ialglib.h"
#include "tsort.h"
#include "ratint.h"
#include "apserv.h"
#include "polint.h"
#include "spline1d.h"
#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "linmin.h"
#include "fbls.h"
#include "minlbfgs.h"
#include "minlm.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "lsfit.h"


/*$ Declarations $*/


/*$ Body $*/


ae_bool testlsfit(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

