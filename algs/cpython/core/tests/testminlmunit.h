
#ifndef _testminlmunit_h
#define _testminlmunit_h

#include "aenv.h"
#include "ialglib.h"
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
#include "apserv.h"
#include "matinv.h"
#include "linmin.h"
#include "fbls.h"
#include "minlbfgs.h"
#include "minlm.h"


/*$ Declarations $*/


/*$ Body $*/


ae_bool testminlm(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

