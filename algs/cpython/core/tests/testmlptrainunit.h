
#ifndef _testmlptrainunit_h
#define _testmlptrainunit_h

#include "aenv.h"
#include "ialglib.h"
#include "mlpbase.h"
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
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "mlptrain.h"


/*$ Declarations $*/


/*$ Body $*/


ae_bool testmlptrain(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

