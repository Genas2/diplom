
#ifndef _testevdunit_h
#define _testevdunit_h

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
#include "evd.h"


/*$ Declarations $*/


/*$ Body $*/


/*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************/
ae_bool testevd(ae_bool silent, ae_state *_state);


/*$ End $*/
#endif

