
#ifndef _testalglibbasicsunit_h
#define _testalglibbasicsunit_h

#include "aenv.h"
#include "ialglib.h"
#include "alglibbasics.h"


/*$ Declarations $*/


typedef struct
{
    ae_bool bfield;
    double rfield;
    ae_int_t ifield;
    ae_complex cfield;
    ae_vector b1field;
    ae_vector r1field;
    ae_vector i1field;
    ae_vector c1field;
    ae_matrix b2field;
    ae_matrix r2field;
    ae_matrix i2field;
    ae_matrix c2field;
} rec1;


/*$ Body $*/


ae_bool testalglibbasics(ae_bool silent, ae_state *_state);
ae_bool _rec1_init(rec1* p, ae_state *_state, ae_bool make_automatic);
ae_bool _rec1_init_copy(rec1* dst, rec1* src, ae_state *_state, ae_bool make_automatic);
void _rec1_clear(rec1* p);


/*$ End $*/
#endif

