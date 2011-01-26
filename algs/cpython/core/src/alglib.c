/*

*/

#include <stdlib.h>
#include "aenv.h"
#include "alglib.h"

#ifdef X_FOR_WINDOWS
#define DLLEXPORT __declspec(dllexport)
#endif
#ifdef X_FOR_LINUX
#define DLLEXPORT
#endif

#pragma pack(1)
#define x_nb 16

enum
{
    X_OK = 0,
    X_MALLOC_ERROR = 1,
    X_DIV_BY_ZERO_ERROR = 2,
    X_32_64_ERROR = 3,
    X_ARRAY_TOO_LARGE = 4,
    X_ASSERTION_FAILED = 5
};

DLLEXPORT int x_malloc(void **p, ae_int64_t size)
{
    size_t volatile tmp;
    tmp = (size_t)size;
    if( tmp!=size )
        return X_MALLOC_ERROR;
    *p = aligned_malloc(tmp,16);
    if( *p || tmp==0 )
        return X_OK;
    return X_MALLOC_ERROR;
}
DLLEXPORT int x_free(void *p)
{
    aligned_free(p);
    return X_OK;
}
DLLEXPORT ae_int64_t x_alloc_counter()
{
    return _alloc_counter;
}
DLLEXPORT ae_bool x_is_symmetric_e_(x_matrix *a)
{
    return x_is_symmetric(a);
}
DLLEXPORT ae_bool x_is_hermitian_e_(x_matrix *a)
{
    return x_is_hermitian(a);
}
DLLEXPORT ae_bool x_force_symmetric_e_(x_matrix *a)
{
    return x_force_symmetric(a);
}
DLLEXPORT ae_bool x_force_hermitian_e_(x_matrix *a)
{
    return x_force_hermitian(a);
}




typedef struct
{
    hqrndstate obj;
} x_hqrndstate;
x_hqrndstate* x_obj_alloc_hqrndstate(ae_state *_state)
{
    x_hqrndstate *result;
    result = ae_malloc(sizeof(x_hqrndstate), _state);
    _hqrndstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_hqrndstate(x_hqrndstate *obj)
{
    if( obj==NULL )
        return;
    _hqrndstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_hqrndrandomize(const char **errormsg, x_hqrndstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *state = x_obj_alloc_hqrndstate(&_alglib_env_state);
    hqrndrandomize(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrndseed(const char **errormsg, ae_int_t* s1, ae_int_t* s2, x_hqrndstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *state = x_obj_alloc_hqrndstate(&_alglib_env_state);
    hqrndseed(*s1, *s2, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrnduniformr(const char **errormsg, double* result, x_hqrndstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = hqrnduniformr(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrnduniformi(const char **errormsg, ae_int_t* result, x_hqrndstate** state, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = hqrnduniformi(&(*state)->obj, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrndnormal(const char **errormsg, double* result, x_hqrndstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = hqrndnormal(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrndunit2(const char **errormsg, x_hqrndstate** state, double* x, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    hqrndunit2(&(*state)->obj, x, y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrndnormal2(const char **errormsg, x_hqrndstate** state, double* x1, double* x2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    hqrndnormal2(&(*state)->obj, x1, x2, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hqrndexponential(const char **errormsg, double* result, x_hqrndstate** state, double* lambdav)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = hqrndexponential(&(*state)->obj, *lambdav, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    kdtree obj;
} x_kdtree;
x_kdtree* x_obj_alloc_kdtree(ae_state *_state)
{
    x_kdtree *result;
    result = ae_malloc(sizeof(x_kdtree), _state);
    _kdtree_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_kdtree(x_kdtree *obj)
{
    if( obj==NULL )
        return;
    _kdtree_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_kdtreebuild(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* nx, ae_int_t* ny, ae_int_t* normtype, x_kdtree** kdt)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *kdt = x_obj_alloc_kdtree(&_alglib_env_state);
    kdtreebuild(&_xy, *n, *nx, *ny, *normtype, &(*kdt)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreebuildtagged(const char **errormsg, x_matrix* xy, x_vector* tags, ae_int_t* n, ae_int_t* nx, ae_int_t* ny, ae_int_t* normtype, x_kdtree** kdt)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _tags;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tags, tags, &_alglib_env_state, ae_true);
    *kdt = x_obj_alloc_kdtree(&_alglib_env_state);
    kdtreebuildtagged(&_xy, &_tags, *n, *nx, *ny, *normtype, &(*kdt)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryknn(const char **errormsg, ae_int_t* result, x_kdtree** kdt, x_vector* x, ae_int_t* k, ae_bool* selfmatch)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *result = kdtreequeryknn(&(*kdt)->obj, &_x, *k, *selfmatch, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryrnn(const char **errormsg, ae_int_t* result, x_kdtree** kdt, x_vector* x, double* r, ae_bool* selfmatch)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *result = kdtreequeryrnn(&(*kdt)->obj, &_x, *r, *selfmatch, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryaknn(const char **errormsg, ae_int_t* result, x_kdtree** kdt, x_vector* x, ae_int_t* k, ae_bool* selfmatch, double* eps)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *result = kdtreequeryaknn(&(*kdt)->obj, &_x, *k, *selfmatch, *eps, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsx(const char **errormsg, x_kdtree** kdt, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    kdtreequeryresultsx(&(*kdt)->obj, &_x, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsxy(const char **errormsg, x_kdtree** kdt, x_matrix* xy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    kdtreequeryresultsxy(&(*kdt)->obj, &_xy, &_alglib_env_state);
    ae_x_set_matrix(xy, &_xy, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultstags(const char **errormsg, x_kdtree** kdt, x_vector* tags)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _tags;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_tags, tags, &_alglib_env_state, ae_true);
    kdtreequeryresultstags(&(*kdt)->obj, &_tags, &_alglib_env_state);
    ae_x_set_vector(tags, &_tags, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsdistances(const char **errormsg, x_kdtree** kdt, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_r, r, &_alglib_env_state, ae_true);
    kdtreequeryresultsdistances(&(*kdt)->obj, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsxi(const char **errormsg, x_kdtree** kdt, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    kdtreequeryresultsxi(&(*kdt)->obj, &_x, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsxyi(const char **errormsg, x_kdtree** kdt, x_matrix* xy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_xy, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    kdtreequeryresultsxyi(&(*kdt)->obj, &_xy, &_alglib_env_state);
    ae_x_set_matrix(xy, &_xy, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultstagsi(const char **errormsg, x_kdtree** kdt, x_vector* tags)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _tags;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_tags, 0, DT_INT, &_alglib_env_state, ae_true);
    kdtreequeryresultstagsi(&(*kdt)->obj, &_tags, &_alglib_env_state);
    ae_x_set_vector(tags, &_tags, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kdtreequeryresultsdistancesi(const char **errormsg, x_kdtree** kdt, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_r, 0, DT_REAL, &_alglib_env_state, ae_true);
    kdtreequeryresultsdistancesi(&(*kdt)->obj, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixtranspose(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_matrix* b, ae_int_t* ib, ae_int_t* jb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    cmatrixtranspose(*m, *n, &_a, *ia, *ja, &_b, *ib, *jb, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixtranspose(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_matrix* b, ae_int_t* ib, ae_int_t* jb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    rmatrixtranspose(*m, *n, &_a, *ia, *ja, &_b, *ib, *jb, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixcopy(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_matrix* b, ae_int_t* ib, ae_int_t* jb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    cmatrixcopy(*m, *n, &_a, *ia, *ja, &_b, *ib, *jb, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixcopy(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_matrix* b, ae_int_t* ib, ae_int_t* jb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    rmatrixcopy(*m, *n, &_a, *ia, *ja, &_b, *ib, *jb, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrank1(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_vector* u, ae_int_t* iu, x_vector* v, ae_int_t* iv)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _u;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_u, u, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_v, v, &_alglib_env_state, ae_true);
    cmatrixrank1(*m, *n, &_a, *ia, *ja, &_u, *iu, &_v, *iv, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(u, &_u, &_alglib_env_state);
    ae_x_set_vector(v, &_v, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrank1(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, x_vector* u, ae_int_t* iu, x_vector* v, ae_int_t* iv)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _u;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_u, u, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_v, v, &_alglib_env_state, ae_true);
    rmatrixrank1(*m, *n, &_a, *ia, *ja, &_u, *iu, &_v, *iv, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(u, &_u, &_alglib_env_state);
    ae_x_set_vector(v, &_v, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixmv(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* opa, x_vector* x, ae_int_t* ix, x_vector* y, ae_int_t* iy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    cmatrixmv(*m, *n, &_a, *ia, *ja, *opa, &_x, *ix, &_y, *iy, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixmv(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* opa, x_vector* x, ae_int_t* ix, x_vector* y, ae_int_t* iy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    rmatrixmv(*m, *n, &_a, *ia, *ja, *opa, &_x, *ix, &_y, *iy, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrighttrsm(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* i1, ae_int_t* j1, ae_bool* isupper, ae_bool* isunit, ae_int_t* optype, x_matrix* x, ae_int_t* i2, ae_int_t* j2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    cmatrixrighttrsm(*m, *n, &_a, *i1, *j1, *isupper, *isunit, *optype, &_x, *i2, *j2, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlefttrsm(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* i1, ae_int_t* j1, ae_bool* isupper, ae_bool* isunit, ae_int_t* optype, x_matrix* x, ae_int_t* i2, ae_int_t* j2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    cmatrixlefttrsm(*m, *n, &_a, *i1, *j1, *isupper, *isunit, *optype, &_x, *i2, *j2, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrighttrsm(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* i1, ae_int_t* j1, ae_bool* isupper, ae_bool* isunit, ae_int_t* optype, x_matrix* x, ae_int_t* i2, ae_int_t* j2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    rmatrixrighttrsm(*m, *n, &_a, *i1, *j1, *isupper, *isunit, *optype, &_x, *i2, *j2, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlefttrsm(const char **errormsg, ae_int_t* m, ae_int_t* n, x_matrix* a, ae_int_t* i1, ae_int_t* j1, ae_bool* isupper, ae_bool* isunit, ae_int_t* optype, x_matrix* x, ae_int_t* i2, ae_int_t* j2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    rmatrixlefttrsm(*m, *n, &_a, *i1, *j1, *isupper, *isunit, *optype, &_x, *i2, *j2, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixsyrk(const char **errormsg, ae_int_t* n, ae_int_t* k, double* alpha, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* optypea, double* beta, x_matrix* c, ae_int_t* ic, ae_int_t* jc, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    cmatrixsyrk(*n, *k, *alpha, &_a, *ia, *ja, *optypea, *beta, &_c, *ic, *jc, *isupper, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixsyrk(const char **errormsg, ae_int_t* n, ae_int_t* k, double* alpha, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* optypea, double* beta, x_matrix* c, ae_int_t* ic, ae_int_t* jc, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    rmatrixsyrk(*n, *k, *alpha, &_a, *ia, *ja, *optypea, *beta, &_c, *ic, *jc, *isupper, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixgemm(const char **errormsg, ae_int_t* m, ae_int_t* n, ae_int_t* k, ae_complex* alpha, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* optypea, x_matrix* b, ae_int_t* ib, ae_int_t* jb, ae_int_t* optypeb, ae_complex* beta, x_matrix* c, ae_int_t* ic, ae_int_t* jc)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    cmatrixgemm(*m, *n, *k, *alpha, &_a, *ia, *ja, *optypea, &_b, *ib, *jb, *optypeb, *beta, &_c, *ic, *jc, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixgemm(const char **errormsg, ae_int_t* m, ae_int_t* n, ae_int_t* k, double* alpha, x_matrix* a, ae_int_t* ia, ae_int_t* ja, ae_int_t* optypea, x_matrix* b, ae_int_t* ib, ae_int_t* jb, ae_int_t* optypeb, double* beta, x_matrix* c, ae_int_t* ic, ae_int_t* jc)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    rmatrixgemm(*m, *n, *k, *alpha, &_a, *ia, *ja, *optypea, &_b, *ib, *jb, *optypeb, *beta, &_c, *ic, *jc, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_samplemoments(const char **errormsg, x_vector* x, ae_int_t* n, double* mean, double* variance, double* skewness, double* kurtosis)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    samplemoments(&_x, *n, mean, variance, skewness, kurtosis, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_sampleadev(const char **errormsg, x_vector* x, ae_int_t* n, double* adev)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    sampleadev(&_x, *n, adev, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_samplemedian(const char **errormsg, x_vector* x, ae_int_t* n, double* median)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    samplemedian(&_x, *n, median, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_samplepercentile(const char **errormsg, x_vector* x, ae_int_t* n, double* p, double* v)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    samplepercentile(&_x, *n, *p, v, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cov2(const char **errormsg, double* result, x_vector* x, x_vector* y, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *result = cov2(&_x, &_y, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pearsoncorr2(const char **errormsg, double* result, x_vector* x, x_vector* y, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *result = pearsoncorr2(&_x, &_y, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spearmancorr2(const char **errormsg, double* result, x_vector* x, x_vector* y, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *result = spearmancorr2(&_x, &_y, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_covm(const char **errormsg, x_matrix* x, ae_int_t* n, ae_int_t* m, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    covm(&_x, *n, *m, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pearsoncorrm(const char **errormsg, x_matrix* x, ae_int_t* n, ae_int_t* m, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    pearsoncorrm(&_x, *n, *m, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spearmancorrm(const char **errormsg, x_matrix* x, ae_int_t* n, ae_int_t* m, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spearmancorrm(&_x, *n, *m, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_covm2(const char **errormsg, x_matrix* x, x_matrix* y, ae_int_t* n, ae_int_t* m1, ae_int_t* m2, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _y;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    covm2(&_x, &_y, *n, *m1, *m2, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pearsoncorrm2(const char **errormsg, x_matrix* x, x_matrix* y, ae_int_t* n, ae_int_t* m1, ae_int_t* m2, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _y;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    pearsoncorrm2(&_x, &_y, *n, *m1, *m2, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spearmancorrm2(const char **errormsg, x_matrix* x, x_matrix* y, ae_int_t* n, ae_int_t* m1, ae_int_t* m2, x_matrix* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_matrix _y;
    ae_matrix _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spearmancorrm2(&_x, &_y, *n, *m1, *m2, &_c, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pearsoncorrelation(const char **errormsg, double* result, x_vector* x, x_vector* y, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *result = pearsoncorrelation(&_x, &_y, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spearmanrankcorrelation(const char **errormsg, double* result, x_vector* x, x_vector* y, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *result = spearmanrankcorrelation(&_x, &_y, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dsoptimalsplit2(const char **errormsg, x_vector* a, x_vector* c, ae_int_t* n, ae_int_t* info, double* threshold, double* pal, double* pbl, double* par, double* pbr, double* cve)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    dsoptimalsplit2(&_a, &_c, *n, info, threshold, pal, pbl, par, pbr, cve, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dsoptimalsplit2fast(const char **errormsg, x_vector* a, x_vector* c, x_vector* tiesbuf, x_vector* cntbuf, x_vector* bufr, x_vector* bufi, ae_int_t* n, ae_int_t* nc, double* alpha, ae_int_t* info, double* threshold, double* rms, double* cvrms)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _c;
    ae_vector _tiesbuf;
    ae_vector _cntbuf;
    ae_vector _bufr;
    ae_vector _bufi;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tiesbuf, tiesbuf, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_cntbuf, cntbuf, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bufr, bufr, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bufi, bufi, &_alglib_env_state, ae_true);
    dsoptimalsplit2fast(&_a, &_c, &_tiesbuf, &_cntbuf, &_bufr, &_bufi, *n, *nc, *alpha, info, threshold, rms, cvrms, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_x_set_vector(tiesbuf, &_tiesbuf, &_alglib_env_state);
    ae_x_set_vector(cntbuf, &_cntbuf, &_alglib_env_state);
    ae_x_set_vector(bufr, &_bufr, &_alglib_env_state);
    ae_x_set_vector(bufi, &_bufi, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    decisionforest obj;
} x_decisionforest;
x_decisionforest* x_obj_alloc_decisionforest(ae_state *_state)
{
    x_decisionforest *result;
    result = ae_malloc(sizeof(x_decisionforest), _state);
    _decisionforest_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_decisionforest(x_decisionforest *obj)
{
    if( obj==NULL )
        return;
    _decisionforest_clear(&obj->obj);
    ae_free(obj);
    return;
}
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double oobrelclserror;
    double oobavgce;
    double oobrmserror;
    double oobavgerror;
    double oobavgrelerror;
} x_dfreport;
void x_set_dfreport(x_dfreport *dst, dfreport *src, ae_state *_state)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->oobrelclserror = src->oobrelclserror;
    dst->oobavgce = src->oobavgce;
    dst->oobrmserror = src->oobrmserror;
    dst->oobavgerror = src->oobavgerror;
    dst->oobavgrelerror = src->oobavgrelerror;
}
void dfreport_init_from_x(dfreport *dst, x_dfreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->oobrelclserror = src->oobrelclserror;
    dst->oobavgce = src->oobavgce;
    dst->oobrmserror = src->oobrmserror;
    dst->oobavgerror = src->oobavgerror;
    dst->oobavgrelerror = src->oobavgrelerror;
}
DLLEXPORT int alglib_dfbuildrandomdecisionforest(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* nclasses, ae_int_t* ntrees, double* r, ae_int_t* info, x_decisionforest** df, x_dfreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    dfreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *df = x_obj_alloc_decisionforest(&_alglib_env_state);
    _dfreport_init(&_rep, &_alglib_env_state, ae_true);
    dfbuildrandomdecisionforest(&_xy, *npoints, *nvars, *nclasses, *ntrees, *r, info, &(*df)->obj, &_rep, &_alglib_env_state);
    x_set_dfreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfprocess(const char **errormsg, x_decisionforest** df, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    dfprocess(&(*df)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfprocessi(const char **errormsg, x_decisionforest** df, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init(&_y, 0, DT_REAL, &_alglib_env_state, ae_true);
    dfprocessi(&(*df)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfrelclserror(const char **errormsg, double* result, x_decisionforest** df, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = dfrelclserror(&(*df)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfavgce(const char **errormsg, double* result, x_decisionforest** df, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = dfavgce(&(*df)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfrmserror(const char **errormsg, double* result, x_decisionforest** df, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = dfrmserror(&(*df)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfavgerror(const char **errormsg, double* result, x_decisionforest** df, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = dfavgerror(&(*df)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dfavgrelerror(const char **errormsg, double* result, x_decisionforest** df, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = dfavgrelerror(&(*df)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_kmeansgenerate(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* k, ae_int_t* restarts, ae_int_t* info, x_matrix* c, x_vector* xyc)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_matrix _c;
    ae_vector _xyc;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_matrix_init(&_c, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_xyc, 0, DT_INT, &_alglib_env_state, ae_true);
    kmeansgenerate(&_xy, *npoints, *nvars, *k, *restarts, info, &_c, &_xyc, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_x_set_vector(xyc, &_xyc, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixqr(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixqr(&_a, *m, *n, &_tau, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixlq(&_a, *m, *n, &_tau, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixqr(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixqr(&_a, *m, *n, &_tau, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixlq(&_a, *m, *n, &_tau, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixqrunpackq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau, ae_int_t* qcolumns, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixqrunpackq(&_a, *m, *n, &_tau, *qcolumns, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixqrunpackr(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_matrix* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_r, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixqrunpackr(&_a, *m, *n, &_r, &_alglib_env_state);
    ae_x_set_matrix(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlqunpackq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau, ae_int_t* qrows, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixlqunpackq(&_a, *m, *n, &_tau, *qrows, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlqunpackl(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_matrix* l)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _l;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_l, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixlqunpackl(&_a, *m, *n, &_l, &_alglib_env_state);
    ae_x_set_matrix(l, &_l, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixqrunpackq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau, ae_int_t* qcolumns, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixqrunpackq(&_a, *m, *n, &_tau, *qcolumns, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixqrunpackr(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_matrix* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_r, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixqrunpackr(&_a, *m, *n, &_r, &_alglib_env_state);
    ae_x_set_matrix(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlqunpackq(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tau, ae_int_t* qrows, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixlqunpackq(&_a, *m, *n, &_tau, *qrows, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlqunpackl(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_matrix* l)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _l;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_l, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixlqunpackl(&_a, *m, *n, &_l, &_alglib_env_state);
    ae_x_set_matrix(l, &_l, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbd(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* tauq, x_vector* taup)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tauq;
    ae_vector _taup;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tauq, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_taup, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixbd(&_a, *m, *n, &_tauq, &_taup, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tauq, &_tauq, &_alglib_env_state);
    ae_x_set_vector(taup, &_taup, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdunpackq(const char **errormsg, x_matrix* qp, ae_int_t* m, ae_int_t* n, x_vector* tauq, ae_int_t* qcolumns, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _qp;
    ae_vector _tauq;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_qp, qp, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tauq, tauq, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixbdunpackq(&_qp, *m, *n, &_tauq, *qcolumns, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdmultiplybyq(const char **errormsg, x_matrix* qp, ae_int_t* m, ae_int_t* n, x_vector* tauq, x_matrix* z, ae_int_t* zrows, ae_int_t* zcolumns, ae_bool* fromtheright, ae_bool* dotranspose)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _qp;
    ae_vector _tauq;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_qp, qp, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tauq, tauq, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_z, z, &_alglib_env_state, ae_true);
    rmatrixbdmultiplybyq(&_qp, *m, *n, &_tauq, &_z, *zrows, *zcolumns, *fromtheright, *dotranspose, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdunpackpt(const char **errormsg, x_matrix* qp, ae_int_t* m, ae_int_t* n, x_vector* taup, ae_int_t* ptrows, x_matrix* pt)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _qp;
    ae_vector _taup;
    ae_matrix _pt;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_qp, qp, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_taup, taup, &_alglib_env_state, ae_true);
    ae_matrix_init(&_pt, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixbdunpackpt(&_qp, *m, *n, &_taup, *ptrows, &_pt, &_alglib_env_state);
    ae_x_set_matrix(pt, &_pt, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdmultiplybyp(const char **errormsg, x_matrix* qp, ae_int_t* m, ae_int_t* n, x_vector* taup, x_matrix* z, ae_int_t* zrows, ae_int_t* zcolumns, ae_bool* fromtheright, ae_bool* dotranspose)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _qp;
    ae_vector _taup;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_qp, qp, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_taup, taup, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_z, z, &_alglib_env_state, ae_true);
    rmatrixbdmultiplybyp(&_qp, *m, *n, &_taup, &_z, *zrows, *zcolumns, *fromtheright, *dotranspose, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdunpackdiagonals(const char **errormsg, x_matrix* b, ae_int_t* m, ae_int_t* n, ae_bool* isupper, x_vector* d, x_vector* e)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _b;
    ae_vector _d;
    ae_vector _e;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_e, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixbdunpackdiagonals(&_b, *m, *n, isupper, &_d, &_e, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_vector(e, &_e, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixhessenberg(const char **errormsg, x_matrix* a, ae_int_t* n, x_vector* tau)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixhessenberg(&_a, *n, &_tau, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixhessenbergunpackq(const char **errormsg, x_matrix* a, ae_int_t* n, x_vector* tau, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixhessenbergunpackq(&_a, *n, &_tau, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixhessenbergunpackh(const char **errormsg, x_matrix* a, ae_int_t* n, x_matrix* h)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _h;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_h, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixhessenbergunpackh(&_a, *n, &_h, &_alglib_env_state);
    ae_x_set_matrix(h, &_h, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixtd(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* tau, x_vector* d, x_vector* e)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_vector _d;
    ae_vector _e;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_e, 0, DT_REAL, &_alglib_env_state, ae_true);
    smatrixtd(&_a, *n, *isupper, &_tau, &_d, &_e, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_vector(e, &_e, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixtdunpackq(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* tau, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    smatrixtdunpackq(&_a, *n, *isupper, &_tau, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixtd(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* tau, x_vector* d, x_vector* e)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_vector _d;
    ae_vector _e;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_tau, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_e, 0, DT_REAL, &_alglib_env_state, ae_true);
    hmatrixtd(&_a, *n, *isupper, &_tau, &_d, &_e, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(tau, &_tau, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_vector(e, &_e, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixtdunpackq(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* tau, x_matrix* q)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _tau;
    ae_matrix _q;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_tau, tau, &_alglib_env_state, ae_true);
    ae_matrix_init(&_q, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hmatrixtdunpackq(&_a, *n, *isupper, &_tau, &_q, &_alglib_env_state);
    ae_x_set_matrix(q, &_q, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixevd(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, x_vector* d, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _d;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = smatrixevd(&_a, *n, *zneeded, *isupper, &_d, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixevdr(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, double* b1, double* b2, ae_int_t* m, x_vector* w, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _w;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = smatrixevdr(&_a, *n, *zneeded, *isupper, *b1, *b2, m, &_w, &_z, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixevdi(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, ae_int_t* i1, ae_int_t* i2, x_vector* w, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _w;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = smatrixevdi(&_a, *n, *zneeded, *isupper, *i1, *i2, &_w, &_z, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixevd(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, x_vector* d, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _d;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    *result = hmatrixevd(&_a, *n, *zneeded, *isupper, &_d, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixevdr(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, double* b1, double* b2, ae_int_t* m, x_vector* w, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _w;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    *result = hmatrixevdr(&_a, *n, *zneeded, *isupper, *b1, *b2, m, &_w, &_z, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixevdi(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* zneeded, ae_bool* isupper, ae_int_t* i1, ae_int_t* i2, x_vector* w, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _w;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    *result = hmatrixevdi(&_a, *n, *zneeded, *isupper, *i1, *i2, &_w, &_z, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixtdevd(const char **errormsg, ae_bool* result, x_vector* d, x_vector* e, ae_int_t* n, ae_int_t* zneeded, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _d;
    ae_vector _e;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_d, d, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_e, e, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_z, z, &_alglib_env_state, ae_true);
    *result = smatrixtdevd(&_d, &_e, *n, *zneeded, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixtdevdr(const char **errormsg, ae_bool* result, x_vector* d, x_vector* e, ae_int_t* n, ae_int_t* zneeded, double* a, double* b, ae_int_t* m, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _d;
    ae_vector _e;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_d, d, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_e, e, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_z, z, &_alglib_env_state, ae_true);
    *result = smatrixtdevdr(&_d, &_e, *n, *zneeded, *a, *b, m, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixtdevdi(const char **errormsg, ae_bool* result, x_vector* d, x_vector* e, ae_int_t* n, ae_int_t* zneeded, ae_int_t* i1, ae_int_t* i2, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _d;
    ae_vector _e;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_d, d, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_e, e, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_z, z, &_alglib_env_state, ae_true);
    *result = smatrixtdevdi(&_d, &_e, *n, *zneeded, *i1, *i2, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixevd(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_int_t* vneeded, x_vector* wr, x_vector* wi, x_matrix* vl, x_matrix* vr)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _wr;
    ae_vector _wi;
    ae_matrix _vl;
    ae_matrix _vr;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_wr, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wi, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_vl, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_vr, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = rmatrixevd(&_a, *n, *vneeded, &_wr, &_wi, &_vl, &_vr, &_alglib_env_state);
    ae_x_set_vector(wr, &_wr, &_alglib_env_state);
    ae_x_set_vector(wi, &_wi, &_alglib_env_state);
    ae_x_set_matrix(vl, &_vl, &_alglib_env_state);
    ae_x_set_matrix(vr, &_vr, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrndorthogonal(const char **errormsg, ae_int_t* n, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixrndorthogonal(*n, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrndorthogonal(const char **errormsg, ae_int_t* n, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixrndorthogonal(*n, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    smatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spdmatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hmatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixrndcond(const char **errormsg, ae_int_t* n, double* c, x_matrix* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hpdmatrixrndcond(*n, *c, &_a, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrndorthogonalfromtheright(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    rmatrixrndorthogonalfromtheright(&_a, *m, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrndorthogonalfromtheleft(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    rmatrixrndorthogonalfromtheleft(&_a, *m, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrndorthogonalfromtheright(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    cmatrixrndorthogonalfromtheright(&_a, *m, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrndorthogonalfromtheleft(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    cmatrixrndorthogonalfromtheleft(&_a, *m, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixrndmultiply(const char **errormsg, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    smatrixrndmultiply(&_a, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hmatrixrndmultiply(const char **errormsg, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    hmatrixrndmultiply(&_a, *n, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlu(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* pivots)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_pivots, 0, DT_INT, &_alglib_env_state, ae_true);
    rmatrixlu(&_a, *m, *n, &_pivots, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(pivots, &_pivots, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlu(const char **errormsg, x_matrix* a, ae_int_t* m, ae_int_t* n, x_vector* pivots)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_pivots, 0, DT_INT, &_alglib_env_state, ae_true);
    cmatrixlu(&_a, *m, *n, &_pivots, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_vector(pivots, &_pivots, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixcholesky(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = hpdmatrixcholesky(&_a, *n, *isupper, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholesky(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = spdmatrixcholesky(&_a, *n, *isupper, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrcond1(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = rmatrixrcond1(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixrcondinf(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = rmatrixrcondinf(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixrcond(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = spdmatrixrcond(&_a, *n, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixtrrcond1(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = rmatrixtrrcond1(&_a, *n, *isupper, *isunit, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixtrrcondinf(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = rmatrixtrrcondinf(&_a, *n, *isupper, *isunit, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixrcond(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = hpdmatrixrcond(&_a, *n, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrcond1(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = cmatrixrcond1(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixrcondinf(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = cmatrixrcondinf(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlurcond1(const char **errormsg, double* result, x_matrix* lua, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    *result = rmatrixlurcond1(&_lua, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlurcondinf(const char **errormsg, double* result, x_matrix* lua, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    *result = rmatrixlurcondinf(&_lua, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholeskyrcond(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = spdmatrixcholeskyrcond(&_a, *n, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixcholeskyrcond(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = hpdmatrixcholeskyrcond(&_a, *n, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlurcond1(const char **errormsg, double* result, x_matrix* lua, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    *result = cmatrixlurcond1(&_lua, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlurcondinf(const char **errormsg, double* result, x_matrix* lua, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    *result = cmatrixlurcondinf(&_lua, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixtrrcond1(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = cmatrixtrrcond1(&_a, *n, *isupper, *isunit, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixtrrcondinf(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = cmatrixtrrcondinf(&_a, *n, *isupper, *isunit, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    double r1;
    double rinf;
} x_matinvreport;
void x_set_matinvreport(x_matinvreport *dst, matinvreport *src, ae_state *_state)
{
    dst->r1 = src->r1;
    dst->rinf = src->rinf;
}
void matinvreport_init_from_x(matinvreport *dst, x_matinvreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->r1 = src->r1;
    dst->rinf = src->rinf;
}
DLLEXPORT int alglib_rmatrixluinverse(const char **errormsg, x_matrix* a, x_vector* pivots, ae_int_t* n, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pivots, pivots, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    rmatrixluinverse(&_a, &_pivots, *n, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    rmatrixinverse(&_a, *n, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixluinverse(const char **errormsg, x_matrix* a, x_vector* pivots, ae_int_t* n, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pivots, pivots, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    cmatrixluinverse(&_a, &_pivots, *n, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    cmatrixinverse(&_a, *n, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholeskyinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    spdmatrixcholeskyinverse(&_a, *n, *isupper, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    spdmatrixinverse(&_a, *n, *isupper, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixcholeskyinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    hpdmatrixcholeskyinverse(&_a, *n, *isupper, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    hpdmatrixinverse(&_a, *n, *isupper, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixtrinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    rmatrixtrinverse(&_a, *n, *isupper, *isunit, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixtrinverse(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, ae_bool* isunit, ae_int_t* info, x_matinvreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    matinvreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    _matinvreport_init(&_rep, &_alglib_env_state, ae_true);
    cmatrixtrinverse(&_a, *n, *isupper, *isunit, info, &_rep, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    x_set_matinvreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fisherlda(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* nclasses, ae_int_t* info, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    fisherlda(&_xy, *npoints, *nvars, *nclasses, info, &_w, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fisherldan(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* nclasses, ae_int_t* info, x_matrix* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_matrix _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_matrix_init(&_w, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    fisherldan(&_xy, *npoints, *nvars, *nclasses, info, &_w, &_alglib_env_state);
    ae_x_set_matrix(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gammafunction(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = gammafunction(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lngamma(const char **errormsg, double* result, double* x, double* sgngam)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = lngamma(*x, sgngam, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_errorfunction(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = errorfunction(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_errorfunctionc(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = errorfunctionc(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_normaldistribution(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = normaldistribution(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_inverf(const char **errormsg, double* result, double* e)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = inverf(*e, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invnormaldistribution(const char **errormsg, double* result, double* y0)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invnormaldistribution(*y0, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_incompletegamma(const char **errormsg, double* result, double* a, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = incompletegamma(*a, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_incompletegammac(const char **errormsg, double* result, double* a, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = incompletegammac(*a, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invincompletegammac(const char **errormsg, double* result, double* a, double* y0)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invincompletegammac(*a, *y0, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixbdsvd(const char **errormsg, ae_bool* result, x_vector* d, x_vector* e, ae_int_t* n, ae_bool* isupper, ae_bool* isfractionalaccuracyrequired, x_matrix* u, ae_int_t* nru, x_matrix* c, ae_int_t* ncc, x_matrix* vt, ae_int_t* ncvt)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _d;
    ae_vector _e;
    ae_matrix _u;
    ae_matrix _c;
    ae_matrix _vt;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_d, d, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_e, e, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_u, u, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_vt, vt, &_alglib_env_state, ae_true);
    *result = rmatrixbdsvd(&_d, &_e, *n, *isupper, *isfractionalaccuracyrequired, &_u, *nru, &_c, *ncc, &_vt, *ncvt, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(u, &_u, &_alglib_env_state);
    ae_x_set_matrix(c, &_c, &_alglib_env_state);
    ae_x_set_matrix(vt, &_vt, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixsvd(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* m, ae_int_t* n, ae_int_t* uneeded, ae_int_t* vtneeded, ae_int_t* additionalmemory, x_vector* w, x_matrix* u, x_matrix* vt)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _w;
    ae_matrix _u;
    ae_matrix _vt;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_u, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_vt, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = rmatrixsvd(&_a, *m, *n, *uneeded, *vtneeded, *additionalmemory, &_w, &_u, &_vt, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_x_set_matrix(u, &_u, &_alglib_env_state);
    ae_x_set_matrix(vt, &_vt, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    linearmodel obj;
} x_linearmodel;
x_linearmodel* x_obj_alloc_linearmodel(ae_state *_state)
{
    x_linearmodel *result;
    result = ae_malloc(sizeof(x_linearmodel), _state);
    _linearmodel_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_linearmodel(x_linearmodel *obj)
{
    if( obj==NULL )
        return;
    _linearmodel_clear(&obj->obj);
    ae_free(obj);
    return;
}
typedef struct
{
    x_matrix c;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    ae_int_t ncvdefects;
    x_vector cvdefects;
} x_lrreport;
void x_set_lrreport(x_lrreport *dst, lrreport *src, ae_state *_state)
{
    ae_x_set_matrix(&dst->c, &src->c, _state);
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->cvrmserror = src->cvrmserror;
    dst->cvavgerror = src->cvavgerror;
    dst->cvavgrelerror = src->cvavgrelerror;
    dst->ncvdefects = src->ncvdefects;
    ae_x_set_vector(&dst->cvdefects, &src->cvdefects, _state);
}
void lrreport_init_from_x(lrreport *dst, x_lrreport *src, ae_state *_state, ae_bool make_automatic)
{
    ae_matrix_init_from_x(&dst->c, &src->c, _state, make_automatic);
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->cvrmserror = src->cvrmserror;
    dst->cvavgerror = src->cvavgerror;
    dst->cvavgrelerror = src->cvavgrelerror;
    dst->ncvdefects = src->ncvdefects;
    ae_vector_init_from_x(&dst->cvdefects, &src->cvdefects, _state, make_automatic);
}
DLLEXPORT int alglib_lrbuild(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* info, x_linearmodel** lm, x_lrreport* ar)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    lrreport _ar;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_linearmodel(&_alglib_env_state);
    _lrreport_init(&_ar, &_alglib_env_state, ae_true);
    lrbuild(&_xy, *npoints, *nvars, info, &(*lm)->obj, &_ar, &_alglib_env_state);
    x_set_lrreport(ar, &_ar, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrbuilds(const char **errormsg, x_matrix* xy, x_vector* s, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* info, x_linearmodel** lm, x_lrreport* ar)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _s;
    lrreport _ar;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_s, s, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_linearmodel(&_alglib_env_state);
    _lrreport_init(&_ar, &_alglib_env_state, ae_true);
    lrbuilds(&_xy, &_s, *npoints, *nvars, info, &(*lm)->obj, &_ar, &_alglib_env_state);
    x_set_lrreport(ar, &_ar, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrbuildzs(const char **errormsg, x_matrix* xy, x_vector* s, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* info, x_linearmodel** lm, x_lrreport* ar)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _s;
    lrreport _ar;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_s, s, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_linearmodel(&_alglib_env_state);
    _lrreport_init(&_ar, &_alglib_env_state, ae_true);
    lrbuildzs(&_xy, &_s, *npoints, *nvars, info, &(*lm)->obj, &_ar, &_alglib_env_state);
    x_set_lrreport(ar, &_ar, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrbuildz(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* info, x_linearmodel** lm, x_lrreport* ar)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    lrreport _ar;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_linearmodel(&_alglib_env_state);
    _lrreport_init(&_ar, &_alglib_env_state, ae_true);
    lrbuildz(&_xy, *npoints, *nvars, info, &(*lm)->obj, &_ar, &_alglib_env_state);
    x_set_lrreport(ar, &_ar, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrunpack(const char **errormsg, x_linearmodel** lm, x_vector* v, ae_int_t* nvars)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_v, 0, DT_REAL, &_alglib_env_state, ae_true);
    lrunpack(&(*lm)->obj, &_v, nvars, &_alglib_env_state);
    ae_x_set_vector(v, &_v, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrpack(const char **errormsg, x_vector* v, ae_int_t* nvars, x_linearmodel** lm)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_v, v, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_linearmodel(&_alglib_env_state);
    lrpack(&_v, *nvars, &(*lm)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrprocess(const char **errormsg, double* result, x_linearmodel** lm, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *result = lrprocess(&(*lm)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lrrmserror(const char **errormsg, double* result, x_linearmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = lrrmserror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lravgerror(const char **errormsg, double* result, x_linearmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = lravgerror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lravgrelerror(const char **errormsg, double* result, x_linearmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = lravgrelerror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    multilayerperceptron obj;
} x_multilayerperceptron;
x_multilayerperceptron* x_obj_alloc_multilayerperceptron(ae_state *_state)
{
    x_multilayerperceptron *result;
    result = ae_malloc(sizeof(x_multilayerperceptron), _state);
    _multilayerperceptron_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_multilayerperceptron(x_multilayerperceptron *obj)
{
    if( obj==NULL )
        return;
    _multilayerperceptron_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_mlpcreate0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreate0(*nin, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreate1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreate1(*nin, *nhid, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreate2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreate2(*nin, *nhid1, *nhid2, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreateb0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, double* b, double* d, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreateb0(*nin, *nout, *b, *d, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreateb1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, double* b, double* d, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreateb1(*nin, *nhid, *nout, *b, *d, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreateb2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, double* b, double* d, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreateb2(*nin, *nhid1, *nhid2, *nout, *b, *d, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreater0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, double* a, double* b, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreater0(*nin, *nout, *a, *b, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreater1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, double* a, double* b, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreater1(*nin, *nhid, *nout, *a, *b, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreater2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, double* a, double* b, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreater2(*nin, *nhid1, *nhid2, *nout, *a, *b, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreatec0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreatec0(*nin, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreatec1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreatec1(*nin, *nhid, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpcreatec2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *network = x_obj_alloc_multilayerperceptron(&_alglib_env_state);
    mlpcreatec2(*nin, *nhid1, *nhid2, *nout, &(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlprandomize(const char **errormsg, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mlprandomize(&(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlprandomizefull(const char **errormsg, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mlprandomizefull(&(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpproperties(const char **errormsg, x_multilayerperceptron** network, ae_int_t* nin, ae_int_t* nout, ae_int_t* wcount)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mlpproperties(&(*network)->obj, nin, nout, wcount, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpissoftmax(const char **errormsg, ae_bool* result, x_multilayerperceptron** network)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = mlpissoftmax(&(*network)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpprocess(const char **errormsg, x_multilayerperceptron** network, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    mlpprocess(&(*network)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpprocessi(const char **errormsg, x_multilayerperceptron** network, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init(&_y, 0, DT_REAL, &_alglib_env_state, ae_true);
    mlpprocessi(&(*network)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlperror(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlperror(&(*network)->obj, &_xy, *ssize, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlperrorn(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlperrorn(&(*network)->obj, &_xy, *ssize, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpclserror(const char **errormsg, ae_int_t* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpclserror(&(*network)->obj, &_xy, *ssize, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlprelclserror(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlprelclserror(&(*network)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpavgce(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpavgce(&(*network)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlprmserror(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlprmserror(&(*network)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpavgerror(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpavgerror(&(*network)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpavgrelerror(const char **errormsg, double* result, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpavgrelerror(&(*network)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpgrad(const char **errormsg, x_multilayerperceptron** network, x_vector* x, x_vector* desiredy, double* e, x_vector* grad)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _desiredy;
    ae_vector _grad;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_desiredy, desiredy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    mlpgrad(&(*network)->obj, &_x, &_desiredy, e, &_grad, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpgradn(const char **errormsg, x_multilayerperceptron** network, x_vector* x, x_vector* desiredy, double* e, x_vector* grad)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _desiredy;
    ae_vector _grad;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_desiredy, desiredy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    mlpgradn(&(*network)->obj, &_x, &_desiredy, e, &_grad, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpgradbatch(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize, double* e, x_vector* grad)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _grad;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    mlpgradbatch(&(*network)->obj, &_xy, *ssize, e, &_grad, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpgradnbatch(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize, double* e, x_vector* grad)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _grad;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    mlpgradnbatch(&(*network)->obj, &_xy, *ssize, e, &_grad, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlphessiannbatch(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize, double* e, x_vector* grad, x_matrix* h)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _grad;
    ae_matrix _h;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_h, h, &_alglib_env_state, ae_true);
    mlphessiannbatch(&(*network)->obj, &_xy, *ssize, e, &_grad, &_h, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_x_set_matrix(h, &_h, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlphessianbatch(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* ssize, double* e, x_vector* grad, x_matrix* h)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_vector _grad;
    ae_matrix _h;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_grad, grad, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_h, h, &_alglib_env_state, ae_true);
    mlphessianbatch(&(*network)->obj, &_xy, *ssize, e, &_grad, &_h, &_alglib_env_state);
    ae_x_set_vector(grad, &_grad, &_alglib_env_state);
    ae_x_set_matrix(h, &_h, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    double r1;
    double rinf;
} x_densesolverreport;
void x_set_densesolverreport(x_densesolverreport *dst, densesolverreport *src, ae_state *_state)
{
    dst->r1 = src->r1;
    dst->rinf = src->rinf;
}
void densesolverreport_init_from_x(densesolverreport *dst, x_densesolverreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->r1 = src->r1;
    dst->rinf = src->rinf;
}
typedef struct
{
    double r2;
    x_matrix cx;
    ae_int_t n;
    ae_int_t k;
} x_densesolverlsreport;
void x_set_densesolverlsreport(x_densesolverlsreport *dst, densesolverlsreport *src, ae_state *_state)
{
    dst->r2 = src->r2;
    ae_x_set_matrix(&dst->cx, &src->cx, _state);
    dst->n = src->n;
    dst->k = src->k;
}
void densesolverlsreport_init_from_x(densesolverlsreport *dst, x_densesolverlsreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->r2 = src->r2;
    ae_matrix_init_from_x(&dst->cx, &src->cx, _state, make_automatic);
    dst->n = src->n;
    dst->k = src->k;
}
DLLEXPORT int alglib_rmatrixsolve(const char **errormsg, x_matrix* a, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixsolve(&_a, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixsolvem(const char **errormsg, x_matrix* a, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_bool* rfs, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixsolvem(&_a, *n, &_b, *m, *rfs, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlusolve(const char **errormsg, x_matrix* lua, x_vector* p, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_vector _p;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixlusolve(&_lua, &_p, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixlusolvem(const char **errormsg, x_matrix* lua, x_vector* p, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_vector _p;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixlusolvem(&_lua, &_p, *n, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixmixedsolve(const char **errormsg, x_matrix* a, x_matrix* lua, x_vector* p, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _lua;
    ae_vector _p;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixmixedsolve(&_a, &_lua, &_p, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixmixedsolvem(const char **errormsg, x_matrix* a, x_matrix* lua, x_vector* p, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _lua;
    ae_vector _p;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixmixedsolvem(&_a, &_lua, &_p, *n, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixsolvem(const char **errormsg, x_matrix* a, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_bool* rfs, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixsolvem(&_a, *n, &_b, *m, *rfs, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixsolve(const char **errormsg, x_matrix* a, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixsolve(&_a, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlusolvem(const char **errormsg, x_matrix* lua, x_vector* p, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_vector _p;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixlusolvem(&_lua, &_p, *n, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixlusolve(const char **errormsg, x_matrix* lua, x_vector* p, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _lua;
    ae_vector _p;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixlusolve(&_lua, &_p, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixmixedsolvem(const char **errormsg, x_matrix* a, x_matrix* lua, x_vector* p, ae_int_t* n, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _lua;
    ae_vector _p;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixmixedsolvem(&_a, &_lua, &_p, *n, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixmixedsolve(const char **errormsg, x_matrix* a, x_matrix* lua, x_vector* p, ae_int_t* n, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _lua;
    ae_vector _p;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_lua, lua, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    cmatrixmixedsolve(&_a, &_lua, &_p, *n, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixsolvem(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spdmatrixsolvem(&_a, *n, *isupper, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixsolve(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    spdmatrixsolve(&_a, *n, *isupper, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholeskysolvem(const char **errormsg, x_matrix* cha, ae_int_t* n, ae_bool* isupper, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _cha;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_cha, cha, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spdmatrixcholeskysolvem(&_cha, *n, *isupper, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholeskysolve(const char **errormsg, x_matrix* cha, ae_int_t* n, ae_bool* isupper, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _cha;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_cha, cha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    spdmatrixcholeskysolve(&_cha, *n, *isupper, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixsolvem(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hpdmatrixsolvem(&_a, *n, *isupper, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixsolve(const char **errormsg, x_matrix* a, ae_int_t* n, ae_bool* isupper, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hpdmatrixsolve(&_a, *n, *isupper, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixcholeskysolvem(const char **errormsg, x_matrix* cha, ae_int_t* n, ae_bool* isupper, x_matrix* b, ae_int_t* m, ae_int_t* info, x_densesolverreport* rep, x_matrix* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _cha;
    ae_matrix _b;
    densesolverreport _rep;
    ae_matrix _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_cha, cha, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_matrix_init(&_x, 0, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hpdmatrixcholeskysolvem(&_cha, *n, *isupper, &_b, *m, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_matrix(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hpdmatrixcholeskysolve(const char **errormsg, x_matrix* cha, ae_int_t* n, ae_bool* isupper, x_vector* b, ae_int_t* info, x_densesolverreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _cha;
    ae_vector _b;
    densesolverreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_cha, cha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    hpdmatrixcholeskysolve(&_cha, *n, *isupper, &_b, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixsolvels(const char **errormsg, x_matrix* a, ae_int_t* nrows, ae_int_t* ncols, x_vector* b, double* threshold, ae_int_t* info, x_densesolverlsreport* rep, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _b;
    densesolverlsreport _rep;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    _densesolverlsreport_init(&_rep, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    rmatrixsolvels(&_a, *nrows, *ncols, &_b, *threshold, info, &_rep, &_x, &_alglib_env_state);
    x_set_densesolverlsreport(rep, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    logitmodel obj;
} x_logitmodel;
x_logitmodel* x_obj_alloc_logitmodel(ae_state *_state)
{
    x_logitmodel *result;
    result = ae_malloc(sizeof(x_logitmodel), _state);
    _logitmodel_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_logitmodel(x_logitmodel *obj)
{
    if( obj==NULL )
        return;
    _logitmodel_clear(&obj->obj);
    ae_free(obj);
    return;
}
typedef struct
{
    ae_int_t ngrad;
    ae_int_t nhess;
} x_mnlreport;
void x_set_mnlreport(x_mnlreport *dst, mnlreport *src, ae_state *_state)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
}
void mnlreport_init_from_x(mnlreport *dst, x_mnlreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
}
DLLEXPORT int alglib_mnltrainh(const char **errormsg, x_matrix* xy, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* nclasses, ae_int_t* info, x_logitmodel** lm, x_mnlreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mnlreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_logitmodel(&_alglib_env_state);
    _mnlreport_init(&_rep, &_alglib_env_state, ae_true);
    mnltrainh(&_xy, *npoints, *nvars, *nclasses, info, &(*lm)->obj, &_rep, &_alglib_env_state);
    x_set_mnlreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlprocess(const char **errormsg, x_logitmodel** lm, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    mnlprocess(&(*lm)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlprocessi(const char **errormsg, x_logitmodel** lm, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init(&_y, 0, DT_REAL, &_alglib_env_state, ae_true);
    mnlprocessi(&(*lm)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlunpack(const char **errormsg, x_logitmodel** lm, x_matrix* a, ae_int_t* nvars, ae_int_t* nclasses)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_a, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    mnlunpack(&(*lm)->obj, &_a, nvars, nclasses, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlpack(const char **errormsg, x_matrix* a, ae_int_t* nvars, ae_int_t* nclasses, x_logitmodel** lm)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *lm = x_obj_alloc_logitmodel(&_alglib_env_state);
    mnlpack(&_a, *nvars, *nclasses, &(*lm)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlavgce(const char **errormsg, double* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlavgce(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlrelclserror(const char **errormsg, double* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlrelclserror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlrmserror(const char **errormsg, double* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlrmserror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlavgerror(const char **errormsg, double* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlavgerror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlavgrelerror(const char **errormsg, double* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* ssize)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlavgrelerror(&(*lm)->obj, &_xy, *ssize, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mnlclserror(const char **errormsg, ae_int_t* result, x_logitmodel** lm, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mnlclserror(&(*lm)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    minlbfgsstate obj;
} x_minlbfgsstate;
x_minlbfgsstate* x_obj_alloc_minlbfgsstate(ae_state *_state)
{
    x_minlbfgsstate *result;
    result = ae_malloc(sizeof(x_minlbfgsstate), _state);
    _minlbfgsstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_minlbfgsstate(x_minlbfgsstate *obj)
{
    if( obj==NULL )
        return;
    _minlbfgsstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_minlbfgsstate_get_needfg(x_minlbfgsstate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_minlbfgsstate_set_needfg(x_minlbfgsstate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_minlbfgsstate_get_xupdated(x_minlbfgsstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_minlbfgsstate_set_xupdated(x_minlbfgsstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_minlbfgsstate_get_f(x_minlbfgsstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_minlbfgsstate_set_f(x_minlbfgsstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_minlbfgsstate_get_g(x_minlbfgsstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_minlbfgsstate_get_x(x_minlbfgsstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
} x_minlbfgsreport;
void x_set_minlbfgsreport(x_minlbfgsreport *dst, minlbfgsreport *src, ae_state *_state)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
void minlbfgsreport_init_from_x(minlbfgsreport *dst, x_minlbfgsreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
DLLEXPORT int alglib_minlbfgscreate(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_minlbfgsstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlbfgsstate(&_alglib_env_state);
    minlbfgscreate(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgssetcond(const char **errormsg, x_minlbfgsstate** state, double* epsg, double* epsf, double* epsx, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlbfgssetcond(&(*state)->obj, *epsg, *epsf, *epsx, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgssetxrep(const char **errormsg, x_minlbfgsstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlbfgssetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgssetstpmax(const char **errormsg, x_minlbfgsstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlbfgssetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgssetdefaultpreconditioner(const char **errormsg, x_minlbfgsstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlbfgssetdefaultpreconditioner(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgssetcholeskypreconditioner(const char **errormsg, x_minlbfgsstate** state, x_matrix* p, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _p;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_p, p, &_alglib_env_state, ae_true);
    minlbfgssetcholeskypreconditioner(&(*state)->obj, &_p, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgsiteration(const char **errormsg, ae_bool* result, x_minlbfgsstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = minlbfgsiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgsresults(const char **errormsg, x_minlbfgsstate** state, x_vector* x, x_minlbfgsreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minlbfgsreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _minlbfgsreport_init(&_rep, &_alglib_env_state, ae_true);
    minlbfgsresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minlbfgsreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgsresultsbuf(const char **errormsg, x_minlbfgsstate** state, x_vector* x, x_minlbfgsreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minlbfgsreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minlbfgsreport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    minlbfgsresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minlbfgsreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlbfgsrestartfrom(const char **errormsg, x_minlbfgsstate** state, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minlbfgsrestartfrom(&(*state)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
} x_mlpreport;
void x_set_mlpreport(x_mlpreport *dst, mlpreport *src, ae_state *_state)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
}
void mlpreport_init_from_x(mlpreport *dst, x_mlpreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
}
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} x_mlpcvreport;
void x_set_mlpcvreport(x_mlpcvreport *dst, mlpcvreport *src, ae_state *_state)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
}
void mlpcvreport_init_from_x(mlpcvreport *dst, x_mlpcvreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
}
DLLEXPORT int alglib_mlptrainlm(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, ae_int_t* info, x_mlpreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    mlptrainlm(&(*network)->obj, &_xy, *npoints, *decay, *restarts, info, &_rep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlptrainlbfgs(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, double* wstep, ae_int_t* maxits, ae_int_t* info, x_mlpreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    mlptrainlbfgs(&(*network)->obj, &_xy, *npoints, *decay, *restarts, *wstep, *maxits, info, &_rep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlptraines(const char **errormsg, x_multilayerperceptron** network, x_matrix* trnxy, ae_int_t* trnsize, x_matrix* valxy, ae_int_t* valsize, double* decay, ae_int_t* restarts, ae_int_t* info, x_mlpreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _trnxy;
    ae_matrix _valxy;
    mlpreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_trnxy, trnxy, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_valxy, valxy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    mlptraines(&(*network)->obj, &_trnxy, *trnsize, &_valxy, *valsize, *decay, *restarts, info, &_rep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpkfoldcvlbfgs(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, double* wstep, ae_int_t* maxits, ae_int_t* foldscount, ae_int_t* info, x_mlpreport* rep, x_mlpcvreport* cvrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    mlpcvreport _cvrep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    _mlpcvreport_init(&_cvrep, &_alglib_env_state, ae_true);
    mlpkfoldcvlbfgs(&(*network)->obj, &_xy, *npoints, *decay, *restarts, *wstep, *maxits, *foldscount, info, &_rep, &_cvrep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    x_set_mlpcvreport(cvrep, &_cvrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpkfoldcvlm(const char **errormsg, x_multilayerperceptron** network, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, ae_int_t* foldscount, ae_int_t* info, x_mlpreport* rep, x_mlpcvreport* cvrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    mlpcvreport _cvrep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    _mlpcvreport_init(&_cvrep, &_alglib_env_state, ae_true);
    mlpkfoldcvlm(&(*network)->obj, &_xy, *npoints, *decay, *restarts, *foldscount, info, &_rep, &_cvrep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    x_set_mlpcvreport(cvrep, &_cvrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    mlpensemble obj;
} x_mlpensemble;
x_mlpensemble* x_obj_alloc_mlpensemble(ae_state *_state)
{
    x_mlpensemble *result;
    result = ae_malloc(sizeof(x_mlpensemble), _state);
    _mlpensemble_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_mlpensemble(x_mlpensemble *obj)
{
    if( obj==NULL )
        return;
    _mlpensemble_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_mlpecreate0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreate0(*nin, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreate1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreate1(*nin, *nhid, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreate2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreate2(*nin, *nhid1, *nhid2, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreateb0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, double* b, double* d, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreateb0(*nin, *nout, *b, *d, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreateb1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, double* b, double* d, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreateb1(*nin, *nhid, *nout, *b, *d, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreateb2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, double* b, double* d, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreateb2(*nin, *nhid1, *nhid2, *nout, *b, *d, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreater0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, double* a, double* b, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreater0(*nin, *nout, *a, *b, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreater1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, double* a, double* b, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreater1(*nin, *nhid, *nout, *a, *b, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreater2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, double* a, double* b, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreater2(*nin, *nhid1, *nhid2, *nout, *a, *b, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreatec0(const char **errormsg, ae_int_t* nin, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreatec0(*nin, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreatec1(const char **errormsg, ae_int_t* nin, ae_int_t* nhid, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreatec1(*nin, *nhid, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreatec2(const char **errormsg, ae_int_t* nin, ae_int_t* nhid1, ae_int_t* nhid2, ae_int_t* nout, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreatec2(*nin, *nhid1, *nhid2, *nout, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpecreatefromnetwork(const char **errormsg, x_multilayerperceptron** network, ae_int_t* ensemblesize, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *ensemble = x_obj_alloc_mlpensemble(&_alglib_env_state);
    mlpecreatefromnetwork(&(*network)->obj, *ensemblesize, &(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlperandomize(const char **errormsg, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mlperandomize(&(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeproperties(const char **errormsg, x_mlpensemble** ensemble, ae_int_t* nin, ae_int_t* nout)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mlpeproperties(&(*ensemble)->obj, nin, nout, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeissoftmax(const char **errormsg, ae_bool* result, x_mlpensemble** ensemble)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = mlpeissoftmax(&(*ensemble)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeprocess(const char **errormsg, x_mlpensemble** ensemble, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    mlpeprocess(&(*ensemble)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeprocessi(const char **errormsg, x_mlpensemble** ensemble, x_vector* x, x_vector* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init(&_y, 0, DT_REAL, &_alglib_env_state, ae_true);
    mlpeprocessi(&(*ensemble)->obj, &_x, &_y, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlperelclserror(const char **errormsg, double* result, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlperelclserror(&(*ensemble)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeavgce(const char **errormsg, double* result, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpeavgce(&(*ensemble)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpermserror(const char **errormsg, double* result, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpermserror(&(*ensemble)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeavgerror(const char **errormsg, double* result, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpeavgerror(&(*ensemble)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpeavgrelerror(const char **errormsg, double* result, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *result = mlpeavgrelerror(&(*ensemble)->obj, &_xy, *npoints, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpebagginglm(const char **errormsg, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, ae_int_t* info, x_mlpreport* rep, x_mlpcvreport* ooberrors)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    mlpcvreport _ooberrors;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    _mlpcvreport_init(&_ooberrors, &_alglib_env_state, ae_true);
    mlpebagginglm(&(*ensemble)->obj, &_xy, *npoints, *decay, *restarts, info, &_rep, &_ooberrors, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    x_set_mlpcvreport(ooberrors, &_ooberrors, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpebagginglbfgs(const char **errormsg, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, double* wstep, ae_int_t* maxits, ae_int_t* info, x_mlpreport* rep, x_mlpcvreport* ooberrors)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    mlpcvreport _ooberrors;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    _mlpcvreport_init(&_ooberrors, &_alglib_env_state, ae_true);
    mlpebagginglbfgs(&(*ensemble)->obj, &_xy, *npoints, *decay, *restarts, *wstep, *maxits, info, &_rep, &_ooberrors, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    x_set_mlpcvreport(ooberrors, &_ooberrors, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mlpetraines(const char **errormsg, x_mlpensemble** ensemble, x_matrix* xy, ae_int_t* npoints, double* decay, ae_int_t* restarts, ae_int_t* info, x_mlpreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    mlpreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    _mlpreport_init(&_rep, &_alglib_env_state, ae_true);
    mlpetraines(&(*ensemble)->obj, &_xy, *npoints, *decay, *restarts, info, &_rep, &_alglib_env_state);
    x_set_mlpreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pcabuildbasis(const char **errormsg, x_matrix* x, ae_int_t* npoints, ae_int_t* nvars, ae_int_t* info, x_vector* s2, x_matrix* v)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _s2;
    ae_matrix _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init(&_s2, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_v, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    pcabuildbasis(&_x, *npoints, *nvars, info, &_s2, &_v, &_alglib_env_state);
    ae_x_set_vector(s2, &_s2, &_alglib_env_state);
    ae_x_set_matrix(v, &_v, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    odesolverstate obj;
} x_odesolverstate;
x_odesolverstate* x_obj_alloc_odesolverstate(ae_state *_state)
{
    x_odesolverstate *result;
    result = ae_malloc(sizeof(x_odesolverstate), _state);
    _odesolverstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_odesolverstate(x_odesolverstate *obj)
{
    if( obj==NULL )
        return;
    _odesolverstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_odesolverstate_get_needdy(x_odesolverstate *obj, ae_bool *result)
{
    *result = obj->obj.needdy;
}
DLLEXPORT void x_odesolverstate_set_needdy(x_odesolverstate *obj, ae_bool *result)
{
    obj->obj.needdy = *result;
}
DLLEXPORT void x_odesolverstate_get_y(x_odesolverstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.y);
}
DLLEXPORT void x_odesolverstate_get_dy(x_odesolverstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.dy);
}
DLLEXPORT void x_odesolverstate_get_x(x_odesolverstate *obj, double *result)
{
    *result = obj->obj.x;
}
DLLEXPORT void x_odesolverstate_set_x(x_odesolverstate *obj, double *result)
{
    obj->obj.x = *result;
}
typedef struct
{
    ae_int_t nfev;
    ae_int_t terminationtype;
} x_odesolverreport;
void x_set_odesolverreport(x_odesolverreport *dst, odesolverreport *src, ae_state *_state)
{
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
void odesolverreport_init_from_x(odesolverreport *dst, x_odesolverreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
DLLEXPORT int alglib_odesolverrkck(const char **errormsg, x_vector* y, ae_int_t* n, x_vector* x, ae_int_t* m, double* eps, double* h, x_odesolverstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_odesolverstate(&_alglib_env_state);
    odesolverrkck(&_y, *n, &_x, *m, *eps, *h, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_odesolveriteration(const char **errormsg, ae_bool* result, x_odesolverstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = odesolveriteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_odesolverresults(const char **errormsg, x_odesolverstate** state, ae_int_t* m, x_vector* xtbl, x_matrix* ytbl, x_odesolverreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _xtbl;
    ae_matrix _ytbl;
    odesolverreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_xtbl, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_ytbl, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    _odesolverreport_init(&_rep, &_alglib_env_state, ae_true);
    odesolverresults(&(*state)->obj, m, &_xtbl, &_ytbl, &_rep, &_alglib_env_state);
    ae_x_set_vector(xtbl, &_xtbl, &_alglib_env_state);
    ae_x_set_matrix(ytbl, &_ytbl, &_alglib_env_state);
    x_set_odesolverreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fftc1d(const char **errormsg, x_vector* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    fftc1d(&_a, *n, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fftc1dinv(const char **errormsg, x_vector* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    fftc1dinv(&_a, *n, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fftr1d(const char **errormsg, x_vector* a, ae_int_t* n, x_vector* f)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_f, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    fftr1d(&_a, *n, &_f, &_alglib_env_state);
    ae_x_set_vector(f, &_f, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fftr1dinv(const char **errormsg, x_vector* f, ae_int_t* n, x_vector* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _f;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    ae_vector_init(&_a, 0, DT_REAL, &_alglib_env_state, ae_true);
    fftr1dinv(&_f, *n, &_a, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convc1d(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    convc1d(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convc1dinv(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    convc1dinv(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convc1dcircular(const char **errormsg, x_vector* s, ae_int_t* m, x_vector* r, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _s;
    ae_vector _r;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_s, s, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_r, r, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    convc1dcircular(&_s, *m, &_r, *n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convc1dcircularinv(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    convc1dcircularinv(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convr1d(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_REAL, &_alglib_env_state, ae_true);
    convr1d(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convr1dinv(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_REAL, &_alglib_env_state, ae_true);
    convr1dinv(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convr1dcircular(const char **errormsg, x_vector* s, ae_int_t* m, x_vector* r, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _s;
    ae_vector _r;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_s, s, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_r, r, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    convr1dcircular(&_s, *m, &_r, *n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_convr1dcircularinv(const char **errormsg, x_vector* a, ae_int_t* m, x_vector* b, ae_int_t* n, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_REAL, &_alglib_env_state, ae_true);
    convr1dcircularinv(&_a, *m, &_b, *n, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_corrc1d(const char **errormsg, x_vector* signal, ae_int_t* n, x_vector* pattern, ae_int_t* m, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _signal;
    ae_vector _pattern;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_signal, signal, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pattern, pattern, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    corrc1d(&_signal, *n, &_pattern, *m, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_corrc1dcircular(const char **errormsg, x_vector* signal, ae_int_t* m, x_vector* pattern, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _signal;
    ae_vector _pattern;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_signal, signal, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pattern, pattern, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_COMPLEX, &_alglib_env_state, ae_true);
    corrc1dcircular(&_signal, *m, &_pattern, *n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_corrr1d(const char **errormsg, x_vector* signal, ae_int_t* n, x_vector* pattern, ae_int_t* m, x_vector* r)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _signal;
    ae_vector _pattern;
    ae_vector _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_signal, signal, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pattern, pattern, &_alglib_env_state, ae_true);
    ae_vector_init(&_r, 0, DT_REAL, &_alglib_env_state, ae_true);
    corrr1d(&_signal, *n, &_pattern, *m, &_r, &_alglib_env_state);
    ae_x_set_vector(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_corrr1dcircular(const char **errormsg, x_vector* signal, ae_int_t* m, x_vector* pattern, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _signal;
    ae_vector _pattern;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_signal, signal, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pattern, pattern, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    corrr1dcircular(&_signal, *m, &_pattern, *n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fhtr1d(const char **errormsg, x_vector* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    fhtr1d(&_a, *n, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fhtr1dinv(const char **errormsg, x_vector* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    fhtr1dinv(&_a, *n, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgeneraterec(const char **errormsg, x_vector* alpha, x_vector* beta, double* mu0, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _alpha;
    ae_vector _beta;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_alpha, alpha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_beta, beta, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgeneraterec(&_alpha, &_beta, *mu0, *n, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategausslobattorec(const char **errormsg, x_vector* alpha, x_vector* beta, double* mu0, double* a, double* b, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _alpha;
    ae_vector _beta;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_alpha, alpha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_beta, beta, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategausslobattorec(&_alpha, &_beta, *mu0, *a, *b, *n, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategaussradaurec(const char **errormsg, x_vector* alpha, x_vector* beta, double* mu0, double* a, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _alpha;
    ae_vector _beta;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_alpha, alpha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_beta, beta, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategaussradaurec(&_alpha, &_beta, *mu0, *a, *n, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategausslegendre(const char **errormsg, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategausslegendre(*n, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategaussjacobi(const char **errormsg, ae_int_t* n, double* alpha, double* beta, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategaussjacobi(*n, *alpha, *beta, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategausslaguerre(const char **errormsg, ae_int_t* n, double* alpha, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategausslaguerre(*n, *alpha, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gqgenerategausshermite(const char **errormsg, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    gqgenerategausshermite(*n, info, &_x, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gkqgeneraterec(const char **errormsg, x_vector* alpha, x_vector* beta, double* mu0, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* wkronrod, x_vector* wgauss)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _alpha;
    ae_vector _beta;
    ae_vector _x;
    ae_vector _wkronrod;
    ae_vector _wgauss;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_alpha, alpha, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_beta, beta, &_alglib_env_state, ae_true);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wkronrod, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wgauss, 0, DT_REAL, &_alglib_env_state, ae_true);
    gkqgeneraterec(&_alpha, &_beta, *mu0, *n, info, &_x, &_wkronrod, &_wgauss, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(wkronrod, &_wkronrod, &_alglib_env_state);
    ae_x_set_vector(wgauss, &_wgauss, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gkqgenerategausslegendre(const char **errormsg, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* wkronrod, x_vector* wgauss)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _wkronrod;
    ae_vector _wgauss;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wkronrod, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wgauss, 0, DT_REAL, &_alglib_env_state, ae_true);
    gkqgenerategausslegendre(*n, info, &_x, &_wkronrod, &_wgauss, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(wkronrod, &_wkronrod, &_alglib_env_state);
    ae_x_set_vector(wgauss, &_wgauss, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gkqgenerategaussjacobi(const char **errormsg, ae_int_t* n, double* alpha, double* beta, ae_int_t* info, x_vector* x, x_vector* wkronrod, x_vector* wgauss)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _wkronrod;
    ae_vector _wgauss;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wkronrod, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wgauss, 0, DT_REAL, &_alglib_env_state, ae_true);
    gkqgenerategaussjacobi(*n, *alpha, *beta, info, &_x, &_wkronrod, &_wgauss, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(wkronrod, &_wkronrod, &_alglib_env_state);
    ae_x_set_vector(wgauss, &_wgauss, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gkqlegendrecalc(const char **errormsg, ae_int_t* n, ae_int_t* info, x_vector* x, x_vector* wkronrod, x_vector* wgauss)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _wkronrod;
    ae_vector _wgauss;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wkronrod, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wgauss, 0, DT_REAL, &_alglib_env_state, ae_true);
    gkqlegendrecalc(*n, info, &_x, &_wkronrod, &_wgauss, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(wkronrod, &_wkronrod, &_alglib_env_state);
    ae_x_set_vector(wgauss, &_wgauss, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_gkqlegendretbl(const char **errormsg, ae_int_t* n, x_vector* x, x_vector* wkronrod, x_vector* wgauss, double* eps)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _wkronrod;
    ae_vector _wgauss;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wkronrod, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_wgauss, 0, DT_REAL, &_alglib_env_state, ae_true);
    gkqlegendretbl(*n, &_x, &_wkronrod, &_wgauss, eps, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(wkronrod, &_wkronrod, &_alglib_env_state);
    ae_x_set_vector(wgauss, &_wgauss, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    ae_int_t terminationtype;
    ae_int_t nfev;
    ae_int_t nintervals;
} x_autogkreport;
void x_set_autogkreport(x_autogkreport *dst, autogkreport *src, ae_state *_state)
{
    dst->terminationtype = src->terminationtype;
    dst->nfev = src->nfev;
    dst->nintervals = src->nintervals;
}
void autogkreport_init_from_x(autogkreport *dst, x_autogkreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->terminationtype = src->terminationtype;
    dst->nfev = src->nfev;
    dst->nintervals = src->nintervals;
}
typedef struct
{
    autogkstate obj;
} x_autogkstate;
x_autogkstate* x_obj_alloc_autogkstate(ae_state *_state)
{
    x_autogkstate *result;
    result = ae_malloc(sizeof(x_autogkstate), _state);
    _autogkstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_autogkstate(x_autogkstate *obj)
{
    if( obj==NULL )
        return;
    _autogkstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_autogkstate_get_needf(x_autogkstate *obj, ae_bool *result)
{
    *result = obj->obj.needf;
}
DLLEXPORT void x_autogkstate_set_needf(x_autogkstate *obj, ae_bool *result)
{
    obj->obj.needf = *result;
}
DLLEXPORT void x_autogkstate_get_x(x_autogkstate *obj, double *result)
{
    *result = obj->obj.x;
}
DLLEXPORT void x_autogkstate_set_x(x_autogkstate *obj, double *result)
{
    obj->obj.x = *result;
}
DLLEXPORT void x_autogkstate_get_xminusa(x_autogkstate *obj, double *result)
{
    *result = obj->obj.xminusa;
}
DLLEXPORT void x_autogkstate_set_xminusa(x_autogkstate *obj, double *result)
{
    obj->obj.xminusa = *result;
}
DLLEXPORT void x_autogkstate_get_bminusx(x_autogkstate *obj, double *result)
{
    *result = obj->obj.bminusx;
}
DLLEXPORT void x_autogkstate_set_bminusx(x_autogkstate *obj, double *result)
{
    obj->obj.bminusx = *result;
}
DLLEXPORT void x_autogkstate_get_f(x_autogkstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_autogkstate_set_f(x_autogkstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT int alglib_autogksmooth(const char **errormsg, double* a, double* b, x_autogkstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *state = x_obj_alloc_autogkstate(&_alglib_env_state);
    autogksmooth(*a, *b, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_autogksmoothw(const char **errormsg, double* a, double* b, double* xwidth, x_autogkstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *state = x_obj_alloc_autogkstate(&_alglib_env_state);
    autogksmoothw(*a, *b, *xwidth, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_autogksingular(const char **errormsg, double* a, double* b, double* alpha, double* beta, x_autogkstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *state = x_obj_alloc_autogkstate(&_alglib_env_state);
    autogksingular(*a, *b, *alpha, *beta, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_autogkiteration(const char **errormsg, ae_bool* result, x_autogkstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = autogkiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_autogkresults(const char **errormsg, x_autogkstate** state, double* v, x_autogkreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    autogkreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    _autogkreport_init(&_rep, &_alglib_env_state, ae_true);
    autogkresults(&(*state)->obj, v, &_rep, &_alglib_env_state);
    x_set_autogkreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    idwinterpolant obj;
} x_idwinterpolant;
x_idwinterpolant* x_obj_alloc_idwinterpolant(ae_state *_state)
{
    x_idwinterpolant *result;
    result = ae_malloc(sizeof(x_idwinterpolant), _state);
    _idwinterpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_idwinterpolant(x_idwinterpolant *obj)
{
    if( obj==NULL )
        return;
    _idwinterpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_idwcalc(const char **errormsg, double* result, x_idwinterpolant** z, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *result = idwcalc(&(*z)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_idwbuildmodifiedshepard(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* nx, ae_int_t* d, ae_int_t* nq, ae_int_t* nw, x_idwinterpolant** z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *z = x_obj_alloc_idwinterpolant(&_alglib_env_state);
    idwbuildmodifiedshepard(&_xy, *n, *nx, *d, *nq, *nw, &(*z)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_idwbuildmodifiedshepardr(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* nx, double* r, x_idwinterpolant** z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *z = x_obj_alloc_idwinterpolant(&_alglib_env_state);
    idwbuildmodifiedshepardr(&_xy, *n, *nx, *r, &(*z)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_idwbuildnoisy(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* nx, ae_int_t* d, ae_int_t* nq, ae_int_t* nw, x_idwinterpolant** z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *z = x_obj_alloc_idwinterpolant(&_alglib_env_state);
    idwbuildnoisy(&_xy, *n, *nx, *d, *nq, *nw, &(*z)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    barycentricinterpolant obj;
} x_barycentricinterpolant;
x_barycentricinterpolant* x_obj_alloc_barycentricinterpolant(ae_state *_state)
{
    x_barycentricinterpolant *result;
    result = ae_malloc(sizeof(x_barycentricinterpolant), _state);
    _barycentricinterpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_barycentricinterpolant(x_barycentricinterpolant *obj)
{
    if( obj==NULL )
        return;
    _barycentricinterpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_barycentriccalc(const char **errormsg, double* result, x_barycentricinterpolant** b, double* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = barycentriccalc(&(*b)->obj, *t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricdiff1(const char **errormsg, x_barycentricinterpolant** b, double* t, double* f, double* df)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    barycentricdiff1(&(*b)->obj, *t, f, df, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricdiff2(const char **errormsg, x_barycentricinterpolant** b, double* t, double* f, double* df, double* d2f)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    barycentricdiff2(&(*b)->obj, *t, f, df, d2f, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentriclintransx(const char **errormsg, x_barycentricinterpolant** b, double* ca, double* cb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    barycentriclintransx(&(*b)->obj, *ca, *cb, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentriclintransy(const char **errormsg, x_barycentricinterpolant** b, double* ca, double* cb)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    barycentriclintransy(&(*b)->obj, *ca, *cb, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricunpack(const char **errormsg, x_barycentricinterpolant** b, ae_int_t* n, x_vector* x, x_vector* y, x_vector* w)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_y, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_w, 0, DT_REAL, &_alglib_env_state, ae_true);
    barycentricunpack(&(*b)->obj, n, &_x, &_y, &_w, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    ae_x_set_vector(y, &_y, &_alglib_env_state);
    ae_x_set_vector(w, &_w, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricbuildxyw(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, x_barycentricinterpolant** b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    *b = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    barycentricbuildxyw(&_x, &_y, &_w, *n, &(*b)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricbuildfloaterhormann(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* d, x_barycentricinterpolant** b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *b = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    barycentricbuildfloaterhormann(&_x, &_y, *n, *d, &(*b)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbar2cheb(const char **errormsg, x_barycentricinterpolant** p, double* a, double* b, x_vector* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _t;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_t, 0, DT_REAL, &_alglib_env_state, ae_true);
    polynomialbar2cheb(&(*p)->obj, *a, *b, &_t, &_alglib_env_state);
    ae_x_set_vector(t, &_t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialcheb2bar(const char **errormsg, x_vector* t, ae_int_t* n, double* a, double* b, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _t;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_t, t, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialcheb2bar(&_t, *n, *a, *b, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbar2pow(const char **errormsg, x_barycentricinterpolant** p, double* c, double* s, x_vector* a)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_a, 0, DT_REAL, &_alglib_env_state, ae_true);
    polynomialbar2pow(&(*p)->obj, *c, *s, &_a, &_alglib_env_state);
    ae_x_set_vector(a, &_a, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialpow2bar(const char **errormsg, x_vector* a, ae_int_t* n, double* c, double* s, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialpow2bar(&_a, *n, *c, *s, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbuild(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialbuild(&_x, &_y, *n, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbuildeqdist(const char **errormsg, double* a, double* b, x_vector* y, ae_int_t* n, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialbuildeqdist(*a, *b, &_y, *n, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbuildcheb1(const char **errormsg, double* a, double* b, x_vector* y, ae_int_t* n, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialbuildcheb1(*a, *b, &_y, *n, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialbuildcheb2(const char **errormsg, double* a, double* b, x_vector* y, ae_int_t* n, x_barycentricinterpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    polynomialbuildcheb2(*a, *b, &_y, *n, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialcalceqdist(const char **errormsg, double* result, double* a, double* b, x_vector* f, ae_int_t* n, double* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    *result = polynomialcalceqdist(*a, *b, &_f, *n, *t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialcalccheb1(const char **errormsg, double* result, double* a, double* b, x_vector* f, ae_int_t* n, double* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    *result = polynomialcalccheb1(*a, *b, &_f, *n, *t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialcalccheb2(const char **errormsg, double* result, double* a, double* b, x_vector* f, ae_int_t* n, double* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    *result = polynomialcalccheb2(*a, *b, &_f, *n, *t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    spline1dinterpolant obj;
} x_spline1dinterpolant;
x_spline1dinterpolant* x_obj_alloc_spline1dinterpolant(ae_state *_state)
{
    x_spline1dinterpolant *result;
    result = ae_malloc(sizeof(x_spline1dinterpolant), _state);
    _spline1dinterpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_spline1dinterpolant(x_spline1dinterpolant *obj)
{
    if( obj==NULL )
        return;
    _spline1dinterpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_spline1dbuildlinear(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, x_spline1dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    spline1dbuildlinear(&_x, &_y, *n, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dbuildcubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_spline1dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    spline1dbuildcubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dgriddiffcubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_vector* d)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _d;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dgriddiffcubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &_d, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dgriddiff2cubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_vector* d1, x_vector* d2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _d1;
    ae_vector _d2;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init(&_d1, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_d2, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dgriddiff2cubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &_d1, &_d2, &_alglib_env_state);
    ae_x_set_vector(d1, &_d1, &_alglib_env_state);
    ae_x_set_vector(d2, &_d2, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dconvcubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_vector* x2, ae_int_t* n2, x_vector* y2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _x2;
    ae_vector _y2;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x2, x2, &_alglib_env_state, ae_true);
    ae_vector_init(&_y2, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dconvcubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &_x2, *n2, &_y2, &_alglib_env_state);
    ae_x_set_vector(y2, &_y2, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dconvdiffcubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_vector* x2, ae_int_t* n2, x_vector* y2, x_vector* d2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _x2;
    ae_vector _y2;
    ae_vector _d2;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x2, x2, &_alglib_env_state, ae_true);
    ae_vector_init(&_y2, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_d2, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dconvdiffcubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &_x2, *n2, &_y2, &_d2, &_alglib_env_state);
    ae_x_set_vector(y2, &_y2, &_alglib_env_state);
    ae_x_set_vector(d2, &_d2, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dconvdiff2cubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundltype, double* boundl, ae_int_t* boundrtype, double* boundr, x_vector* x2, ae_int_t* n2, x_vector* y2, x_vector* d2, x_vector* dd2)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _x2;
    ae_vector _y2;
    ae_vector _d2;
    ae_vector _dd2;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_x2, x2, &_alglib_env_state, ae_true);
    ae_vector_init(&_y2, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_d2, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_vector_init(&_dd2, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dconvdiff2cubic(&_x, &_y, *n, *boundltype, *boundl, *boundrtype, *boundr, &_x2, *n2, &_y2, &_d2, &_dd2, &_alglib_env_state);
    ae_x_set_vector(y2, &_y2, &_alglib_env_state);
    ae_x_set_vector(d2, &_d2, &_alglib_env_state);
    ae_x_set_vector(dd2, &_dd2, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dbuildcatmullrom(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* boundtype, double* tension, x_spline1dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    spline1dbuildcatmullrom(&_x, &_y, *n, *boundtype, *tension, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dbuildhermite(const char **errormsg, x_vector* x, x_vector* y, x_vector* d, ae_int_t* n, x_spline1dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _d;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_d, d, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    spline1dbuildhermite(&_x, &_y, &_d, *n, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dbuildakima(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, x_spline1dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    spline1dbuildakima(&_x, &_y, *n, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dcalc(const char **errormsg, double* result, x_spline1dinterpolant** c, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = spline1dcalc(&(*c)->obj, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1ddiff(const char **errormsg, x_spline1dinterpolant** c, double* x, double* s, double* ds, double* d2s)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline1ddiff(&(*c)->obj, *x, s, ds, d2s, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dunpack(const char **errormsg, x_spline1dinterpolant** c, ae_int_t* n, x_matrix* tbl)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _tbl;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_tbl, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline1dunpack(&(*c)->obj, n, &_tbl, &_alglib_env_state);
    ae_x_set_matrix(tbl, &_tbl, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dlintransx(const char **errormsg, x_spline1dinterpolant** c, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline1dlintransx(&(*c)->obj, *a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dlintransy(const char **errormsg, x_spline1dinterpolant** c, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline1dlintransy(&(*c)->obj, *a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dintegrate(const char **errormsg, double* result, x_spline1dinterpolant** c, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = spline1dintegrate(&(*c)->obj, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    minlmstate obj;
} x_minlmstate;
x_minlmstate* x_obj_alloc_minlmstate(ae_state *_state)
{
    x_minlmstate *result;
    result = ae_malloc(sizeof(x_minlmstate), _state);
    _minlmstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_minlmstate(x_minlmstate *obj)
{
    if( obj==NULL )
        return;
    _minlmstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_minlmstate_get_needf(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.needf;
}
DLLEXPORT void x_minlmstate_set_needf(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.needf = *result;
}
DLLEXPORT void x_minlmstate_get_needfg(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_minlmstate_set_needfg(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_minlmstate_get_needfgh(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.needfgh;
}
DLLEXPORT void x_minlmstate_set_needfgh(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.needfgh = *result;
}
DLLEXPORT void x_minlmstate_get_needfi(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.needfi;
}
DLLEXPORT void x_minlmstate_set_needfi(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.needfi = *result;
}
DLLEXPORT void x_minlmstate_get_needfij(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.needfij;
}
DLLEXPORT void x_minlmstate_set_needfij(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.needfij = *result;
}
DLLEXPORT void x_minlmstate_get_xupdated(x_minlmstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_minlmstate_set_xupdated(x_minlmstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_minlmstate_get_f(x_minlmstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_minlmstate_set_f(x_minlmstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_minlmstate_get_fi(x_minlmstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.fi);
}
DLLEXPORT void x_minlmstate_get_g(x_minlmstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_minlmstate_get_h(x_minlmstate *obj, x_matrix *result)
{
    ae_x_attach_to_matrix(result, &obj->obj.h);
}
DLLEXPORT void x_minlmstate_get_j(x_minlmstate *obj, x_matrix *result)
{
    ae_x_attach_to_matrix(result, &obj->obj.j);
}
DLLEXPORT void x_minlmstate_get_x(x_minlmstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t terminationtype;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
} x_minlmreport;
void x_set_minlmreport(x_minlmreport *dst, minlmreport *src, ae_state *_state)
{
    dst->iterationscount = src->iterationscount;
    dst->terminationtype = src->terminationtype;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
}
void minlmreport_init_from_x(minlmreport *dst, x_minlmreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->terminationtype = src->terminationtype;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
}
DLLEXPORT int alglib_minlmcreatevj(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatevj(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmcreatev(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, double* diffstep, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatev(*n, *m, &_x, *diffstep, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmcreatefgh(const char **errormsg, ae_int_t* n, x_vector* x, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatefgh(*n, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmcreatevgj(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatevgj(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmcreatefgj(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatefgj(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmcreatefj(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minlmstate(&_alglib_env_state);
    minlmcreatefj(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmsetcond(const char **errormsg, x_minlmstate** state, double* epsg, double* epsf, double* epsx, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlmsetcond(&(*state)->obj, *epsg, *epsf, *epsx, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmsetxrep(const char **errormsg, x_minlmstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlmsetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmsetstpmax(const char **errormsg, x_minlmstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlmsetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmsetacctype(const char **errormsg, x_minlmstate** state, ae_int_t* acctype)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minlmsetacctype(&(*state)->obj, *acctype, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmiteration(const char **errormsg, ae_bool* result, x_minlmstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = minlmiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmresults(const char **errormsg, x_minlmstate** state, x_vector* x, x_minlmreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minlmreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _minlmreport_init(&_rep, &_alglib_env_state, ae_true);
    minlmresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minlmreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmresultsbuf(const char **errormsg, x_minlmstate** state, x_vector* x, x_minlmreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minlmreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minlmreport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    minlmresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minlmreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minlmrestartfrom(const char **errormsg, x_minlmstate** state, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minlmrestartfrom(&(*state)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} x_polynomialfitreport;
void x_set_polynomialfitreport(x_polynomialfitreport *dst, polynomialfitreport *src, ae_state *_state)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
void polynomialfitreport_init_from_x(polynomialfitreport *dst, x_polynomialfitreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
typedef struct
{
    double taskrcond;
    ae_int_t dbest;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} x_barycentricfitreport;
void x_set_barycentricfitreport(x_barycentricfitreport *dst, barycentricfitreport *src, ae_state *_state)
{
    dst->taskrcond = src->taskrcond;
    dst->dbest = src->dbest;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
void barycentricfitreport_init_from_x(barycentricfitreport *dst, x_barycentricfitreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->dbest = src->dbest;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} x_spline1dfitreport;
void x_set_spline1dfitreport(x_spline1dfitreport *dst, spline1dfitreport *src, ae_state *_state)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
void spline1dfitreport_init_from_x(spline1dfitreport *dst, x_spline1dfitreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
typedef struct
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
} x_lsfitreport;
void x_set_lsfitreport(x_lsfitreport *dst, lsfitreport *src, ae_state *_state)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
void lsfitreport_init_from_x(lsfitreport *dst, x_lsfitreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->taskrcond = src->taskrcond;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->maxerror = src->maxerror;
}
typedef struct
{
    lsfitstate obj;
} x_lsfitstate;
x_lsfitstate* x_obj_alloc_lsfitstate(ae_state *_state)
{
    x_lsfitstate *result;
    result = ae_malloc(sizeof(x_lsfitstate), _state);
    _lsfitstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_lsfitstate(x_lsfitstate *obj)
{
    if( obj==NULL )
        return;
    _lsfitstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_lsfitstate_get_needf(x_lsfitstate *obj, ae_bool *result)
{
    *result = obj->obj.needf;
}
DLLEXPORT void x_lsfitstate_set_needf(x_lsfitstate *obj, ae_bool *result)
{
    obj->obj.needf = *result;
}
DLLEXPORT void x_lsfitstate_get_needfg(x_lsfitstate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_lsfitstate_set_needfg(x_lsfitstate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_lsfitstate_get_needfgh(x_lsfitstate *obj, ae_bool *result)
{
    *result = obj->obj.needfgh;
}
DLLEXPORT void x_lsfitstate_set_needfgh(x_lsfitstate *obj, ae_bool *result)
{
    obj->obj.needfgh = *result;
}
DLLEXPORT void x_lsfitstate_get_xupdated(x_lsfitstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_lsfitstate_set_xupdated(x_lsfitstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_lsfitstate_get_c(x_lsfitstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.c);
}
DLLEXPORT void x_lsfitstate_get_f(x_lsfitstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_lsfitstate_set_f(x_lsfitstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_lsfitstate_get_g(x_lsfitstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_lsfitstate_get_h(x_lsfitstate *obj, x_matrix *result)
{
    ae_x_attach_to_matrix(result, &obj->obj.h);
}
DLLEXPORT void x_lsfitstate_get_x(x_lsfitstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
DLLEXPORT int alglib_polynomialfit(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_barycentricinterpolant** p, x_polynomialfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    polynomialfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    _polynomialfitreport_init(&_rep, &_alglib_env_state, ae_true);
    polynomialfit(&_x, &_y, *n, *m, info, &(*p)->obj, &_rep, &_alglib_env_state);
    x_set_polynomialfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_polynomialfitwc(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, x_vector* xc, x_vector* yc, x_vector* dc, ae_int_t* k, ae_int_t* m, ae_int_t* info, x_barycentricinterpolant** p, x_polynomialfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_vector _dc;
    polynomialfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_xc, xc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_yc, yc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_dc, dc, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    _polynomialfitreport_init(&_rep, &_alglib_env_state, ae_true);
    polynomialfitwc(&_x, &_y, &_w, *n, &_xc, &_yc, &_dc, *k, *m, info, &(*p)->obj, &_rep, &_alglib_env_state);
    x_set_polynomialfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricfitfloaterhormannwc(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, x_vector* xc, x_vector* yc, x_vector* dc, ae_int_t* k, ae_int_t* m, ae_int_t* info, x_barycentricinterpolant** b, x_barycentricfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_vector _dc;
    barycentricfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_xc, xc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_yc, yc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_dc, dc, &_alglib_env_state, ae_true);
    *b = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    _barycentricfitreport_init(&_rep, &_alglib_env_state, ae_true);
    barycentricfitfloaterhormannwc(&_x, &_y, &_w, *n, &_xc, &_yc, &_dc, *k, *m, info, &(*b)->obj, &_rep, &_alglib_env_state);
    x_set_barycentricfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_barycentricfitfloaterhormann(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_barycentricinterpolant** b, x_barycentricfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    barycentricfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *b = x_obj_alloc_barycentricinterpolant(&_alglib_env_state);
    _barycentricfitreport_init(&_rep, &_alglib_env_state, ae_true);
    barycentricfitfloaterhormann(&_x, &_y, *n, *m, info, &(*b)->obj, &_rep, &_alglib_env_state);
    x_set_barycentricfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfitpenalized(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* m, double* rho, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfitpenalized(&_x, &_y, *n, *m, *rho, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfitpenalizedw(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, ae_int_t* m, double* rho, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfitpenalizedw(&_x, &_y, &_w, *n, *m, *rho, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfitcubicwc(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, x_vector* xc, x_vector* yc, x_vector* dc, ae_int_t* k, ae_int_t* m, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_vector _dc;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_xc, xc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_yc, yc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_dc, dc, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfitcubicwc(&_x, &_y, &_w, *n, &_xc, &_yc, &_dc, *k, *m, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfithermitewc(const char **errormsg, x_vector* x, x_vector* y, x_vector* w, ae_int_t* n, x_vector* xc, x_vector* yc, x_vector* dc, ae_int_t* k, ae_int_t* m, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _xc;
    ae_vector _yc;
    ae_vector _dc;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_xc, xc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_yc, yc, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_dc, dc, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfithermitewc(&_x, &_y, &_w, *n, &_xc, &_yc, &_dc, *k, *m, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfitcubic(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfitcubic(&_x, &_y, *n, *m, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline1dfithermite(const char **errormsg, x_vector* x, x_vector* y, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_spline1dinterpolant** s, x_spline1dfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    spline1dfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    *s = x_obj_alloc_spline1dinterpolant(&_alglib_env_state);
    _spline1dfitreport_init(&_rep, &_alglib_env_state, ae_true);
    spline1dfithermite(&_x, &_y, *n, *m, info, &(*s)->obj, &_rep, &_alglib_env_state);
    x_set_spline1dfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitlinearw(const char **errormsg, x_vector* y, x_vector* w, x_matrix* fmatrix, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_vector* c, x_lsfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_vector _w;
    ae_matrix _fmatrix;
    ae_vector _c;
    lsfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_fmatrix, fmatrix, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    _lsfitreport_init(&_rep, &_alglib_env_state, ae_true);
    lsfitlinearw(&_y, &_w, &_fmatrix, *n, *m, info, &_c, &_rep, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    x_set_lsfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitlinearwc(const char **errormsg, x_vector* y, x_vector* w, x_matrix* fmatrix, x_matrix* cmatrix, ae_int_t* n, ae_int_t* m, ae_int_t* k, ae_int_t* info, x_vector* c, x_lsfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_vector _w;
    ae_matrix _fmatrix;
    ae_matrix _cmatrix;
    ae_vector _c;
    lsfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_fmatrix, fmatrix, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_cmatrix, cmatrix, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    _lsfitreport_init(&_rep, &_alglib_env_state, ae_true);
    lsfitlinearwc(&_y, &_w, &_fmatrix, &_cmatrix, *n, *m, *k, info, &_c, &_rep, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    x_set_lsfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitlinear(const char **errormsg, x_vector* y, x_matrix* fmatrix, ae_int_t* n, ae_int_t* m, ae_int_t* info, x_vector* c, x_lsfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_matrix _fmatrix;
    ae_vector _c;
    lsfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_fmatrix, fmatrix, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    _lsfitreport_init(&_rep, &_alglib_env_state, ae_true);
    lsfitlinear(&_y, &_fmatrix, *n, *m, info, &_c, &_rep, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    x_set_lsfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitlinearc(const char **errormsg, x_vector* y, x_matrix* fmatrix, x_matrix* cmatrix, ae_int_t* n, ae_int_t* m, ae_int_t* k, ae_int_t* info, x_vector* c, x_lsfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _y;
    ae_matrix _fmatrix;
    ae_matrix _cmatrix;
    ae_vector _c;
    lsfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_fmatrix, fmatrix, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_cmatrix, cmatrix, &_alglib_env_state, ae_true);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    _lsfitreport_init(&_rep, &_alglib_env_state, ae_true);
    lsfitlinearc(&_y, &_fmatrix, &_cmatrix, *n, *m, *k, info, &_c, &_rep, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    x_set_lsfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatewf(const char **errormsg, x_matrix* x, x_vector* y, x_vector* w, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, double* diffstep, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatewf(&_x, &_y, &_w, &_c, *n, *m, *k, *diffstep, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatef(const char **errormsg, x_matrix* x, x_vector* y, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, double* diffstep, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatef(&_x, &_y, &_c, *n, *m, *k, *diffstep, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatewfg(const char **errormsg, x_matrix* x, x_vector* y, x_vector* w, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, ae_bool* cheapfg, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatewfg(&_x, &_y, &_w, &_c, *n, *m, *k, *cheapfg, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatefg(const char **errormsg, x_matrix* x, x_vector* y, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, ae_bool* cheapfg, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatefg(&_x, &_y, &_c, *n, *m, *k, *cheapfg, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatewfgh(const char **errormsg, x_matrix* x, x_vector* y, x_vector* w, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _w;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_w, w, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatewfgh(&_x, &_y, &_w, &_c, *n, *m, *k, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitcreatefgh(const char **errormsg, x_matrix* x, x_vector* y, x_vector* c, ae_int_t* n, ae_int_t* m, ae_int_t* k, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _x;
    ae_vector _y;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_lsfitstate(&_alglib_env_state);
    lsfitcreatefgh(&_x, &_y, &_c, *n, *m, *k, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitsetcond(const char **errormsg, x_lsfitstate** state, double* epsf, double* epsx, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    lsfitsetcond(&(*state)->obj, *epsf, *epsx, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitsetstpmax(const char **errormsg, x_lsfitstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    lsfitsetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitsetxrep(const char **errormsg, x_lsfitstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    lsfitsetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfititeration(const char **errormsg, ae_bool* result, x_lsfitstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = lsfititeration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_lsfitresults(const char **errormsg, x_lsfitstate** state, ae_int_t* info, x_vector* c, x_lsfitreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    lsfitreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    _lsfitreport_init(&_rep, &_alglib_env_state, ae_true);
    lsfitresults(&(*state)->obj, info, &_c, &_rep, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    x_set_lsfitreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    pspline2interpolant obj;
} x_pspline2interpolant;
x_pspline2interpolant* x_obj_alloc_pspline2interpolant(ae_state *_state)
{
    x_pspline2interpolant *result;
    result = ae_malloc(sizeof(x_pspline2interpolant), _state);
    _pspline2interpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_pspline2interpolant(x_pspline2interpolant *obj)
{
    if( obj==NULL )
        return;
    _pspline2interpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
typedef struct
{
    pspline3interpolant obj;
} x_pspline3interpolant;
x_pspline3interpolant* x_obj_alloc_pspline3interpolant(ae_state *_state)
{
    x_pspline3interpolant *result;
    result = ae_malloc(sizeof(x_pspline3interpolant), _state);
    _pspline3interpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_pspline3interpolant(x_pspline3interpolant *obj)
{
    if( obj==NULL )
        return;
    _pspline3interpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_pspline2build(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* st, ae_int_t* pt, x_pspline2interpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_pspline2interpolant(&_alglib_env_state);
    pspline2build(&_xy, *n, *st, *pt, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3build(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* st, ae_int_t* pt, x_pspline3interpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_pspline3interpolant(&_alglib_env_state);
    pspline3build(&_xy, *n, *st, *pt, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2buildperiodic(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* st, ae_int_t* pt, x_pspline2interpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_pspline2interpolant(&_alglib_env_state);
    pspline2buildperiodic(&_xy, *n, *st, *pt, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3buildperiodic(const char **errormsg, x_matrix* xy, ae_int_t* n, ae_int_t* st, ae_int_t* pt, x_pspline3interpolant** p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _xy;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_xy, xy, &_alglib_env_state, ae_true);
    *p = x_obj_alloc_pspline3interpolant(&_alglib_env_state);
    pspline3buildperiodic(&_xy, *n, *st, *pt, &(*p)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2parametervalues(const char **errormsg, x_pspline2interpolant** p, ae_int_t* n, x_vector* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _t;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_t, 0, DT_REAL, &_alglib_env_state, ae_true);
    pspline2parametervalues(&(*p)->obj, n, &_t, &_alglib_env_state);
    ae_x_set_vector(t, &_t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3parametervalues(const char **errormsg, x_pspline3interpolant** p, ae_int_t* n, x_vector* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _t;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_t, 0, DT_REAL, &_alglib_env_state, ae_true);
    pspline3parametervalues(&(*p)->obj, n, &_t, &_alglib_env_state);
    ae_x_set_vector(t, &_t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2calc(const char **errormsg, x_pspline2interpolant** p, double* t, double* x, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline2calc(&(*p)->obj, *t, x, y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3calc(const char **errormsg, x_pspline3interpolant** p, double* t, double* x, double* y, double* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline3calc(&(*p)->obj, *t, x, y, z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2tangent(const char **errormsg, x_pspline2interpolant** p, double* t, double* x, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline2tangent(&(*p)->obj, *t, x, y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3tangent(const char **errormsg, x_pspline3interpolant** p, double* t, double* x, double* y, double* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline3tangent(&(*p)->obj, *t, x, y, z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2diff(const char **errormsg, x_pspline2interpolant** p, double* t, double* x, double* dx, double* y, double* dy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline2diff(&(*p)->obj, *t, x, dx, y, dy, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3diff(const char **errormsg, x_pspline3interpolant** p, double* t, double* x, double* dx, double* y, double* dy, double* z, double* dz)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline3diff(&(*p)->obj, *t, x, dx, y, dy, z, dz, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2diff2(const char **errormsg, x_pspline2interpolant** p, double* t, double* x, double* dx, double* d2x, double* y, double* dy, double* d2y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline2diff2(&(*p)->obj, *t, x, dx, d2x, y, dy, d2y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3diff2(const char **errormsg, x_pspline3interpolant** p, double* t, double* x, double* dx, double* d2x, double* y, double* dy, double* d2y, double* z, double* dz, double* d2z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pspline3diff2(&(*p)->obj, *t, x, dx, d2x, y, dy, d2y, z, dz, d2z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline2arclength(const char **errormsg, double* result, x_pspline2interpolant** p, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = pspline2arclength(&(*p)->obj, *a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pspline3arclength(const char **errormsg, double* result, x_pspline3interpolant** p, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = pspline3arclength(&(*p)->obj, *a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    spline2dinterpolant obj;
} x_spline2dinterpolant;
x_spline2dinterpolant* x_obj_alloc_spline2dinterpolant(ae_state *_state)
{
    x_spline2dinterpolant *result;
    result = ae_malloc(sizeof(x_spline2dinterpolant), _state);
    _spline2dinterpolant_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_spline2dinterpolant(x_spline2dinterpolant *obj)
{
    if( obj==NULL )
        return;
    _spline2dinterpolant_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT int alglib_spline2dbuildbilinear(const char **errormsg, x_vector* x, x_vector* y, x_matrix* f, ae_int_t* m, ae_int_t* n, x_spline2dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_matrix _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline2dinterpolant(&_alglib_env_state);
    spline2dbuildbilinear(&_x, &_y, &_f, *m, *n, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dbuildbicubic(const char **errormsg, x_vector* x, x_vector* y, x_matrix* f, ae_int_t* m, ae_int_t* n, x_spline2dinterpolant** c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_matrix _f;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_f, f, &_alglib_env_state, ae_true);
    *c = x_obj_alloc_spline2dinterpolant(&_alglib_env_state);
    spline2dbuildbicubic(&_x, &_y, &_f, *m, *n, &(*c)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dcalc(const char **errormsg, double* result, x_spline2dinterpolant** c, double* x, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = spline2dcalc(&(*c)->obj, *x, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2ddiff(const char **errormsg, x_spline2dinterpolant** c, double* x, double* y, double* f, double* fx, double* fy, double* fxy)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline2ddiff(&(*c)->obj, *x, *y, f, fx, fy, fxy, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dunpack(const char **errormsg, x_spline2dinterpolant** c, ae_int_t* m, ae_int_t* n, x_matrix* tbl)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _tbl;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init(&_tbl, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline2dunpack(&(*c)->obj, m, n, &_tbl, &_alglib_env_state);
    ae_x_set_matrix(tbl, &_tbl, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dlintransxy(const char **errormsg, x_spline2dinterpolant** c, double* ax, double* bx, double* ay, double* by)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline2dlintransxy(&(*c)->obj, *ax, *bx, *ay, *by, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dlintransf(const char **errormsg, x_spline2dinterpolant** c, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spline2dlintransf(&(*c)->obj, *a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dresamplebicubic(const char **errormsg, x_matrix* a, ae_int_t* oldheight, ae_int_t* oldwidth, x_matrix* b, ae_int_t* newheight, ae_int_t* newwidth)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_b, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline2dresamplebicubic(&_a, *oldheight, *oldwidth, &_b, *newheight, *newwidth, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spline2dresamplebilinear(const char **errormsg, x_matrix* a, ae_int_t* oldheight, ae_int_t* oldwidth, x_matrix* b, ae_int_t* newheight, ae_int_t* newwidth)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_b, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    spline2dresamplebilinear(&_a, *oldheight, *oldwidth, &_b, *newheight, *newwidth, &_alglib_env_state);
    ae_x_set_matrix(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixludet(const char **errormsg, double* result, x_matrix* a, x_vector* pivots, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pivots, pivots, &_alglib_env_state, ae_true);
    *result = rmatrixludet(&_a, &_pivots, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixdet(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = rmatrixdet(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixludet(const char **errormsg, ae_complex* result, x_matrix* a, x_vector* pivots, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_vector _pivots;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_pivots, pivots, &_alglib_env_state, ae_true);
    *result = cmatrixludet(&_a, &_pivots, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_cmatrixdet(const char **errormsg, ae_complex* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = cmatrixdet(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixcholeskydet(const char **errormsg, double* result, x_matrix* a, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = spdmatrixcholeskydet(&_a, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spdmatrixdet(const char **errormsg, double* result, x_matrix* a, ae_int_t* n, ae_bool* isupper)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    *result = spdmatrixdet(&_a, *n, *isupper, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixgevd(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_bool* isuppera, x_matrix* b, ae_bool* isupperb, ae_int_t* zneeded, ae_int_t* problemtype, x_vector* d, x_matrix* z)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_vector _d;
    ae_matrix _z;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_vector_init(&_d, 0, DT_REAL, &_alglib_env_state, ae_true);
    ae_matrix_init(&_z, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = smatrixgevd(&_a, *n, *isuppera, &_b, *isupperb, *zneeded, *problemtype, &_d, &_z, &_alglib_env_state);
    ae_x_set_vector(d, &_d, &_alglib_env_state);
    ae_x_set_matrix(z, &_z, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_smatrixgevdreduce(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, ae_bool* isuppera, x_matrix* b, ae_bool* isupperb, ae_int_t* problemtype, x_matrix* r, ae_bool* isupperr)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _b;
    ae_matrix _r;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init_from_x(&_b, b, &_alglib_env_state, ae_true);
    ae_matrix_init(&_r, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = smatrixgevdreduce(&_a, *n, *isuppera, &_b, *isupperb, *problemtype, &_r, isupperr, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_matrix(r, &_r, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixinvupdatesimple(const char **errormsg, x_matrix* inva, ae_int_t* n, ae_int_t* updrow, ae_int_t* updcolumn, double* updval)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _inva;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_inva, inva, &_alglib_env_state, ae_true);
    rmatrixinvupdatesimple(&_inva, *n, *updrow, *updcolumn, *updval, &_alglib_env_state);
    ae_x_set_matrix(inva, &_inva, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixinvupdaterow(const char **errormsg, x_matrix* inva, ae_int_t* n, ae_int_t* updrow, x_vector* v)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _inva;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_inva, inva, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_v, v, &_alglib_env_state, ae_true);
    rmatrixinvupdaterow(&_inva, *n, *updrow, &_v, &_alglib_env_state);
    ae_x_set_matrix(inva, &_inva, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixinvupdatecolumn(const char **errormsg, x_matrix* inva, ae_int_t* n, ae_int_t* updcolumn, x_vector* u)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _inva;
    ae_vector _u;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_inva, inva, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_u, u, &_alglib_env_state, ae_true);
    rmatrixinvupdatecolumn(&_inva, *n, *updcolumn, &_u, &_alglib_env_state);
    ae_x_set_matrix(inva, &_inva, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixinvupdateuv(const char **errormsg, x_matrix* inva, ae_int_t* n, x_vector* u, x_vector* v)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _inva;
    ae_vector _u;
    ae_vector _v;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_inva, inva, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_u, u, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_v, v, &_alglib_env_state, ae_true);
    rmatrixinvupdateuv(&_inva, *n, &_u, &_v, &_alglib_env_state);
    ae_x_set_matrix(inva, &_inva, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_rmatrixschur(const char **errormsg, ae_bool* result, x_matrix* a, ae_int_t* n, x_matrix* s)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _a;
    ae_matrix _s;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_matrix_init(&_s, 0, 0, DT_REAL, &_alglib_env_state, ae_true);
    *result = rmatrixschur(&_a, *n, &_s, &_alglib_env_state);
    ae_x_set_matrix(a, &_a, &_alglib_env_state);
    ae_x_set_matrix(s, &_s, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    minasastate obj;
} x_minasastate;
x_minasastate* x_obj_alloc_minasastate(ae_state *_state)
{
    x_minasastate *result;
    result = ae_malloc(sizeof(x_minasastate), _state);
    _minasastate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_minasastate(x_minasastate *obj)
{
    if( obj==NULL )
        return;
    _minasastate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_minasastate_get_needfg(x_minasastate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_minasastate_set_needfg(x_minasastate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_minasastate_get_xupdated(x_minasastate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_minasastate_set_xupdated(x_minasastate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_minasastate_get_f(x_minasastate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_minasastate_set_f(x_minasastate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_minasastate_get_g(x_minasastate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_minasastate_get_x(x_minasastate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
    ae_int_t activeconstraints;
} x_minasareport;
void x_set_minasareport(x_minasareport *dst, minasareport *src, ae_state *_state)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->activeconstraints = src->activeconstraints;
}
void minasareport_init_from_x(minasareport *dst, x_minasareport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->activeconstraints = src->activeconstraints;
}
DLLEXPORT int alglib_minasacreate(const char **errormsg, ae_int_t* n, x_vector* x, x_vector* bndl, x_vector* bndu, x_minasastate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _bndl;
    ae_vector _bndu;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bndl, bndl, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bndu, bndu, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minasastate(&_alglib_env_state);
    minasacreate(*n, &_x, &_bndl, &_bndu, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasasetcond(const char **errormsg, x_minasastate** state, double* epsg, double* epsf, double* epsx, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minasasetcond(&(*state)->obj, *epsg, *epsf, *epsx, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasasetxrep(const char **errormsg, x_minasastate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minasasetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasasetalgorithm(const char **errormsg, x_minasastate** state, ae_int_t* algotype)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minasasetalgorithm(&(*state)->obj, *algotype, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasasetstpmax(const char **errormsg, x_minasastate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minasasetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasaiteration(const char **errormsg, ae_bool* result, x_minasastate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = minasaiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasaresults(const char **errormsg, x_minasastate** state, x_vector* x, x_minasareport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minasareport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _minasareport_init(&_rep, &_alglib_env_state, ae_true);
    minasaresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minasareport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasaresultsbuf(const char **errormsg, x_minasastate** state, x_vector* x, x_minasareport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minasareport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minasareport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    minasaresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minasareport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minasarestartfrom(const char **errormsg, x_minasastate** state, x_vector* x, x_vector* bndl, x_vector* bndu)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _bndl;
    ae_vector _bndu;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bndl, bndl, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bndu, bndu, &_alglib_env_state, ae_true);
    minasarestartfrom(&(*state)->obj, &_x, &_bndl, &_bndu, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    mincgstate obj;
} x_mincgstate;
x_mincgstate* x_obj_alloc_mincgstate(ae_state *_state)
{
    x_mincgstate *result;
    result = ae_malloc(sizeof(x_mincgstate), _state);
    _mincgstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_mincgstate(x_mincgstate *obj)
{
    if( obj==NULL )
        return;
    _mincgstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_mincgstate_get_needfg(x_mincgstate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_mincgstate_set_needfg(x_mincgstate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_mincgstate_get_xupdated(x_mincgstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_mincgstate_set_xupdated(x_mincgstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_mincgstate_get_f(x_mincgstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_mincgstate_set_f(x_mincgstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_mincgstate_get_g(x_mincgstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_mincgstate_get_x(x_mincgstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
} x_mincgreport;
void x_set_mincgreport(x_mincgreport *dst, mincgreport *src, ae_state *_state)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
void mincgreport_init_from_x(mincgreport *dst, x_mincgreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}
DLLEXPORT int alglib_mincgcreate(const char **errormsg, ae_int_t* n, x_vector* x, x_mincgstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_mincgstate(&_alglib_env_state);
    mincgcreate(*n, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgsetcond(const char **errormsg, x_mincgstate** state, double* epsg, double* epsf, double* epsx, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mincgsetcond(&(*state)->obj, *epsg, *epsf, *epsx, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgsetxrep(const char **errormsg, x_mincgstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mincgsetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgsetcgtype(const char **errormsg, x_mincgstate** state, ae_int_t* cgtype)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mincgsetcgtype(&(*state)->obj, *cgtype, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgsetstpmax(const char **errormsg, x_mincgstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mincgsetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgsuggeststep(const char **errormsg, x_mincgstate** state, double* stp)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    mincgsuggeststep(&(*state)->obj, *stp, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgiteration(const char **errormsg, ae_bool* result, x_mincgstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = mincgiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgresults(const char **errormsg, x_mincgstate** state, x_vector* x, x_mincgreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    mincgreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _mincgreport_init(&_rep, &_alglib_env_state, ae_true);
    mincgresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_mincgreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgresultsbuf(const char **errormsg, x_mincgstate** state, x_vector* x, x_mincgreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    mincgreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    mincgreport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    mincgresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_mincgreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mincgrestartfrom(const char **errormsg, x_mincgstate** state, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    mincgrestartfrom(&(*state)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    minbleicstate obj;
} x_minbleicstate;
x_minbleicstate* x_obj_alloc_minbleicstate(ae_state *_state)
{
    x_minbleicstate *result;
    result = ae_malloc(sizeof(x_minbleicstate), _state);
    _minbleicstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_minbleicstate(x_minbleicstate *obj)
{
    if( obj==NULL )
        return;
    _minbleicstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_minbleicstate_get_needfg(x_minbleicstate *obj, ae_bool *result)
{
    *result = obj->obj.needfg;
}
DLLEXPORT void x_minbleicstate_set_needfg(x_minbleicstate *obj, ae_bool *result)
{
    obj->obj.needfg = *result;
}
DLLEXPORT void x_minbleicstate_get_xupdated(x_minbleicstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_minbleicstate_set_xupdated(x_minbleicstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_minbleicstate_get_f(x_minbleicstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_minbleicstate_set_f(x_minbleicstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_minbleicstate_get_g(x_minbleicstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.g);
}
DLLEXPORT void x_minbleicstate_get_x(x_minbleicstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
    double debugeqerr;
    double debugfs;
    double debugff;
    double debugdx;
} x_minbleicreport;
void x_set_minbleicreport(x_minbleicreport *dst, minbleicreport *src, ae_state *_state)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->debugeqerr = src->debugeqerr;
    dst->debugfs = src->debugfs;
    dst->debugff = src->debugff;
    dst->debugdx = src->debugdx;
}
void minbleicreport_init_from_x(minbleicreport *dst, x_minbleicreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->debugeqerr = src->debugeqerr;
    dst->debugfs = src->debugfs;
    dst->debugff = src->debugff;
    dst->debugdx = src->debugdx;
}
DLLEXPORT int alglib_minbleiccreate(const char **errormsg, ae_int_t* n, x_vector* x, x_minbleicstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_minbleicstate(&_alglib_env_state);
    minbleiccreate(*n, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetbc(const char **errormsg, x_minbleicstate** state, x_vector* bndl, x_vector* bndu)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _bndl;
    ae_vector _bndu;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_bndl, bndl, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_bndu, bndu, &_alglib_env_state, ae_true);
    minbleicsetbc(&(*state)->obj, &_bndl, &_bndu, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetlc(const char **errormsg, x_minbleicstate** state, x_matrix* c, x_vector* ct, ae_int_t* k)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_matrix _c;
    ae_vector _ct;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_matrix_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_ct, ct, &_alglib_env_state, ae_true);
    minbleicsetlc(&(*state)->obj, &_c, &_ct, *k, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetinnercond(const char **errormsg, x_minbleicstate** state, double* epsg, double* epsf, double* epsx)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetinnercond(&(*state)->obj, *epsg, *epsf, *epsx, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetoutercond(const char **errormsg, x_minbleicstate** state, double* epsx, double* epsi)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetoutercond(&(*state)->obj, *epsx, *epsi, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetbarrierwidth(const char **errormsg, x_minbleicstate** state, double* mu)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetbarrierwidth(&(*state)->obj, *mu, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetbarrierdecay(const char **errormsg, x_minbleicstate** state, double* mudecay)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetbarrierdecay(&(*state)->obj, *mudecay, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetmaxits(const char **errormsg, x_minbleicstate** state, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetmaxits(&(*state)->obj, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetxrep(const char **errormsg, x_minbleicstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicsetstpmax(const char **errormsg, x_minbleicstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    minbleicsetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleiciteration(const char **errormsg, ae_bool* result, x_minbleicstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = minbleiciteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicresults(const char **errormsg, x_minbleicstate** state, x_vector* x, x_minbleicreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minbleicreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _minbleicreport_init(&_rep, &_alglib_env_state, ae_true);
    minbleicresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minbleicreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicresultsbuf(const char **errormsg, x_minbleicstate** state, x_vector* x, x_minbleicreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    minbleicreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minbleicreport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    minbleicresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_minbleicreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_minbleicrestartfrom(const char **errormsg, x_minbleicstate** state, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    minbleicrestartfrom(&(*state)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
typedef struct
{
    nleqstate obj;
} x_nleqstate;
x_nleqstate* x_obj_alloc_nleqstate(ae_state *_state)
{
    x_nleqstate *result;
    result = ae_malloc(sizeof(x_nleqstate), _state);
    _nleqstate_init(&result->obj, _state, ae_false);
    return result;
}
DLLEXPORT void x_obj_free_nleqstate(x_nleqstate *obj)
{
    if( obj==NULL )
        return;
    _nleqstate_clear(&obj->obj);
    ae_free(obj);
    return;
}
DLLEXPORT void x_nleqstate_get_needf(x_nleqstate *obj, ae_bool *result)
{
    *result = obj->obj.needf;
}
DLLEXPORT void x_nleqstate_set_needf(x_nleqstate *obj, ae_bool *result)
{
    obj->obj.needf = *result;
}
DLLEXPORT void x_nleqstate_get_needfij(x_nleqstate *obj, ae_bool *result)
{
    *result = obj->obj.needfij;
}
DLLEXPORT void x_nleqstate_set_needfij(x_nleqstate *obj, ae_bool *result)
{
    obj->obj.needfij = *result;
}
DLLEXPORT void x_nleqstate_get_xupdated(x_nleqstate *obj, ae_bool *result)
{
    *result = obj->obj.xupdated;
}
DLLEXPORT void x_nleqstate_set_xupdated(x_nleqstate *obj, ae_bool *result)
{
    obj->obj.xupdated = *result;
}
DLLEXPORT void x_nleqstate_get_f(x_nleqstate *obj, double *result)
{
    *result = obj->obj.f;
}
DLLEXPORT void x_nleqstate_set_f(x_nleqstate *obj, double *result)
{
    obj->obj.f = *result;
}
DLLEXPORT void x_nleqstate_get_fi(x_nleqstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.fi);
}
DLLEXPORT void x_nleqstate_get_j(x_nleqstate *obj, x_matrix *result)
{
    ae_x_attach_to_matrix(result, &obj->obj.j);
}
DLLEXPORT void x_nleqstate_get_x(x_nleqstate *obj, x_vector *result)
{
    ae_x_attach_to_vector(result, &obj->obj.x);
}
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t terminationtype;
} x_nleqreport;
void x_set_nleqreport(x_nleqreport *dst, nleqreport *src, ae_state *_state)
{
    dst->iterationscount = src->iterationscount;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->terminationtype = src->terminationtype;
}
void nleqreport_init_from_x(nleqreport *dst, x_nleqreport *src, ae_state *_state, ae_bool make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->terminationtype = src->terminationtype;
}
DLLEXPORT int alglib_nleqcreatelm(const char **errormsg, ae_int_t* n, ae_int_t* m, x_vector* x, x_nleqstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    *state = x_obj_alloc_nleqstate(&_alglib_env_state);
    nleqcreatelm(*n, *m, &_x, &(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqsetcond(const char **errormsg, x_nleqstate** state, double* epsf, ae_int_t* maxits)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    nleqsetcond(&(*state)->obj, *epsf, *maxits, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqsetxrep(const char **errormsg, x_nleqstate** state, ae_bool* needxrep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    nleqsetxrep(&(*state)->obj, *needxrep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqsetstpmax(const char **errormsg, x_nleqstate** state, double* stpmax)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    nleqsetstpmax(&(*state)->obj, *stpmax, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqiteration(const char **errormsg, ae_bool* result, x_nleqstate** state)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = nleqiteration(&(*state)->obj, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqresults(const char **errormsg, x_nleqstate** state, x_vector* x, x_nleqreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    nleqreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_x, 0, DT_REAL, &_alglib_env_state, ae_true);
    _nleqreport_init(&_rep, &_alglib_env_state, ae_true);
    nleqresults(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_nleqreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqresultsbuf(const char **errormsg, x_nleqstate** state, x_vector* x, x_nleqreport* rep)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    nleqreport _rep;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    nleqreport_init_from_x(&_rep, rep, &_alglib_env_state, ae_true);
    nleqresultsbuf(&(*state)->obj, &_x, &_rep, &_alglib_env_state);
    ae_x_set_vector(x, &_x, &_alglib_env_state);
    x_set_nleqreport(rep, &_rep, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_nleqrestartfrom(const char **errormsg, x_nleqstate** state, x_vector* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    nleqrestartfrom(&(*state)->obj, &_x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_airy(const char **errormsg, double* x, double* ai, double* aip, double* bi, double* bip)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    airy(*x, ai, aip, bi, bip, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselj0(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselj0(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselj1(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselj1(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besseljn(const char **errormsg, double* result, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besseljn(*n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_bessely0(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = bessely0(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_bessely1(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = bessely1(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselyn(const char **errormsg, double* result, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselyn(*n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besseli0(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besseli0(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besseli1(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besseli1(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselk0(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselk0(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselk1(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselk1(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_besselkn(const char **errormsg, double* result, ae_int_t* nn, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = besselkn(*nn, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_beta(const char **errormsg, double* result, double* a, double* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = beta(*a, *b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_incompletebeta(const char **errormsg, double* result, double* a, double* b, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = incompletebeta(*a, *b, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invincompletebeta(const char **errormsg, double* result, double* a, double* b, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invincompletebeta(*a, *b, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_binomialdistribution(const char **errormsg, double* result, ae_int_t* k, ae_int_t* n, double* p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = binomialdistribution(*k, *n, *p, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_binomialcdistribution(const char **errormsg, double* result, ae_int_t* k, ae_int_t* n, double* p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = binomialcdistribution(*k, *n, *p, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invbinomialdistribution(const char **errormsg, double* result, ae_int_t* k, ae_int_t* n, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invbinomialdistribution(*k, *n, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_chebyshevcalculate(const char **errormsg, double* result, ae_int_t* r, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = chebyshevcalculate(*r, *n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_chebyshevsum(const char **errormsg, double* result, x_vector* c, ae_int_t* r, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *result = chebyshevsum(&_c, *r, *n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_chebyshevcoefficients(const char **errormsg, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    chebyshevcoefficients(*n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fromchebyshev(const char **errormsg, x_vector* a, ae_int_t* n, x_vector* b)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _a;
    ae_vector _b;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_a, a, &_alglib_env_state, ae_true);
    ae_vector_init(&_b, 0, DT_REAL, &_alglib_env_state, ae_true);
    fromchebyshev(&_a, *n, &_b, &_alglib_env_state);
    ae_x_set_vector(b, &_b, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_chisquaredistribution(const char **errormsg, double* result, double* v, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = chisquaredistribution(*v, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_chisquarecdistribution(const char **errormsg, double* result, double* v, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = chisquarecdistribution(*v, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invchisquaredistribution(const char **errormsg, double* result, double* v, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invchisquaredistribution(*v, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_dawsonintegral(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = dawsonintegral(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_ellipticintegralk(const char **errormsg, double* result, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = ellipticintegralk(*m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_ellipticintegralkhighprecision(const char **errormsg, double* result, double* m1)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = ellipticintegralkhighprecision(*m1, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_incompleteellipticintegralk(const char **errormsg, double* result, double* phi, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = incompleteellipticintegralk(*phi, *m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_ellipticintegrale(const char **errormsg, double* result, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = ellipticintegrale(*m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_incompleteellipticintegrale(const char **errormsg, double* result, double* phi, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = incompleteellipticintegrale(*phi, *m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_exponentialintegralei(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = exponentialintegralei(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_exponentialintegralen(const char **errormsg, double* result, double* x, ae_int_t* n)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = exponentialintegralen(*x, *n, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fdistribution(const char **errormsg, double* result, ae_int_t* a, ae_int_t* b, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = fdistribution(*a, *b, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fcdistribution(const char **errormsg, double* result, ae_int_t* a, ae_int_t* b, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = fcdistribution(*a, *b, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invfdistribution(const char **errormsg, double* result, ae_int_t* a, ae_int_t* b, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invfdistribution(*a, *b, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_fresnelintegral(const char **errormsg, double* x, double* c, double* s)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    fresnelintegral(*x, c, s, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hermitecalculate(const char **errormsg, double* result, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = hermitecalculate(*n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hermitesum(const char **errormsg, double* result, x_vector* c, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *result = hermitesum(&_c, *n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hermitecoefficients(const char **errormsg, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    hermitecoefficients(*n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_jacobianellipticfunctions(const char **errormsg, double* u, double* m, double* sn, double* cn, double* dn, double* ph)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    jacobianellipticfunctions(*u, *m, sn, cn, dn, ph, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_laguerrecalculate(const char **errormsg, double* result, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = laguerrecalculate(*n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_laguerresum(const char **errormsg, double* result, x_vector* c, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *result = laguerresum(&_c, *n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_laguerrecoefficients(const char **errormsg, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    laguerrecoefficients(*n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_legendrecalculate(const char **errormsg, double* result, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = legendrecalculate(*n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_legendresum(const char **errormsg, double* result, x_vector* c, ae_int_t* n, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_c, c, &_alglib_env_state, ae_true);
    *result = legendresum(&_c, *n, *x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_legendrecoefficients(const char **errormsg, ae_int_t* n, x_vector* c)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _c;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init(&_c, 0, DT_REAL, &_alglib_env_state, ae_true);
    legendrecoefficients(*n, &_c, &_alglib_env_state);
    ae_x_set_vector(c, &_c, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_poissondistribution(const char **errormsg, double* result, ae_int_t* k, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = poissondistribution(*k, *m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_poissoncdistribution(const char **errormsg, double* result, ae_int_t* k, double* m)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = poissoncdistribution(*k, *m, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invpoissondistribution(const char **errormsg, double* result, ae_int_t* k, double* y)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invpoissondistribution(*k, *y, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_psi(const char **errormsg, double* result, double* x)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = psi(*x, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_studenttdistribution(const char **errormsg, double* result, ae_int_t* k, double* t)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = studenttdistribution(*k, *t, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_invstudenttdistribution(const char **errormsg, double* result, ae_int_t* k, double* p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    *result = invstudenttdistribution(*k, *p, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_sinecosineintegrals(const char **errormsg, double* x, double* si, double* ci)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    sinecosineintegrals(*x, si, ci, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_hyperbolicsinecosineintegrals(const char **errormsg, double* x, double* shi, double* chi)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    hyperbolicsinecosineintegrals(*x, shi, chi, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_pearsoncorrelationsignificance(const char **errormsg, double* r, ae_int_t* n, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    pearsoncorrelationsignificance(*r, *n, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_spearmanrankcorrelationsignificance(const char **errormsg, double* r, ae_int_t* n, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    spearmanrankcorrelationsignificance(*r, *n, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_jarqueberatest(const char **errormsg, x_vector* x, ae_int_t* n, double* p)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    jarqueberatest(&_x, *n, p, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_mannwhitneyutest(const char **errormsg, x_vector* x, ae_int_t* n, x_vector* y, ae_int_t* m, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    mannwhitneyutest(&_x, *n, &_y, *m, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_onesamplesigntest(const char **errormsg, x_vector* x, ae_int_t* n, double* median, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    onesamplesigntest(&_x, *n, *median, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_studentttest1(const char **errormsg, x_vector* x, ae_int_t* n, double* mean, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    studentttest1(&_x, *n, *mean, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_studentttest2(const char **errormsg, x_vector* x, ae_int_t* n, x_vector* y, ae_int_t* m, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    studentttest2(&_x, *n, &_y, *m, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_unequalvariancettest(const char **errormsg, x_vector* x, ae_int_t* n, x_vector* y, ae_int_t* m, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    unequalvariancettest(&_x, *n, &_y, *m, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_ftest(const char **errormsg, x_vector* x, ae_int_t* n, x_vector* y, ae_int_t* m, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_vector _y;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    ae_vector_init_from_x(&_y, y, &_alglib_env_state, ae_true);
    ftest(&_x, *n, &_y, *m, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_onesamplevariancetest(const char **errormsg, x_vector* x, ae_int_t* n, double* variance, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    onesamplevariancetest(&_x, *n, *variance, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
DLLEXPORT int alglib_wilcoxonsignedranktest(const char **errormsg, x_vector* x, ae_int_t* n, double* e, double* bothtails, double* lefttail, double* righttail)
{
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_vector _x;
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        if( _alglib_env_state.last_error==ERR_OUT_OF_MEMORY )    { *errormsg = "ALGLIB: malloc error"; return X_MALLOC_ERROR; }
        if( _alglib_env_state.last_error==ERR_XARRAY_TOO_LARGE ) { *errormsg = "ALGLIB: array too large"; return X_ARRAY_TOO_LARGE; }
        if( _alglib_env_state.last_error==ERR_ASSERTION_FAILED ) { *errormsg = _alglib_env_state.error_msg; return X_ASSERTION_FAILED; }
        return -1;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    ae_vector_init_from_x(&_x, x, &_alglib_env_state, ae_true);
    wilcoxonsignedranktest(&_x, *n, *e, bothtails, lefttail, righttail, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    return X_OK;
}
