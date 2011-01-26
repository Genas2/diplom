#########################################################################
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (www.fsf.org); either version 2 of the 
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# A copy of the GNU General Public License is available at
# http://www.fsf.org/licensing/licenses
#########################################################################
import ctypes

import ctypes
import sys
import os

if ctypes.sizeof(ctypes.c_void_p)==4:
    c_ptrint_t = ctypes.c_int32
else:
    c_ptrint_t = ctypes.c_int64

DT_BOOL = 1
DT_INT = 2
DT_REAL = 3
DT_COMPLEX = 4

SIZE_BOOL = 1
SIZE_INT = 4
SIZE_REAL = 8
SIZE_COMPLEX = 16

OWN_CALLER = 1
OWN_AE = 2

X_ASSERTION_FAILED = 5

X_SET     = 1 # data are copied into x-vector/matrix; previous contents of x-structure is freed
X_CREATE  = 2 # x-vector/matrix is created, its previous contents is ignored
X_REWRITE = 3 # data are copied into x-structure; size of Python structure must be equal to the x-structure size

class x_complex(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double)]

class x_multiptr(ctypes.Union):
    _fields_ = [("p_ptr",       ctypes.c_void_p),
                ("p_char",      ctypes.c_char_p),
                ("p_bool",      ctypes.POINTER(ctypes.c_int8)),
                ("p_int",       ctypes.POINTER(c_ptrint_t)),
                ("p_real",      ctypes.POINTER(ctypes.c_double)),
                ("p_complex",   ctypes.POINTER(x_complex))]

class x_string(ctypes.Structure):
    _fields_ = [("owner",       ctypes.c_uint64),
                ("last_action", ctypes.c_uint64),
                ("ptr",         x_multiptr)]

class x_vector(ctypes.Structure):
    _fields_ = [("cnt",         ctypes.c_uint64),
                ("datatype",    ctypes.c_uint64),
                ("owner",       ctypes.c_uint64),
                ("last_action", ctypes.c_uint64),
                ("ptr",         x_multiptr)]

class x_matrix(ctypes.Structure):
    _fields_ = [("rows",        ctypes.c_uint64),
                ("cols",        ctypes.c_uint64),
                ("stride",      ctypes.c_uint64),
                ("datatype",    ctypes.c_uint64),
                ("owner",       ctypes.c_uint64),
                ("last_action", ctypes.c_uint64),
                ("ptr",         x_multiptr)]

#
# Load ALGLIB shared library
#
_so_candidates = []
if sys.platform=="win32" or sys.platform=="cygwin":
    _shared_lib_name = "alglib.dll"
    _so_candidates.append(os.path.join(os.path.dirname(__file__),_shared_lib_name))
    _so_candidates.append(os.path.join(sys.prefix,_shared_lib_name))
else:
    _shared_lib_name = "alglib.so"
    _so_candidates.append(os.path.join(os.path.dirname(__file__),_shared_lib_name))
    _so_candidates.append(os.path.join(sys.prefix,_shared_lib_name))
    _so_candidates.append(os.path.join("usr","local",_shared_lib_name))
for _candidate in _so_candidates:
    if os.path.isfile(_candidate):
        _lib_alglib = ctypes.cdll.LoadLibrary(_candidate)
        break
if _lib_alglib is None:
    raise RuntimeError("ALGLIB: unable to load shared library")

_lib_alglib.x_malloc.argtypes = [ctypes.c_void_p, ctypes.c_uint64]
_lib_alglib.x_malloc.restype = ctypes.c_int32
_lib_alglib.x_free.argtypes = [ctypes.c_void_p]
_lib_alglib.x_free.restype = ctypes.c_int32
_lib_alglib.x_alloc_counter.argtypes = []
_lib_alglib.x_alloc_counter.restype = ctypes.c_int64
_lib_alglib.x_is_symmetric_e_.argtypes = [ctypes.c_void_p]
_lib_alglib.x_is_symmetric_e_.restype = ctypes.c_uint8
_lib_alglib.x_is_hermitian_e_.argtypes = [ctypes.c_void_p]
_lib_alglib.x_is_hermitian_e_.restype = ctypes.c_uint8
_lib_alglib.x_force_symmetric_e_.argtypes = [ctypes.c_void_p]
_lib_alglib.x_force_symmetric_e_.restype = ctypes.c_uint8
_lib_alglib.x_force_hermitian_e_.argtypes = [ctypes.c_void_p]
_lib_alglib.x_force_hermitian_e_.restype = ctypes.c_uint8

def x_malloc(cnt):
    __cnt = ctypes.c_uint64(cnt)
    if __cnt.value!=cnt:
        raise RuntimeError("malloc() argument is too large!")
    __result = ctypes.c_void_p(0)
    if _lib_alglib.x_malloc(ctypes.byref(__result), __cnt)!=0:
        raise RuntimeError("Error while calling x_malloc()")
    return __result.value

def x_free(ptr):
    __ptr = ctypes.c_void_p(ptr)
    _lib_alglib.x_free(__ptr)
    return

def x_alloc_counter():
    return _lib_alglib.x_alloc_counter()

def x_is_symmetric(x):
    return _lib_alglib.x_is_symmetric_e_(ctypes.byref(x))

def x_force_symmetric(x):
    return _lib_alglib.x_force_symmetric_e_(ctypes.byref(x))

def x_is_hermitian(x):
    return _lib_alglib.x_is_hermitian_e_(ctypes.byref(x))

def x_force_hermitian(x):
    return _lib_alglib.x_force_hermitian_e_(ctypes.byref(x))

def x_vector_clear(x):
    if x.owner==OWN_AE:
        x_free(x.ptr.p_ptr)
    x.cnt = 0
    x.ptr.p_ptr = 0
    return

def x_matrix_clear(x):
    if x.owner==OWN_AE:
        x_free(x.ptr.p_ptr)
    x.cols = 0
    x.rows = 0
    x.ptr.p_ptr = 0
    return

#
# safe vector length:
# * returns list length.
# * throws ValueError if 'v' is not list (it uses
#   'msg' parameter to generate error message
#
def safe_len(msg,v):
    if type(v)!=list:
        raise ValueError(msg)
    return len(v)

#
# safe matrix size
# * returns number of columns
# * throws ValueError if 'v' is not rectangular matrix 
#   (list of lists of same size)
#   it uses 'msg' parameter to generate error message
#
def safe_cols(msg,v):
    if type(v)!=list:
        raise ValueError(msg)
    if len(v)==0:
        return 0
    if type(v[0])!=list:
        raise ValueError(msg)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            raise ValueError(msg)
        if len(x)!=cols:
            raise ValueError(msg)
    return cols

#
# safe matrix size
# * returns number of rows
# * throws ValueError if 'v' is not rectangular matrix 
#   (list of lists of same size)
#   it uses 'msg' parameter to generate error message
#
def safe_rows(msg,v):
    if type(v)!=list:
        raise ValueError(msg)
    if len(v)==0:
        return 0
    if type(v[0])!=list:
        raise ValueError(msg)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            raise ValueError(msg)
        if len(x)!=cols:
            raise ValueError(msg)
    return len(v)

def create_real_vector(cnt):
    if cnt<=0:
        return []
    return [0]*cnt

def create_real_matrix(rows, cols):
    if rows<=0 or cols<=0:
        return [[]]
    matrix = []
    row = 0
    while row<rows:
        matrix += [[0]*cols]
        row += 1
    return matrix

def is_bool(v):
    try:
        tmp = bool(v)
    except:
        return False
    return True

def is_int(v):
    try:
        tmp = int(v)
    except:
        return False
    return True

def is_real(v):
    try:
        tmp = float(v)
    except:
        return False
    return True

def is_complex(v):
    try:
        tmp = complex(v)
    except:
        return False
    return True

def is_bool_vector(v):
    if type(v)!=list:
        return False
    for x in v:
        try:
            tmp = bool(x)
        except:
            return False
    return True

def is_bool_matrix(v):
    if type(v)!=list:
        return False
    if len(v)==0:
        return True
    if type(v[0])!=list:
        return False
    rows = len(v)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            return False
        if len(x)!=cols:
            return False
        for y in x:
            try:
                tmp = bool(y)
            except:
                return False
    return True

def is_int_vector(v):
    if type(v)!=list:
        return False
    for x in v:
        try:
            tmp = int(x)
        except:
            return False
    return True

def is_int_matrix(v):
    if type(v)!=list:
        return False
    if len(v)==0:
        return True
    if type(v[0])!=list:
        return False
    rows = len(v)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            return False
        if len(x)!=cols:
            return False
        for y in x:
            try:
                tmp = int(y)
            except:
                return False
    return True

def is_real_vector(v):
    if type(v)!=list:
        return False
    for x in v:
        try:
            tmp = float(x)
        except:
            return False
    return True

def is_real_matrix(v):
    if type(v)!=list:
        return False
    if len(v)==0:
        return True
    if type(v[0])!=list:
        return False
    rows = len(v)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            return False
        if len(x)!=cols:
            return False
        for y in x:
            try:
                tmp = float(y)
            except:
                return False
    return True

def is_complex_vector(v):
    if type(v)!=list:
        return False
    for x in v:
        try:
            tmp = complex(x)
        except:
            return False
    return True

def is_complex_matrix(v):
    if type(v)!=list:
        return False
    if len(v)==0:
        return True
    if type(v[0])!=list:
        return False
    rows = len(v)
    cols = len(v[0])
    for x in v:
        if type(x)!=list:
            return False
        if len(x)!=cols:
            return False
        for y in x:
            try:
                tmp = complex(y)
            except:
                return False
    return True

#
# conversion from list to x-vector:
#
# x     x-vector. 
# v     list
# dt    datatype
# mode  one of the modes:
#       * X_CREATE -  x is assumed to be uninitialized
#       * X_SET -     x is assumed to be initialized; old contents is freed
#       * X_REWRITE - x is assumed to have same size as v, exception is thrown otherwise
#                     it is rewritten without reallocation of memory;
#
def x_from_list(x, v, dt, mode):
    #
    # check types
    #
    if dt==DT_BOOL:
        if not is_bool_vector(v):
            raise ValueError("can't cast to bool_vector")
        elemsize = SIZE_BOOL
    if dt==DT_INT:
        if not is_int_vector(v):
            raise ValueError("can't cast to int_vector")
        elemsize = SIZE_INT
    if dt==DT_REAL:
        if not is_real_vector(v):
            raise ValueError("can't cast to real_vector")
        elemsize = SIZE_REAL
    if dt==DT_COMPLEX:
        if not is_complex_vector(v):
            raise ValueError("can't cast to complex_vector")
        elemsize = SIZE_COMPLEX
    
    #
    # allocation
    #
    if mode==X_CREATE:
        x.datatype = dt
        x.cnt = len(v)
        x.owner = OWN_AE
        x.ptr.p_ptr = x_malloc(elemsize*x.cnt)
        x.last_action = 1
    if mode==X_SET:
        if x.owner==OWN_AE:
            x_free(x.ptr.p_ptr)
        x.datatype = dt
        x.cnt = len(v)
        x.owner = OWN_AE
        x.ptr.p_ptr = x_malloc(elemsize*x.cnt)
        x.last_action = 3
    if mode==X_REWRITE:
        if x.datatype!=dt:
            raise RuntimeError("Trying to rewrite vector - types don't match")
        if len(v)!=x.cnt:
            raise RuntimeError("Trying to rewrite vector - sizes don't match")
        x.last_action = 2
    
    #
    # copy
    #
    if dt==DT_BOOL:
        cnt = x.cnt
        i = 0
        while i<cnt:
            x.ptr.p_bool[i] =   bool(v[i])
            i += 1
    if dt==DT_INT:
        cnt = x.cnt
        i = 0
        while i<cnt:
            x.ptr.p_int[i] =    int(v[i])
            i += 1
    if dt==DT_REAL:
        cnt = x.cnt
        i = 0
        while i<cnt:
            x.ptr.p_real[i] = float(v[i])
            i += 1
    if dt==DT_COMPLEX:
        cnt = x.cnt
        i = 0
        while i<cnt:
            tmp = complex(v[i])
            x.ptr.p_complex[i].x = tmp.real
            x.ptr.p_complex[i].y = tmp.imag
            i += 1
    
    return


#
# conversion from list of lists to x-matrix:
#
# x     x-matrix.
# v     list
# dt    datatype
# mode  one of the modes:
#       * X_CREATE -  x is assumed to be uninitialized
#       * X_SET -     x is assumed to be initialized; old contents is freed
#       * X_REWRITE - x is assumed to have same size as v, exception is thrown otherwise
#                     it is rewritten without reallocation of memory;
#
def x_from_listlist(x, v, dt, mode):
    #
    # check types
    #
    if dt==DT_BOOL:
        if not is_bool_matrix(v):
            raise ValueError("can't cast to bool_matrix")
        elemsize = SIZE_BOOL
    if dt==DT_INT:
        if not is_int_matrix(v):
            raise ValueError("can't cast to int_matrix")
        elemsize = SIZE_INT
    if dt==DT_REAL:
        if not is_real_matrix(v):
            raise ValueError("can't cast to real_matrix")
        elemsize = SIZE_REAL
    if dt==DT_COMPLEX:
        if not is_complex_matrix(v):
            raise ValueError("can't cast to complex_matrix")
        elemsize = SIZE_COMPLEX
    
    #
    # determine size
    #
    rows = len(v)
    if rows>0:
        cols = len(v[0])
    else:
        cols = 0
    if cols==0:
        rows = 0
    
    #
    # allocation
    #
    if mode==X_CREATE:
        x.datatype = dt
        x.cols = cols
        x.rows = rows
        x.stride = cols
        x.owner = OWN_AE
        x.ptr.p_ptr = x_malloc(elemsize*x.stride*x.rows)
        x.last_action = 1
    if mode==X_SET:
        if x.owner==OWN_AE:
            x_free(x.ptr.p_ptr)
        x.cols = cols
        x.rows = rows
        x.stride =cols
        x.owner = OWN_AE
        x.ptr.p_ptr = x_malloc(elemsize*x.stride*x.rows)
        x.last_action = 3
    if mode==X_REWRITE:
        if x.datatype!=dt:
            raise RuntimeError("Trying to rewrite matrix - types don't match")
        if rows!=x.rows or cols!=x.cols:
            raise RuntimeError("Trying to rewrite vector - sizes don't match")
        x.last_action = 2
    
    #
    # copy
    #
    offs = 0
    endoffs = x.stride-x.cols
    for p in v:
        for q in p:
            if dt==DT_BOOL:
                x.ptr.p_bool[offs] = bool(q)
            if dt==DT_INT:
                x.ptr.p_int[offs]  = int(q)
            if dt==DT_REAL:
                x.ptr.p_real[offs] = float(q)
            if dt==DT_COMPLEX:
                tmp = complex(q)
                x.ptr.p_complex[offs].x = tmp.real
                x.ptr.p_complex[offs].y = tmp.imag
            offs = offs+1
        offs = offs+endoffs
    return


#
# conversion from x-vector to Python vector
#
# Function takes only one parameter - x, x-vector,
# which is NOT freed after use.
#
def list_from_x(x):
    if x.cnt==0:
        return []
    r = [0]*x.cnt
    i = 0
    if x.datatype==DT_BOOL:
        while i<x.cnt:
            r[i] = bool(x.ptr.p_bool[i])
            i += 1
    if x.datatype==DT_INT:
        while i<x.cnt:
            r[i] = x.ptr.p_int[i]
            i += 1
    if x.datatype==DT_REAL:
        while i<x.cnt:
            r[i] = x.ptr.p_real[i]
            i += 1
    if x.datatype==DT_COMPLEX:
        while i<x.cnt:
            r[i] = complex(x.ptr.p_complex[i].x, x.ptr.p_complex[i].y)
            i += 1
    return r


#
# conversion from x-matrix to Python matrix
#
# Function takes only one parameter - x, x-matrix,
# which is NOT freed after use.
#
def listlist_from_x(x):
    if x.cols==0 or x.rows==0:
        return [[]]
    r = create_real_matrix(x.rows, x.cols)
    offs = 0
    endoffs = x.stride-x.cols
    m = x.rows
    n = x.cols
    dt = x.datatype
    i = 0
    while i<m:
        j = 0
        while j<n:
            if dt==DT_BOOL:
                r[i][j] = bool(x.ptr.p_bool[offs])
            if dt==DT_INT:
                r[i][j] = x.ptr.p_int[offs]
            if dt==DT_REAL:
                r[i][j] = x.ptr.p_real[offs]
            if dt==DT_COMPLEX:
                r[i][j] = complex(x.ptr.p_complex[offs].x, x.ptr.p_complex[offs].y)
            offs += 1
            j += 1
        offs += endoffs
        i += 1
    return r



#
# this function copies x-vector to previously allocated list 
# which should be large enough to store x-vector.
#
# invalid access to list is generated if list is too small.
# x-vector is not freed after use.
#
def copy_x_to_list(x,r):
    if x.cnt==0:
        return
    i = 0
    if x.datatype==DT_BOOL:
        while i<x.cnt:
            r[i] = bool(x.ptr.p_bool[i])
            i += 1
    if x.datatype==DT_INT:
        while i<x.cnt:
            r[i] = x.ptr.p_int[i]
            i += 1
    if x.datatype==DT_REAL:
        while i<x.cnt:
            r[i] = x.ptr.p_real[i]
            i += 1
    if x.datatype==DT_COMPLEX:
        while i<x.cnt:
            r[i] = complex(x.ptr.p_complex[i].x, x.ptr.p_complex[i].y)
            i += 1
    return

_lib_alglib.x_obj_free_hqrndstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_hqrndstate.restype = None


class hqrndstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_hqrndstate(self.ptr)
_lib_alglib.alglib_hqrndrandomize.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndrandomize.restype = ctypes.c_int32
def hqrndrandomize():
    pass
    __state = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndrandomize(ctypes.byref(_error_msg), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndrandomize'")
        __r__state = hqrndstate(__state)
        return __r__state
    finally:
        pass


_lib_alglib.alglib_hqrndseed.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndseed.restype = ctypes.c_int32
def hqrndseed(s1, s2):
    pass
    __s1 = c_ptrint_t(s1)
    if __s1.value!=s1:
        raise ValueError("Error while converting 's1' parameter to 'c_ptrint_t'")
    __s2 = c_ptrint_t(s2)
    if __s2.value!=s2:
        raise ValueError("Error while converting 's2' parameter to 'c_ptrint_t'")
    __state = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndseed(ctypes.byref(_error_msg), ctypes.byref(__s1), ctypes.byref(__s2), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndseed'")
        __r__state = hqrndstate(__state)
        return __r__state
    finally:
        pass


_lib_alglib.alglib_hqrnduniformr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrnduniformr.restype = ctypes.c_int32
def hqrnduniformr(state):
    pass
    __result = ctypes.c_double(0)
    __state = state.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrnduniformr(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrnduniformr'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_hqrnduniformi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrnduniformi.restype = ctypes.c_int32
def hqrnduniformi(state, n):
    pass
    __result = c_ptrint_t(0)
    __state = state.ptr
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrnduniformi(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__state), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrnduniformi'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_hqrndnormal.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndnormal.restype = ctypes.c_int32
def hqrndnormal(state):
    pass
    __result = ctypes.c_double(0)
    __state = state.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndnormal(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndnormal'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_hqrndunit2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndunit2.restype = ctypes.c_int32
def hqrndunit2(state):
    pass
    __state = state.ptr
    __x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndunit2(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndunit2'")
        __r__x = __x.value
        __r__y = __y.value
        return (__r__x, __r__y)
    finally:
        pass


_lib_alglib.alglib_hqrndnormal2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndnormal2.restype = ctypes.c_int32
def hqrndnormal2(state):
    pass
    __state = state.ptr
    __x1 = ctypes.c_double(0)
    __x2 = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndnormal2(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x1), ctypes.byref(__x2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndnormal2'")
        __r__x1 = __x1.value
        __r__x2 = __x2.value
        return (__r__x1, __r__x2)
    finally:
        pass


_lib_alglib.alglib_hqrndexponential.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hqrndexponential.restype = ctypes.c_int32
def hqrndexponential(state, lambdav):
    pass
    __result = ctypes.c_double(0)
    __state = state.ptr
    __lambdav = ctypes.c_double(lambdav)
    if __lambdav.value!=lambdav:
        raise ValueError("Error while converting 'lambdav' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hqrndexponential(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__state), ctypes.byref(__lambdav))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hqrndexponential'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.x_obj_free_kdtree.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_kdtree.restype = None


class kdtree(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_kdtree(self.ptr)
_lib_alglib.alglib_kdtreebuild.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreebuild.restype = ctypes.c_int32
def kdtreebuild(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        xy,n,nx,ny,normtype = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        xy,nx,ny,normtype = functionargs
        n = safe_rows("'kdtreebuild': incorrect parameters",xy)
    else:
        raise RuntimeError("Error while calling 'kdtreebuild': function must have 4 or 5 parameters")
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nx = c_ptrint_t(nx)
    if __nx.value!=nx:
        raise ValueError("Error while converting 'nx' parameter to 'c_ptrint_t'")
    __ny = c_ptrint_t(ny)
    if __ny.value!=ny:
        raise ValueError("Error while converting 'ny' parameter to 'c_ptrint_t'")
    __normtype = c_ptrint_t(normtype)
    if __normtype.value!=normtype:
        raise ValueError("Error while converting 'normtype' parameter to 'c_ptrint_t'")
    __kdt = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreebuild(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__nx), ctypes.byref(__ny), ctypes.byref(__normtype), ctypes.byref(__kdt))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreebuild'")
        __r__kdt = kdtree(__kdt)
        return __r__kdt
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_kdtreebuildtagged.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreebuildtagged.restype = ctypes.c_int32
def kdtreebuildtagged(*functionargs):
    if len(functionargs)==6:
        __friendly_form = False
        xy,tags,n,nx,ny,normtype = functionargs
    elif len(functionargs)==5:
        __friendly_form = True
        xy,tags,nx,ny,normtype = functionargs
        if safe_rows("'kdtreebuildtagged': incorrect parameters",xy)!=safe_len("'kdtreebuildtagged': incorrect parameters",tags):
            raise RuntimeError("Error while calling 'kdtreebuildtagged': looks like one of arguments has wrong size")
        n = safe_rows("'kdtreebuildtagged': incorrect parameters",xy)
    else:
        raise RuntimeError("Error while calling 'kdtreebuildtagged': function must have 5 or 6 parameters")
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(tags):
        raise ValueError("'tags' parameter can't be cast to int_vector")
    __tags = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nx = c_ptrint_t(nx)
    if __nx.value!=nx:
        raise ValueError("Error while converting 'nx' parameter to 'c_ptrint_t'")
    __ny = c_ptrint_t(ny)
    if __ny.value!=ny:
        raise ValueError("Error while converting 'ny' parameter to 'c_ptrint_t'")
    __normtype = c_ptrint_t(normtype)
    if __normtype.value!=normtype:
        raise ValueError("Error while converting 'normtype' parameter to 'c_ptrint_t'")
    __kdt = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__tags, tags, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreebuildtagged(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__tags), ctypes.byref(__n), ctypes.byref(__nx), ctypes.byref(__ny), ctypes.byref(__normtype), ctypes.byref(__kdt))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreebuildtagged'")
        __r__kdt = kdtree(__kdt)
        return __r__kdt
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__tags)


_lib_alglib.alglib_kdtreequeryknn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryknn.restype = ctypes.c_int32
def kdtreequeryknn(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        kdt,x,k,selfmatch = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        kdt,x,k = functionargs
        selfmatch = True
    else:
        raise RuntimeError("Error while calling 'kdtreequeryknn': function must have 3 or 4 parameters")
    __result = c_ptrint_t(0)
    __kdt = kdt.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __selfmatch = ctypes.c_uint8(selfmatch)
    if __selfmatch.value!=0:
        __selfmatch = ctypes.c_uint8(1)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryknn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__kdt), ctypes.byref(__x), ctypes.byref(__k), ctypes.byref(__selfmatch))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryknn'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_kdtreequeryrnn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryrnn.restype = ctypes.c_int32
def kdtreequeryrnn(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        kdt,x,r,selfmatch = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        kdt,x,r = functionargs
        selfmatch = True
    else:
        raise RuntimeError("Error while calling 'kdtreequeryrnn': function must have 3 or 4 parameters")
    __result = c_ptrint_t(0)
    __kdt = kdt.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __r = ctypes.c_double(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'ctypes.c_double'")
    __selfmatch = ctypes.c_uint8(selfmatch)
    if __selfmatch.value!=0:
        __selfmatch = ctypes.c_uint8(1)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryrnn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__kdt), ctypes.byref(__x), ctypes.byref(__r), ctypes.byref(__selfmatch))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryrnn'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_kdtreequeryaknn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryaknn.restype = ctypes.c_int32
def kdtreequeryaknn(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        kdt,x,k,selfmatch,eps = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        kdt,x,k,eps = functionargs
        selfmatch = True
    else:
        raise RuntimeError("Error while calling 'kdtreequeryaknn': function must have 4 or 5 parameters")
    __result = c_ptrint_t(0)
    __kdt = kdt.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __selfmatch = ctypes.c_uint8(selfmatch)
    if __selfmatch.value!=0:
        __selfmatch = ctypes.c_uint8(1)
    __eps = ctypes.c_double(eps)
    if __eps.value!=eps:
        raise ValueError("Error while converting 'eps' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryaknn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__kdt), ctypes.byref(__x), ctypes.byref(__k), ctypes.byref(__selfmatch), ctypes.byref(__eps))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryaknn'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_kdtreequeryresultsx.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsx.restype = ctypes.c_int32
def kdtreequeryresultsx(kdt, x):
    pass
    __kdt = kdt.ptr
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsx(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsx'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__x)


_lib_alglib.alglib_kdtreequeryresultsxy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsxy.restype = ctypes.c_int32
def kdtreequeryresultsxy(kdt, xy):
    pass
    __kdt = kdt.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsxy(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__xy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsxy'")
        __r__xy = listlist_from_x(__xy)
        return __r__xy
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_kdtreequeryresultstags.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultstags.restype = ctypes.c_int32
def kdtreequeryresultstags(kdt, tags):
    pass
    __kdt = kdt.ptr
    if not is_int_vector(tags):
        raise ValueError("'tags' parameter can't be cast to int_vector")
    __tags = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__tags, tags, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultstags(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__tags))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultstags'")
        __r__tags = list_from_x(__tags)
        return __r__tags
    finally:
        x_vector_clear(__tags)


_lib_alglib.alglib_kdtreequeryresultsdistances.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsdistances.restype = ctypes.c_int32
def kdtreequeryresultsdistances(kdt, r):
    pass
    __kdt = kdt.ptr
    if not is_real_vector(r):
        raise ValueError("'r' parameter can't be cast to real_vector")
    __r = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__r, r, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsdistances(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsdistances'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__r)


_lib_alglib.alglib_kdtreequeryresultsxi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsxi.restype = ctypes.c_int32
def kdtreequeryresultsxi(kdt):
    pass
    __kdt = kdt.ptr
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsxi(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsxi'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__x)


_lib_alglib.alglib_kdtreequeryresultsxyi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsxyi.restype = ctypes.c_int32
def kdtreequeryresultsxyi(kdt):
    pass
    __kdt = kdt.ptr
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsxyi(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__xy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsxyi'")
        __r__xy = listlist_from_x(__xy)
        return __r__xy
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_kdtreequeryresultstagsi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultstagsi.restype = ctypes.c_int32
def kdtreequeryresultstagsi(kdt):
    pass
    __kdt = kdt.ptr
    __tags = x_vector(cnt=0,datatype=DT_INT,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultstagsi(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__tags))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultstagsi'")
        __r__tags = list_from_x(__tags)
        return __r__tags
    finally:
        x_vector_clear(__tags)


_lib_alglib.alglib_kdtreequeryresultsdistancesi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kdtreequeryresultsdistancesi.restype = ctypes.c_int32
def kdtreequeryresultsdistancesi(kdt):
    pass
    __kdt = kdt.ptr
    __r = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kdtreequeryresultsdistancesi(ctypes.byref(_error_msg), ctypes.byref(__kdt), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kdtreequeryresultsdistancesi'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__r)


_lib_alglib.alglib_cmatrixtranspose.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixtranspose.restype = ctypes.c_int32
def cmatrixtranspose(m, n, a, ia, ja, b, ib, jb):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixtranspose(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixtranspose'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_rmatrixtranspose.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixtranspose.restype = ctypes.c_int32
def rmatrixtranspose(m, n, a, ia, ja, b, ib, jb):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixtranspose(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixtranspose'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_cmatrixcopy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixcopy.restype = ctypes.c_int32
def cmatrixcopy(m, n, a, ia, ja, b, ib, jb):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixcopy(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixcopy'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_rmatrixcopy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixcopy.restype = ctypes.c_int32
def rmatrixcopy(m, n, a, ia, ja, b, ib, jb):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixcopy(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixcopy'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_cmatrixrank1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrank1.restype = ctypes.c_int32
def cmatrixrank1(m, n, a, ia, ja, u, iu, v, iv):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_complex_vector(u):
        raise ValueError("'u' parameter can't be cast to complex_vector")
    __u = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iu = c_ptrint_t(iu)
    if __iu.value!=iu:
        raise ValueError("Error while converting 'iu' parameter to 'c_ptrint_t'")
    if not is_complex_vector(v):
        raise ValueError("'v' parameter can't be cast to complex_vector")
    __v = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iv = c_ptrint_t(iv)
    if __iv.value!=iv:
        raise ValueError("Error while converting 'iv' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__u, u, DT_COMPLEX, X_CREATE)
        x_from_list(__v, v, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrank1(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__u), ctypes.byref(__iu), ctypes.byref(__v), ctypes.byref(__iv))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrank1'")
        __r__a = listlist_from_x(__a)
        __r__u = list_from_x(__u)
        __r__v = list_from_x(__v)
        return (__r__a, __r__u, __r__v)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__u)
        x_vector_clear(__v)


_lib_alglib.alglib_rmatrixrank1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrank1.restype = ctypes.c_int32
def rmatrixrank1(m, n, a, ia, ja, u, iu, v, iv):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    if not is_real_vector(u):
        raise ValueError("'u' parameter can't be cast to real_vector")
    __u = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iu = c_ptrint_t(iu)
    if __iu.value!=iu:
        raise ValueError("Error while converting 'iu' parameter to 'c_ptrint_t'")
    if not is_real_vector(v):
        raise ValueError("'v' parameter can't be cast to real_vector")
    __v = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iv = c_ptrint_t(iv)
    if __iv.value!=iv:
        raise ValueError("Error while converting 'iv' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__u, u, DT_REAL, X_CREATE)
        x_from_list(__v, v, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrank1(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__u), ctypes.byref(__iu), ctypes.byref(__v), ctypes.byref(__iv))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrank1'")
        __r__a = listlist_from_x(__a)
        __r__u = list_from_x(__u)
        __r__v = list_from_x(__v)
        return (__r__a, __r__u, __r__v)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__u)
        x_vector_clear(__v)


_lib_alglib.alglib_cmatrixmv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixmv.restype = ctypes.c_int32
def cmatrixmv(m, n, a, ia, ja, opa, x, ix, y, iy):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __opa = c_ptrint_t(opa)
    if __opa.value!=opa:
        raise ValueError("Error while converting 'opa' parameter to 'c_ptrint_t'")
    if not is_complex_vector(x):
        raise ValueError("'x' parameter can't be cast to complex_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ix = c_ptrint_t(ix)
    if __ix.value!=ix:
        raise ValueError("Error while converting 'ix' parameter to 'c_ptrint_t'")
    if not is_complex_vector(y):
        raise ValueError("'y' parameter can't be cast to complex_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iy = c_ptrint_t(iy)
    if __iy.value!=iy:
        raise ValueError("Error while converting 'iy' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__x, x, DT_COMPLEX, X_CREATE)
        x_from_list(__y, y, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixmv(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__opa), ctypes.byref(__x), ctypes.byref(__ix), ctypes.byref(__y), ctypes.byref(__iy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixmv'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_rmatrixmv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixmv.restype = ctypes.c_int32
def rmatrixmv(m, n, a, ia, ja, opa, x, ix, y, iy):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __opa = c_ptrint_t(opa)
    if __opa.value!=opa:
        raise ValueError("Error while converting 'opa' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ix = c_ptrint_t(ix)
    if __ix.value!=ix:
        raise ValueError("Error while converting 'ix' parameter to 'c_ptrint_t'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __iy = c_ptrint_t(iy)
    if __iy.value!=iy:
        raise ValueError("Error while converting 'iy' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixmv(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__opa), ctypes.byref(__x), ctypes.byref(__ix), ctypes.byref(__y), ctypes.byref(__iy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixmv'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_cmatrixrighttrsm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrighttrsm.restype = ctypes.c_int32
def cmatrixrighttrsm(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __j1 = c_ptrint_t(j1)
    if __j1.value!=j1:
        raise ValueError("Error while converting 'j1' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __optype = c_ptrint_t(optype)
    if __optype.value!=optype:
        raise ValueError("Error while converting 'optype' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(x):
        raise ValueError("'x' parameter can't be cast to complex_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __j2 = c_ptrint_t(j2)
    if __j2.value!=j2:
        raise ValueError("Error while converting 'j2' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__x, x, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrighttrsm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__i1), ctypes.byref(__j1), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__optype), ctypes.byref(__x), ctypes.byref(__i2), ctypes.byref(__j2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrighttrsm'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixlefttrsm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlefttrsm.restype = ctypes.c_int32
def cmatrixlefttrsm(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __j1 = c_ptrint_t(j1)
    if __j1.value!=j1:
        raise ValueError("Error while converting 'j1' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __optype = c_ptrint_t(optype)
    if __optype.value!=optype:
        raise ValueError("Error while converting 'optype' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(x):
        raise ValueError("'x' parameter can't be cast to complex_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __j2 = c_ptrint_t(j2)
    if __j2.value!=j2:
        raise ValueError("Error while converting 'j2' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__x, x, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlefttrsm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__i1), ctypes.byref(__j1), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__optype), ctypes.byref(__x), ctypes.byref(__i2), ctypes.byref(__j2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlefttrsm'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__x)


_lib_alglib.alglib_rmatrixrighttrsm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrighttrsm.restype = ctypes.c_int32
def rmatrixrighttrsm(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __j1 = c_ptrint_t(j1)
    if __j1.value!=j1:
        raise ValueError("Error while converting 'j1' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __optype = c_ptrint_t(optype)
    if __optype.value!=optype:
        raise ValueError("Error while converting 'optype' parameter to 'c_ptrint_t'")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __j2 = c_ptrint_t(j2)
    if __j2.value!=j2:
        raise ValueError("Error while converting 'j2' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrighttrsm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__i1), ctypes.byref(__j1), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__optype), ctypes.byref(__x), ctypes.byref(__i2), ctypes.byref(__j2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrighttrsm'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__x)


_lib_alglib.alglib_rmatrixlefttrsm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlefttrsm.restype = ctypes.c_int32
def rmatrixlefttrsm(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __j1 = c_ptrint_t(j1)
    if __j1.value!=j1:
        raise ValueError("Error while converting 'j1' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __optype = c_ptrint_t(optype)
    if __optype.value!=optype:
        raise ValueError("Error while converting 'optype' parameter to 'c_ptrint_t'")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __j2 = c_ptrint_t(j2)
    if __j2.value!=j2:
        raise ValueError("Error while converting 'j2' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlefttrsm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__i1), ctypes.byref(__j1), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__optype), ctypes.byref(__x), ctypes.byref(__i2), ctypes.byref(__j2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlefttrsm'")
        __r__x = listlist_from_x(__x)
        return __r__x
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixsyrk.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixsyrk.restype = ctypes.c_int32
def cmatrixsyrk(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __optypea = c_ptrint_t(optypea)
    if __optypea.value!=optypea:
        raise ValueError("Error while converting 'optypea' parameter to 'c_ptrint_t'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    if not is_complex_matrix(c):
        raise ValueError("'c' parameter can't be cast to complex_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ic = c_ptrint_t(ic)
    if __ic.value!=ic:
        raise ValueError("Error while converting 'ic' parameter to 'c_ptrint_t'")
    __jc = c_ptrint_t(jc)
    if __jc.value!=jc:
        raise ValueError("Error while converting 'jc' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__c, c, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixsyrk(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__k), ctypes.byref(__alpha), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__optypea), ctypes.byref(__beta), ctypes.byref(__c), ctypes.byref(__ic), ctypes.byref(__jc), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixsyrk'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__c)


_lib_alglib.alglib_rmatrixsyrk.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixsyrk.restype = ctypes.c_int32
def rmatrixsyrk(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __optypea = c_ptrint_t(optypea)
    if __optypea.value!=optypea:
        raise ValueError("Error while converting 'optypea' parameter to 'c_ptrint_t'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    if not is_real_matrix(c):
        raise ValueError("'c' parameter can't be cast to real_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ic = c_ptrint_t(ic)
    if __ic.value!=ic:
        raise ValueError("Error while converting 'ic' parameter to 'c_ptrint_t'")
    __jc = c_ptrint_t(jc)
    if __jc.value!=jc:
        raise ValueError("Error while converting 'jc' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixsyrk(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__k), ctypes.byref(__alpha), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__optypea), ctypes.byref(__beta), ctypes.byref(__c), ctypes.byref(__ic), ctypes.byref(__jc), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixsyrk'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__c)


_lib_alglib.alglib_cmatrixgemm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixgemm.restype = ctypes.c_int32
def cmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __tmp__val = complex(alpha)
    __alpha = x_complex(x=__tmp__val.real, y=__tmp__val.imag)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __optypea = c_ptrint_t(optypea)
    if __optypea.value!=optypea:
        raise ValueError("Error while converting 'optypea' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    __optypeb = c_ptrint_t(optypeb)
    if __optypeb.value!=optypeb:
        raise ValueError("Error while converting 'optypeb' parameter to 'c_ptrint_t'")
    __tmp__val = complex(beta)
    __beta = x_complex(x=__tmp__val.real, y=__tmp__val.imag)
    if not is_complex_matrix(c):
        raise ValueError("'c' parameter can't be cast to complex_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ic = c_ptrint_t(ic)
    if __ic.value!=ic:
        raise ValueError("Error while converting 'ic' parameter to 'c_ptrint_t'")
    __jc = c_ptrint_t(jc)
    if __jc.value!=jc:
        raise ValueError("Error while converting 'jc' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        x_from_listlist(__c, c, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixgemm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__k), ctypes.byref(__alpha), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__optypea), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb), ctypes.byref(__optypeb), ctypes.byref(__beta), ctypes.byref(__c), ctypes.byref(__ic), ctypes.byref(__jc))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixgemm'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_matrix_clear(__c)


_lib_alglib.alglib_rmatrixgemm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixgemm.restype = ctypes.c_int32
def rmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc):
    pass
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ia = c_ptrint_t(ia)
    if __ia.value!=ia:
        raise ValueError("Error while converting 'ia' parameter to 'c_ptrint_t'")
    __ja = c_ptrint_t(ja)
    if __ja.value!=ja:
        raise ValueError("Error while converting 'ja' parameter to 'c_ptrint_t'")
    __optypea = c_ptrint_t(optypea)
    if __optypea.value!=optypea:
        raise ValueError("Error while converting 'optypea' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ib = c_ptrint_t(ib)
    if __ib.value!=ib:
        raise ValueError("Error while converting 'ib' parameter to 'c_ptrint_t'")
    __jb = c_ptrint_t(jb)
    if __jb.value!=jb:
        raise ValueError("Error while converting 'jb' parameter to 'c_ptrint_t'")
    __optypeb = c_ptrint_t(optypeb)
    if __optypeb.value!=optypeb:
        raise ValueError("Error while converting 'optypeb' parameter to 'c_ptrint_t'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    if not is_real_matrix(c):
        raise ValueError("'c' parameter can't be cast to real_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ic = c_ptrint_t(ic)
    if __ic.value!=ic:
        raise ValueError("Error while converting 'ic' parameter to 'c_ptrint_t'")
    __jc = c_ptrint_t(jc)
    if __jc.value!=jc:
        raise ValueError("Error while converting 'jc' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        x_from_listlist(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixgemm(ctypes.byref(_error_msg), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__k), ctypes.byref(__alpha), ctypes.byref(__a), ctypes.byref(__ia), ctypes.byref(__ja), ctypes.byref(__optypea), ctypes.byref(__b), ctypes.byref(__ib), ctypes.byref(__jb), ctypes.byref(__optypeb), ctypes.byref(__beta), ctypes.byref(__c), ctypes.byref(__ic), ctypes.byref(__jc))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixgemm'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_matrix_clear(__c)


_lib_alglib.alglib_samplemoments.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_samplemoments.restype = ctypes.c_int32
def samplemoments(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        x,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'samplemoments': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'samplemoments': function must have 1 or 2 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __mean = ctypes.c_double(0)
    __variance = ctypes.c_double(0)
    __skewness = ctypes.c_double(0)
    __kurtosis = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_samplemoments(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__mean), ctypes.byref(__variance), ctypes.byref(__skewness), ctypes.byref(__kurtosis))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'samplemoments'")
        __r__mean = __mean.value
        __r__variance = __variance.value
        __r__skewness = __skewness.value
        __r__kurtosis = __kurtosis.value
        return (__r__mean, __r__variance, __r__skewness, __r__kurtosis)
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_sampleadev.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_sampleadev.restype = ctypes.c_int32
def sampleadev(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        x,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'sampleadev': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'sampleadev': function must have 1 or 2 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __adev = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_sampleadev(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__adev))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'sampleadev'")
        __r__adev = __adev.value
        return __r__adev
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_samplemedian.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_samplemedian.restype = ctypes.c_int32
def samplemedian(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        x,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'samplemedian': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'samplemedian': function must have 1 or 2 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __median = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_samplemedian(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__median))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'samplemedian'")
        __r__median = __median.value
        return __r__median
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_samplepercentile.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_samplepercentile.restype = ctypes.c_int32
def samplepercentile(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,n,p = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,p = functionargs
        n = safe_len("'samplepercentile': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'samplepercentile': function must have 2 or 3 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_double(p)
    if __p.value!=p:
        raise ValueError("Error while converting 'p' parameter to 'ctypes.c_double'")
    __v = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_samplepercentile(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__p), ctypes.byref(__v))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'samplepercentile'")
        __r__v = __v.value
        return __r__v
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_cov2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cov2.restype = ctypes.c_int32
def cov2(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'cov2': incorrect parameters",x)!=safe_len("'cov2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'cov2': looks like one of arguments has wrong size")
        n = safe_len("'cov2': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'cov2': function must have 2 or 3 parameters")
    __result = ctypes.c_double(0)
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cov2(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cov2'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_pearsoncorr2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pearsoncorr2.restype = ctypes.c_int32
def pearsoncorr2(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'pearsoncorr2': incorrect parameters",x)!=safe_len("'pearsoncorr2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'pearsoncorr2': looks like one of arguments has wrong size")
        n = safe_len("'pearsoncorr2': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'pearsoncorr2': function must have 2 or 3 parameters")
    __result = ctypes.c_double(0)
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pearsoncorr2(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pearsoncorr2'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spearmancorr2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spearmancorr2.restype = ctypes.c_int32
def spearmancorr2(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spearmancorr2': incorrect parameters",x)!=safe_len("'spearmancorr2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spearmancorr2': looks like one of arguments has wrong size")
        n = safe_len("'spearmancorr2': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spearmancorr2': function must have 2 or 3 parameters")
    __result = ctypes.c_double(0)
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spearmancorr2(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spearmancorr2'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_covm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_covm.restype = ctypes.c_int32
def covm(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,n,m = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_rows("'covm': incorrect parameters",x)
        m = safe_cols("'covm': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'covm': function must have 1 or 3 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_covm(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'covm'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__c)


_lib_alglib.alglib_pearsoncorrm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pearsoncorrm.restype = ctypes.c_int32
def pearsoncorrm(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,n,m = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_rows("'pearsoncorrm': incorrect parameters",x)
        m = safe_cols("'pearsoncorrm': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'pearsoncorrm': function must have 1 or 3 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pearsoncorrm(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pearsoncorrm'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__c)


_lib_alglib.alglib_spearmancorrm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spearmancorrm.restype = ctypes.c_int32
def spearmancorrm(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,n,m = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_rows("'spearmancorrm': incorrect parameters",x)
        m = safe_cols("'spearmancorrm': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spearmancorrm': function must have 1 or 3 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spearmancorrm(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spearmancorrm'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__c)


_lib_alglib.alglib_covm2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_covm2.restype = ctypes.c_int32
def covm2(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        x,y,n,m1,m2 = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_rows("'covm2': incorrect parameters",x)!=safe_rows("'covm2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'covm2': looks like one of arguments has wrong size")
        n = safe_rows("'covm2': incorrect parameters",x)
        m1 = safe_cols("'covm2': incorrect parameters",x)
        m2 = safe_cols("'covm2': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'covm2': function must have 2 or 5 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(y):
        raise ValueError("'y' parameter can't be cast to real_matrix")
    __y = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m1 = c_ptrint_t(m1)
    if __m1.value!=m1:
        raise ValueError("Error while converting 'm1' parameter to 'c_ptrint_t'")
    __m2 = c_ptrint_t(m2)
    if __m2.value!=m2:
        raise ValueError("Error while converting 'm2' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_listlist(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_covm2(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m1), ctypes.byref(__m2), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'covm2'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__y)
        x_matrix_clear(__c)


_lib_alglib.alglib_pearsoncorrm2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pearsoncorrm2.restype = ctypes.c_int32
def pearsoncorrm2(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        x,y,n,m1,m2 = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_rows("'pearsoncorrm2': incorrect parameters",x)!=safe_rows("'pearsoncorrm2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'pearsoncorrm2': looks like one of arguments has wrong size")
        n = safe_rows("'pearsoncorrm2': incorrect parameters",x)
        m1 = safe_cols("'pearsoncorrm2': incorrect parameters",x)
        m2 = safe_cols("'pearsoncorrm2': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'pearsoncorrm2': function must have 2 or 5 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(y):
        raise ValueError("'y' parameter can't be cast to real_matrix")
    __y = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m1 = c_ptrint_t(m1)
    if __m1.value!=m1:
        raise ValueError("Error while converting 'm1' parameter to 'c_ptrint_t'")
    __m2 = c_ptrint_t(m2)
    if __m2.value!=m2:
        raise ValueError("Error while converting 'm2' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_listlist(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pearsoncorrm2(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m1), ctypes.byref(__m2), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pearsoncorrm2'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__y)
        x_matrix_clear(__c)


_lib_alglib.alglib_spearmancorrm2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spearmancorrm2.restype = ctypes.c_int32
def spearmancorrm2(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        x,y,n,m1,m2 = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_rows("'spearmancorrm2': incorrect parameters",x)!=safe_rows("'spearmancorrm2': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spearmancorrm2': looks like one of arguments has wrong size")
        n = safe_rows("'spearmancorrm2': incorrect parameters",x)
        m1 = safe_cols("'spearmancorrm2': incorrect parameters",x)
        m2 = safe_cols("'spearmancorrm2': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'spearmancorrm2': function must have 2 or 5 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(y):
        raise ValueError("'y' parameter can't be cast to real_matrix")
    __y = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m1 = c_ptrint_t(m1)
    if __m1.value!=m1:
        raise ValueError("Error while converting 'm1' parameter to 'c_ptrint_t'")
    __m2 = c_ptrint_t(m2)
    if __m2.value!=m2:
        raise ValueError("Error while converting 'm2' parameter to 'c_ptrint_t'")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_listlist(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spearmancorrm2(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m1), ctypes.byref(__m2), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spearmancorrm2'")
        __r__c = listlist_from_x(__c)
        return __r__c
    finally:
        x_matrix_clear(__x)
        x_matrix_clear(__y)
        x_matrix_clear(__c)


_lib_alglib.alglib_pearsoncorrelation.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pearsoncorrelation.restype = ctypes.c_int32
def pearsoncorrelation(x, y, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pearsoncorrelation(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pearsoncorrelation'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spearmanrankcorrelation.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spearmanrankcorrelation.restype = ctypes.c_int32
def spearmanrankcorrelation(x, y, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spearmanrankcorrelation(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spearmanrankcorrelation'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_dsoptimalsplit2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dsoptimalsplit2.restype = ctypes.c_int32
def dsoptimalsplit2(a, c, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(c):
        raise ValueError("'c' parameter can't be cast to int_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __threshold = ctypes.c_double(0)
    __pal = ctypes.c_double(0)
    __pbl = ctypes.c_double(0)
    __par = ctypes.c_double(0)
    __pbr = ctypes.c_double(0)
    __cve = ctypes.c_double(0)
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dsoptimalsplit2(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__threshold), ctypes.byref(__pal), ctypes.byref(__pbl), ctypes.byref(__par), ctypes.byref(__pbr), ctypes.byref(__cve))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dsoptimalsplit2'")
        __r__info = __info.value
        __r__threshold = __threshold.value
        __r__pal = __pal.value
        __r__pbl = __pbl.value
        __r__par = __par.value
        __r__pbr = __pbr.value
        __r__cve = __cve.value
        return (__r__info, __r__threshold, __r__pal, __r__pbl, __r__par, __r__pbr, __r__cve)
    finally:
        x_vector_clear(__a)
        x_vector_clear(__c)


_lib_alglib.alglib_dsoptimalsplit2fast.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dsoptimalsplit2fast.restype = ctypes.c_int32
def dsoptimalsplit2fast(a, c, tiesbuf, cntbuf, bufr, bufi, n, nc, alpha):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(c):
        raise ValueError("'c' parameter can't be cast to int_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(tiesbuf):
        raise ValueError("'tiesbuf' parameter can't be cast to int_vector")
    __tiesbuf = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(cntbuf):
        raise ValueError("'cntbuf' parameter can't be cast to int_vector")
    __cntbuf = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bufr):
        raise ValueError("'bufr' parameter can't be cast to real_vector")
    __bufr = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(bufi):
        raise ValueError("'bufi' parameter can't be cast to int_vector")
    __bufi = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nc = c_ptrint_t(nc)
    if __nc.value!=nc:
        raise ValueError("Error while converting 'nc' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __threshold = ctypes.c_double(0)
    __rms = ctypes.c_double(0)
    __cvrms = ctypes.c_double(0)
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_INT, X_CREATE)
        x_from_list(__tiesbuf, tiesbuf, DT_INT, X_CREATE)
        x_from_list(__cntbuf, cntbuf, DT_INT, X_CREATE)
        x_from_list(__bufr, bufr, DT_REAL, X_CREATE)
        x_from_list(__bufi, bufi, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dsoptimalsplit2fast(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__c), ctypes.byref(__tiesbuf), ctypes.byref(__cntbuf), ctypes.byref(__bufr), ctypes.byref(__bufi), ctypes.byref(__n), ctypes.byref(__nc), ctypes.byref(__alpha), ctypes.byref(__info), ctypes.byref(__threshold), ctypes.byref(__rms), ctypes.byref(__cvrms))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dsoptimalsplit2fast'")
        __r__a = list_from_x(__a)
        __r__c = list_from_x(__c)
        __r__tiesbuf = list_from_x(__tiesbuf)
        __r__cntbuf = list_from_x(__cntbuf)
        __r__bufr = list_from_x(__bufr)
        __r__bufi = list_from_x(__bufi)
        __r__info = __info.value
        __r__threshold = __threshold.value
        __r__rms = __rms.value
        __r__cvrms = __cvrms.value
        return (__r__a, __r__c, __r__tiesbuf, __r__cntbuf, __r__bufr, __r__bufi, __r__info, __r__threshold, __r__rms, __r__cvrms)
    finally:
        x_vector_clear(__a)
        x_vector_clear(__c)
        x_vector_clear(__tiesbuf)
        x_vector_clear(__cntbuf)
        x_vector_clear(__bufr)
        x_vector_clear(__bufi)


_lib_alglib.x_obj_free_decisionforest.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_decisionforest.restype = None


class decisionforest(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_decisionforest(self.ptr)


class x_dfreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("relclserror", ctypes.c_double),
        ("avgce", ctypes.c_double),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("oobrelclserror", ctypes.c_double),
        ("oobavgce", ctypes.c_double),
        ("oobrmserror", ctypes.c_double),
        ("oobavgerror", ctypes.c_double),
        ("oobavgrelerror", ctypes.c_double)
        ]




class dfreport(object):
    def __init__(self):
        self.relclserror = 0
        self.avgce = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.oobrelclserror = 0
        self.oobavgce = 0
        self.oobrmserror = 0
        self.oobavgerror = 0
        self.oobavgrelerror = 0


def x_dfreport_zero_fields(x):
    x.relclserror = 0
    x.avgce = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.oobrelclserror = 0
    x.oobavgce = 0
    x.oobrmserror = 0
    x.oobavgerror = 0
    x.oobavgrelerror = 0
    return




def x_dfreport_clear(x):
    x_dfreport_zero_fields(x)
    return




def x_from_dfreport(x,v):
    x.relclserror = float(v.relclserror)
    x.avgce = float(v.avgce)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.oobrelclserror = float(v.oobrelclserror)
    x.oobavgce = float(v.oobavgce)
    x.oobrmserror = float(v.oobrmserror)
    x.oobavgerror = float(v.oobavgerror)
    x.oobavgrelerror = float(v.oobavgrelerror)
    return




def dfreport_from_x(x):
    r = dfreport()
    r.relclserror = x.relclserror
    r.avgce = x.avgce
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.oobrelclserror = x.oobrelclserror
    r.oobavgce = x.oobavgce
    r.oobrmserror = x.oobrmserror
    r.oobavgerror = x.oobavgerror
    r.oobavgrelerror = x.oobavgrelerror
    return r


_lib_alglib.alglib_dfbuildrandomdecisionforest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfbuildrandomdecisionforest.restype = ctypes.c_int32
def dfbuildrandomdecisionforest(xy, npoints, nvars, nclasses, ntrees, r):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __nclasses = c_ptrint_t(nclasses)
    if __nclasses.value!=nclasses:
        raise ValueError("Error while converting 'nclasses' parameter to 'c_ptrint_t'")
    __ntrees = c_ptrint_t(ntrees)
    if __ntrees.value!=ntrees:
        raise ValueError("Error while converting 'ntrees' parameter to 'c_ptrint_t'")
    __r = ctypes.c_double(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __df = ctypes.c_void_p(0)
    __rep = x_dfreport()
    x_dfreport_zero_fields(__rep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfbuildrandomdecisionforest(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__nclasses), ctypes.byref(__ntrees), ctypes.byref(__r), ctypes.byref(__info), ctypes.byref(__df), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfbuildrandomdecisionforest'")
        __r__info = __info.value
        __r__df = decisionforest(__df)
        __r__rep = dfreport_from_x(__rep)
        return (__r__info, __r__df, __r__rep)
    finally:
        x_matrix_clear(__xy)
        x_dfreport_clear(__rep)


_lib_alglib.alglib_dfprocess.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfprocess.restype = ctypes.c_int32
def dfprocess(df, x, y):
    pass
    __df = df.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfprocess(ctypes.byref(_error_msg), ctypes.byref(__df), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfprocess'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_dfprocessi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfprocessi.restype = ctypes.c_int32
def dfprocessi(df, x):
    pass
    __df = df.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __y = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfprocessi(ctypes.byref(_error_msg), ctypes.byref(__df), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfprocessi'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_dfrelclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfrelclserror.restype = ctypes.c_int32
def dfrelclserror(df, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __df = df.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfrelclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__df), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfrelclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_dfavgce.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfavgce.restype = ctypes.c_int32
def dfavgce(df, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __df = df.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfavgce(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__df), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfavgce'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_dfrmserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfrmserror.restype = ctypes.c_int32
def dfrmserror(df, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __df = df.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfrmserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__df), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfrmserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_dfavgerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfavgerror.restype = ctypes.c_int32
def dfavgerror(df, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __df = df.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfavgerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__df), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfavgerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_dfavgrelerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dfavgrelerror.restype = ctypes.c_int32
def dfavgrelerror(df, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __df = df.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dfavgrelerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__df), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dfavgrelerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_kmeansgenerate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_kmeansgenerate.restype = ctypes.c_int32
def kmeansgenerate(xy, npoints, nvars, k, restarts):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __xyc = x_vector(cnt=0,datatype=DT_INT,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_kmeansgenerate(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__k), ctypes.byref(__restarts), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__xyc))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'kmeansgenerate'")
        __r__info = __info.value
        __r__c = listlist_from_x(__c)
        __r__xyc = list_from_x(__xyc)
        return (__r__info, __r__c, __r__xyc)
    finally:
        x_matrix_clear(__xy)
        x_matrix_clear(__c)
        x_vector_clear(__xyc)


_lib_alglib.alglib_rmatrixqr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixqr.restype = ctypes.c_int32
def rmatrixqr(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tau = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixqr(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixqr'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        return (__r__a, __r__tau)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)


_lib_alglib.alglib_rmatrixlq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlq.restype = ctypes.c_int32
def rmatrixlq(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tau = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlq'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        return (__r__a, __r__tau)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)


_lib_alglib.alglib_cmatrixqr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixqr.restype = ctypes.c_int32
def cmatrixqr(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tau = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixqr(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixqr'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        return (__r__a, __r__tau)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)


_lib_alglib.alglib_cmatrixlq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlq.restype = ctypes.c_int32
def cmatrixlq(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tau = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlq'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        return (__r__a, __r__tau)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)


_lib_alglib.alglib_rmatrixqrunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixqrunpackq.restype = ctypes.c_int32
def rmatrixqrunpackq(a, m, n, tau, qcolumns):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(tau):
        raise ValueError("'tau' parameter can't be cast to real_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __qcolumns = c_ptrint_t(qcolumns)
    if __qcolumns.value!=qcolumns:
        raise ValueError("Error while converting 'qcolumns' parameter to 'c_ptrint_t'")
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__tau, tau, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixqrunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau), ctypes.byref(__qcolumns), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixqrunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_rmatrixqrunpackr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixqrunpackr.restype = ctypes.c_int32
def rmatrixqrunpackr(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixqrunpackr(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixqrunpackr'")
        __r__r = listlist_from_x(__r)
        return __r__r
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__r)


_lib_alglib.alglib_rmatrixlqunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlqunpackq.restype = ctypes.c_int32
def rmatrixlqunpackq(a, m, n, tau, qrows):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(tau):
        raise ValueError("'tau' parameter can't be cast to real_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __qrows = c_ptrint_t(qrows)
    if __qrows.value!=qrows:
        raise ValueError("Error while converting 'qrows' parameter to 'c_ptrint_t'")
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__tau, tau, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlqunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau), ctypes.byref(__qrows), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlqunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_rmatrixlqunpackl.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlqunpackl.restype = ctypes.c_int32
def rmatrixlqunpackl(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __l = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlqunpackl(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__l))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlqunpackl'")
        __r__l = listlist_from_x(__l)
        return __r__l
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__l)


_lib_alglib.alglib_cmatrixqrunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixqrunpackq.restype = ctypes.c_int32
def cmatrixqrunpackq(a, m, n, tau, qcolumns):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(tau):
        raise ValueError("'tau' parameter can't be cast to complex_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __qcolumns = c_ptrint_t(qcolumns)
    if __qcolumns.value!=qcolumns:
        raise ValueError("Error while converting 'qcolumns' parameter to 'c_ptrint_t'")
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__tau, tau, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixqrunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau), ctypes.byref(__qcolumns), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixqrunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_cmatrixqrunpackr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixqrunpackr.restype = ctypes.c_int32
def cmatrixqrunpackr(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixqrunpackr(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixqrunpackr'")
        __r__r = listlist_from_x(__r)
        return __r__r
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__r)


_lib_alglib.alglib_cmatrixlqunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlqunpackq.restype = ctypes.c_int32
def cmatrixlqunpackq(a, m, n, tau, qrows):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(tau):
        raise ValueError("'tau' parameter can't be cast to complex_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __qrows = c_ptrint_t(qrows)
    if __qrows.value!=qrows:
        raise ValueError("Error while converting 'qrows' parameter to 'c_ptrint_t'")
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__tau, tau, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlqunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tau), ctypes.byref(__qrows), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlqunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_cmatrixlqunpackl.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlqunpackl.restype = ctypes.c_int32
def cmatrixlqunpackl(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __l = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlqunpackl(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__l))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlqunpackl'")
        __r__l = listlist_from_x(__l)
        return __r__l
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__l)


_lib_alglib.alglib_rmatrixbd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbd.restype = ctypes.c_int32
def rmatrixbd(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tauq = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __taup = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbd(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tauq), ctypes.byref(__taup))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbd'")
        __r__a = listlist_from_x(__a)
        __r__tauq = list_from_x(__tauq)
        __r__taup = list_from_x(__taup)
        return (__r__a, __r__tauq, __r__taup)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tauq)
        x_vector_clear(__taup)


_lib_alglib.alglib_rmatrixbdunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdunpackq.restype = ctypes.c_int32
def rmatrixbdunpackq(qp, m, n, tauq, qcolumns):
    pass
    if not is_real_matrix(qp):
        raise ValueError("'qp' parameter can't be cast to real_matrix")
    __qp = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(tauq):
        raise ValueError("'tauq' parameter can't be cast to real_vector")
    __tauq = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __qcolumns = c_ptrint_t(qcolumns)
    if __qcolumns.value!=qcolumns:
        raise ValueError("Error while converting 'qcolumns' parameter to 'c_ptrint_t'")
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__qp, qp, DT_REAL, X_CREATE)
        x_from_list(__tauq, tauq, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdunpackq(ctypes.byref(_error_msg), ctypes.byref(__qp), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tauq), ctypes.byref(__qcolumns), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__qp)
        x_vector_clear(__tauq)
        x_matrix_clear(__q)


_lib_alglib.alglib_rmatrixbdmultiplybyq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdmultiplybyq.restype = ctypes.c_int32
def rmatrixbdmultiplybyq(qp, m, n, tauq, z, zrows, zcolumns, fromtheright, dotranspose):
    pass
    if not is_real_matrix(qp):
        raise ValueError("'qp' parameter can't be cast to real_matrix")
    __qp = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(tauq):
        raise ValueError("'tauq' parameter can't be cast to real_vector")
    __tauq = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(z):
        raise ValueError("'z' parameter can't be cast to real_matrix")
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __zrows = c_ptrint_t(zrows)
    if __zrows.value!=zrows:
        raise ValueError("Error while converting 'zrows' parameter to 'c_ptrint_t'")
    __zcolumns = c_ptrint_t(zcolumns)
    if __zcolumns.value!=zcolumns:
        raise ValueError("Error while converting 'zcolumns' parameter to 'c_ptrint_t'")
    __fromtheright = ctypes.c_uint8(fromtheright)
    if __fromtheright.value!=0:
        __fromtheright = ctypes.c_uint8(1)
    __dotranspose = ctypes.c_uint8(dotranspose)
    if __dotranspose.value!=0:
        __dotranspose = ctypes.c_uint8(1)
    try:
        x_from_listlist(__qp, qp, DT_REAL, X_CREATE)
        x_from_list(__tauq, tauq, DT_REAL, X_CREATE)
        x_from_listlist(__z, z, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdmultiplybyq(ctypes.byref(_error_msg), ctypes.byref(__qp), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tauq), ctypes.byref(__z), ctypes.byref(__zrows), ctypes.byref(__zcolumns), ctypes.byref(__fromtheright), ctypes.byref(__dotranspose))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdmultiplybyq'")
        __r__z = listlist_from_x(__z)
        return __r__z
    finally:
        x_matrix_clear(__qp)
        x_vector_clear(__tauq)
        x_matrix_clear(__z)


_lib_alglib.alglib_rmatrixbdunpackpt.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdunpackpt.restype = ctypes.c_int32
def rmatrixbdunpackpt(qp, m, n, taup, ptrows):
    pass
    if not is_real_matrix(qp):
        raise ValueError("'qp' parameter can't be cast to real_matrix")
    __qp = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(taup):
        raise ValueError("'taup' parameter can't be cast to real_vector")
    __taup = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ptrows = c_ptrint_t(ptrows)
    if __ptrows.value!=ptrows:
        raise ValueError("Error while converting 'ptrows' parameter to 'c_ptrint_t'")
    __pt = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__qp, qp, DT_REAL, X_CREATE)
        x_from_list(__taup, taup, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdunpackpt(ctypes.byref(_error_msg), ctypes.byref(__qp), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__taup), ctypes.byref(__ptrows), ctypes.byref(__pt))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdunpackpt'")
        __r__pt = listlist_from_x(__pt)
        return __r__pt
    finally:
        x_matrix_clear(__qp)
        x_vector_clear(__taup)
        x_matrix_clear(__pt)


_lib_alglib.alglib_rmatrixbdmultiplybyp.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdmultiplybyp.restype = ctypes.c_int32
def rmatrixbdmultiplybyp(qp, m, n, taup, z, zrows, zcolumns, fromtheright, dotranspose):
    pass
    if not is_real_matrix(qp):
        raise ValueError("'qp' parameter can't be cast to real_matrix")
    __qp = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(taup):
        raise ValueError("'taup' parameter can't be cast to real_vector")
    __taup = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(z):
        raise ValueError("'z' parameter can't be cast to real_matrix")
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __zrows = c_ptrint_t(zrows)
    if __zrows.value!=zrows:
        raise ValueError("Error while converting 'zrows' parameter to 'c_ptrint_t'")
    __zcolumns = c_ptrint_t(zcolumns)
    if __zcolumns.value!=zcolumns:
        raise ValueError("Error while converting 'zcolumns' parameter to 'c_ptrint_t'")
    __fromtheright = ctypes.c_uint8(fromtheright)
    if __fromtheright.value!=0:
        __fromtheright = ctypes.c_uint8(1)
    __dotranspose = ctypes.c_uint8(dotranspose)
    if __dotranspose.value!=0:
        __dotranspose = ctypes.c_uint8(1)
    try:
        x_from_listlist(__qp, qp, DT_REAL, X_CREATE)
        x_from_list(__taup, taup, DT_REAL, X_CREATE)
        x_from_listlist(__z, z, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdmultiplybyp(ctypes.byref(_error_msg), ctypes.byref(__qp), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__taup), ctypes.byref(__z), ctypes.byref(__zrows), ctypes.byref(__zcolumns), ctypes.byref(__fromtheright), ctypes.byref(__dotranspose))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdmultiplybyp'")
        __r__z = listlist_from_x(__z)
        return __r__z
    finally:
        x_matrix_clear(__qp)
        x_vector_clear(__taup)
        x_matrix_clear(__z)


_lib_alglib.alglib_rmatrixbdunpackdiagonals.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdunpackdiagonals.restype = ctypes.c_int32
def rmatrixbdunpackdiagonals(b, m, n):
    pass
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(0)
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __e = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdunpackdiagonals(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__d), ctypes.byref(__e))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdunpackdiagonals'")
        __r__isupper = __isupper.value!=0
        __r__d = list_from_x(__d)
        __r__e = list_from_x(__e)
        return (__r__isupper, __r__d, __r__e)
    finally:
        x_matrix_clear(__b)
        x_vector_clear(__d)
        x_vector_clear(__e)


_lib_alglib.alglib_rmatrixhessenberg.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixhessenberg.restype = ctypes.c_int32
def rmatrixhessenberg(a, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __tau = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixhessenberg(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__tau))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixhessenberg'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        return (__r__a, __r__tau)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)


_lib_alglib.alglib_rmatrixhessenbergunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixhessenbergunpackq.restype = ctypes.c_int32
def rmatrixhessenbergunpackq(a, n, tau):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(tau):
        raise ValueError("'tau' parameter can't be cast to real_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__tau, tau, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixhessenbergunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__tau), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixhessenbergunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_rmatrixhessenbergunpackh.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixhessenbergunpackh.restype = ctypes.c_int32
def rmatrixhessenbergunpackh(a, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __h = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixhessenbergunpackh(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__h))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixhessenbergunpackh'")
        __r__h = listlist_from_x(__h)
        return __r__h
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__h)


_lib_alglib.alglib_smatrixtd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixtd.restype = ctypes.c_int32
def smatrixtd(a, n, isupper):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __tau = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __e = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixtd(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__tau), ctypes.byref(__d), ctypes.byref(__e))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixtd'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        __r__d = list_from_x(__d)
        __r__e = list_from_x(__e)
        return (__r__a, __r__tau, __r__d, __r__e)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_vector_clear(__d)
        x_vector_clear(__e)


_lib_alglib.alglib_smatrixtdunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixtdunpackq.restype = ctypes.c_int32
def smatrixtdunpackq(a, n, isupper, tau):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_real_vector(tau):
        raise ValueError("'tau' parameter can't be cast to real_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__tau, tau, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixtdunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__tau), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixtdunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_hmatrixtd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixtd.restype = ctypes.c_int32
def hmatrixtd(a, n, isupper):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __tau = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __e = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixtd(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__tau), ctypes.byref(__d), ctypes.byref(__e))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixtd'")
        __r__a = listlist_from_x(__a)
        __r__tau = list_from_x(__tau)
        __r__d = list_from_x(__d)
        __r__e = list_from_x(__e)
        return (__r__a, __r__tau, __r__d, __r__e)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_vector_clear(__d)
        x_vector_clear(__e)


_lib_alglib.alglib_hmatrixtdunpackq.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixtdunpackq.restype = ctypes.c_int32
def hmatrixtdunpackq(a, n, isupper, tau):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_complex_vector(tau):
        raise ValueError("'tau' parameter can't be cast to complex_vector")
    __tau = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __q = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__tau, tau, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixtdunpackq(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__tau), ctypes.byref(__q))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixtdunpackq'")
        __r__q = listlist_from_x(__q)
        return __r__q
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__tau)
        x_matrix_clear(__q)


_lib_alglib.alglib_smatrixevd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixevd.restype = ctypes.c_int32
def smatrixevd(a, n, zneeded, isupper):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixevd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__d), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixevd'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__d)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixevdr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixevdr.restype = ctypes.c_int32
def smatrixevdr(a, n, zneeded, isupper, b1, b2):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __b1 = ctypes.c_double(b1)
    if __b1.value!=b1:
        raise ValueError("Error while converting 'b1' parameter to 'ctypes.c_double'")
    __b2 = ctypes.c_double(b2)
    if __b2.value!=b2:
        raise ValueError("Error while converting 'b2' parameter to 'ctypes.c_double'")
    __m = c_ptrint_t(0)
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixevdr(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__b1), ctypes.byref(__b2), ctypes.byref(__m), ctypes.byref(__w), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixevdr'")
        __r__result = __result.value!=0
        __r__m = __m.value
        __r__w = list_from_x(__w)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__m, __r__w, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__w)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixevdi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixevdi.restype = ctypes.c_int32
def smatrixevdi(a, n, zneeded, isupper, i1, i2):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixevdi(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__i1), ctypes.byref(__i2), ctypes.byref(__w), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixevdi'")
        __r__result = __result.value!=0
        __r__w = list_from_x(__w)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__w, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__w)
        x_matrix_clear(__z)


_lib_alglib.alglib_hmatrixevd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixevd.restype = ctypes.c_int32
def hmatrixevd(a, n, zneeded, isupper):
    pass
    __result = ctypes.c_uint8(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixevd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__d), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixevd'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__d)
        x_matrix_clear(__z)


_lib_alglib.alglib_hmatrixevdr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixevdr.restype = ctypes.c_int32
def hmatrixevdr(a, n, zneeded, isupper, b1, b2):
    pass
    __result = ctypes.c_uint8(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __b1 = ctypes.c_double(b1)
    if __b1.value!=b1:
        raise ValueError("Error while converting 'b1' parameter to 'ctypes.c_double'")
    __b2 = ctypes.c_double(b2)
    if __b2.value!=b2:
        raise ValueError("Error while converting 'b2' parameter to 'ctypes.c_double'")
    __m = c_ptrint_t(0)
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixevdr(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__b1), ctypes.byref(__b2), ctypes.byref(__m), ctypes.byref(__w), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixevdr'")
        __r__result = __result.value!=0
        __r__m = __m.value
        __r__w = list_from_x(__w)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__m, __r__w, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__w)
        x_matrix_clear(__z)


_lib_alglib.alglib_hmatrixevdi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixevdi.restype = ctypes.c_int32
def hmatrixevdi(a, n, zneeded, isupper, i1, i2):
    pass
    __result = ctypes.c_uint8(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixevdi(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__isupper), ctypes.byref(__i1), ctypes.byref(__i2), ctypes.byref(__w), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixevdi'")
        __r__result = __result.value!=0
        __r__w = list_from_x(__w)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__w, __r__z)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__w)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixtdevd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixtdevd.restype = ctypes.c_int32
def smatrixtdevd(d, e, n, zneeded, z):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_vector(d):
        raise ValueError("'d' parameter can't be cast to real_vector")
    __d = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(e):
        raise ValueError("'e' parameter can't be cast to real_vector")
    __e = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    if not is_real_matrix(z):
        raise ValueError("'z' parameter can't be cast to real_matrix")
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__d, d, DT_REAL, X_CREATE)
        x_from_list(__e, e, DT_REAL, X_CREATE)
        x_from_listlist(__z, z, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixtdevd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__d), ctypes.byref(__e), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixtdevd'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__z)
    finally:
        x_vector_clear(__d)
        x_vector_clear(__e)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixtdevdr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixtdevdr.restype = ctypes.c_int32
def smatrixtdevdr(d, e, n, zneeded, a, b, z):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_vector(d):
        raise ValueError("'d' parameter can't be cast to real_vector")
    __d = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(e):
        raise ValueError("'e' parameter can't be cast to real_vector")
    __e = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __m = c_ptrint_t(0)
    if not is_real_matrix(z):
        raise ValueError("'z' parameter can't be cast to real_matrix")
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__d, d, DT_REAL, X_CREATE)
        x_from_list(__e, e, DT_REAL, X_CREATE)
        x_from_listlist(__z, z, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixtdevdr(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__d), ctypes.byref(__e), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixtdevdr'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__m = __m.value
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__m, __r__z)
    finally:
        x_vector_clear(__d)
        x_vector_clear(__e)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixtdevdi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixtdevdi.restype = ctypes.c_int32
def smatrixtdevdi(d, e, n, zneeded, i1, i2, z):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_vector(d):
        raise ValueError("'d' parameter can't be cast to real_vector")
    __d = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(e):
        raise ValueError("'e' parameter can't be cast to real_vector")
    __e = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __i1 = c_ptrint_t(i1)
    if __i1.value!=i1:
        raise ValueError("Error while converting 'i1' parameter to 'c_ptrint_t'")
    __i2 = c_ptrint_t(i2)
    if __i2.value!=i2:
        raise ValueError("Error while converting 'i2' parameter to 'c_ptrint_t'")
    if not is_real_matrix(z):
        raise ValueError("'z' parameter can't be cast to real_matrix")
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__d, d, DT_REAL, X_CREATE)
        x_from_list(__e, e, DT_REAL, X_CREATE)
        x_from_listlist(__z, z, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixtdevdi(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__d), ctypes.byref(__e), ctypes.byref(__n), ctypes.byref(__zneeded), ctypes.byref(__i1), ctypes.byref(__i2), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixtdevdi'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__z)
    finally:
        x_vector_clear(__d)
        x_vector_clear(__e)
        x_matrix_clear(__z)


_lib_alglib.alglib_rmatrixevd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixevd.restype = ctypes.c_int32
def rmatrixevd(a, n, vneeded):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __vneeded = c_ptrint_t(vneeded)
    if __vneeded.value!=vneeded:
        raise ValueError("Error while converting 'vneeded' parameter to 'c_ptrint_t'")
    __wr = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wi = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __vl = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __vr = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixevd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__vneeded), ctypes.byref(__wr), ctypes.byref(__wi), ctypes.byref(__vl), ctypes.byref(__vr))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixevd'")
        __r__result = __result.value!=0
        __r__wr = list_from_x(__wr)
        __r__wi = list_from_x(__wi)
        __r__vl = listlist_from_x(__vl)
        __r__vr = listlist_from_x(__vr)
        return (__r__result, __r__wr, __r__wi, __r__vl, __r__vr)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__wr)
        x_vector_clear(__wi)
        x_matrix_clear(__vl)
        x_matrix_clear(__vr)


_lib_alglib.alglib_rmatrixrndorthogonal.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrndorthogonal.restype = ctypes.c_int32
def rmatrixrndorthogonal(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrndorthogonal(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrndorthogonal'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrndcond.restype = ctypes.c_int32
def rmatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrndorthogonal.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrndorthogonal.restype = ctypes.c_int32
def cmatrixrndorthogonal(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrndorthogonal(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrndorthogonal'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrndcond.restype = ctypes.c_int32
def cmatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_smatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixrndcond.restype = ctypes.c_int32
def smatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_spdmatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixrndcond.restype = ctypes.c_int32
def spdmatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_hmatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixrndcond.restype = ctypes.c_int32
def hmatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_hpdmatrixrndcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixrndcond.restype = ctypes.c_int32
def hpdmatrixrndcond(n, c):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixrndcond(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixrndcond'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixrndorthogonalfromtheright.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrndorthogonalfromtheright.restype = ctypes.c_int32
def rmatrixrndorthogonalfromtheright(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrndorthogonalfromtheright(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrndorthogonalfromtheright'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixrndorthogonalfromtheleft.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrndorthogonalfromtheleft.restype = ctypes.c_int32
def rmatrixrndorthogonalfromtheleft(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrndorthogonalfromtheleft(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrndorthogonalfromtheleft'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrndorthogonalfromtheright.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrndorthogonalfromtheright.restype = ctypes.c_int32
def cmatrixrndorthogonalfromtheright(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrndorthogonalfromtheright(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrndorthogonalfromtheright'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrndorthogonalfromtheleft.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrndorthogonalfromtheleft.restype = ctypes.c_int32
def cmatrixrndorthogonalfromtheleft(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrndorthogonalfromtheleft(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrndorthogonalfromtheleft'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_smatrixrndmultiply.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixrndmultiply.restype = ctypes.c_int32
def smatrixrndmultiply(a, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixrndmultiply(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixrndmultiply'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_hmatrixrndmultiply.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hmatrixrndmultiply.restype = ctypes.c_int32
def hmatrixrndmultiply(a, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hmatrixrndmultiply(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hmatrixrndmultiply'")
        __r__a = listlist_from_x(__a)
        return __r__a
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixlu.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlu.restype = ctypes.c_int32
def rmatrixlu(a, m, n):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __pivots = x_vector(cnt=0,datatype=DT_INT,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlu(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__pivots))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlu'")
        __r__a = listlist_from_x(__a)
        __r__pivots = list_from_x(__pivots)
        return (__r__a, __r__pivots)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)


_lib_alglib.alglib_cmatrixlu.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlu.restype = ctypes.c_int32
def cmatrixlu(a, m, n):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __pivots = x_vector(cnt=0,datatype=DT_INT,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlu(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__pivots))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlu'")
        __r__a = listlist_from_x(__a)
        __r__pivots = list_from_x(__pivots)
        return (__r__a, __r__pivots)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)


_lib_alglib.alglib_hpdmatrixcholesky.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixcholesky.restype = ctypes.c_int32
def hpdmatrixcholesky(a, n, isupper):
    pass
    __result = ctypes.c_uint8(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixcholesky(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixcholesky'")
        __r__result = __result.value!=0
        __r__a = listlist_from_x(__a)
        return (__r__result, __r__a)
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_spdmatrixcholesky.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholesky.restype = ctypes.c_int32
def spdmatrixcholesky(a, n, isupper):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholesky(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholesky'")
        __r__result = __result.value!=0
        __r__a = listlist_from_x(__a)
        return (__r__result, __r__a)
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixrcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrcond1.restype = ctypes.c_int32
def rmatrixrcond1(a, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixrcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixrcondinf.restype = ctypes.c_int32
def rmatrixrcondinf(a, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixrcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixrcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_spdmatrixrcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixrcond.restype = ctypes.c_int32
def spdmatrixrcond(a, n, isupper):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixrcond(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixrcond'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixtrrcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixtrrcond1.restype = ctypes.c_int32
def rmatrixtrrcond1(a, n, isupper, isunit):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixtrrcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixtrrcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixtrrcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixtrrcondinf.restype = ctypes.c_int32
def rmatrixtrrcondinf(a, n, isupper, isunit):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixtrrcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixtrrcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_hpdmatrixrcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixrcond.restype = ctypes.c_int32
def hpdmatrixrcond(a, n, isupper):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixrcond(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixrcond'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrcond1.restype = ctypes.c_int32
def cmatrixrcond1(a, n):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixrcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixrcondinf.restype = ctypes.c_int32
def cmatrixrcondinf(a, n):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixrcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixrcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_rmatrixlurcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlurcond1.restype = ctypes.c_int32
def rmatrixlurcond1(lua, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlurcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lua), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlurcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__lua)


_lib_alglib.alglib_rmatrixlurcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlurcondinf.restype = ctypes.c_int32
def rmatrixlurcondinf(lua, n):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlurcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lua), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlurcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__lua)


_lib_alglib.alglib_spdmatrixcholeskyrcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholeskyrcond.restype = ctypes.c_int32
def spdmatrixcholeskyrcond(a, n, isupper):
    pass
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholeskyrcond(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholeskyrcond'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_hpdmatrixcholeskyrcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixcholeskyrcond.restype = ctypes.c_int32
def hpdmatrixcholeskyrcond(a, n, isupper):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixcholeskyrcond(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixcholeskyrcond'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixlurcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlurcond1.restype = ctypes.c_int32
def cmatrixlurcond1(lua, n):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlurcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lua), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlurcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__lua)


_lib_alglib.alglib_cmatrixlurcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlurcondinf.restype = ctypes.c_int32
def cmatrixlurcondinf(lua, n):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlurcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lua), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlurcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__lua)


_lib_alglib.alglib_cmatrixtrrcond1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixtrrcond1.restype = ctypes.c_int32
def cmatrixtrrcond1(a, n, isupper, isunit):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixtrrcond1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixtrrcond1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixtrrcondinf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixtrrcondinf.restype = ctypes.c_int32
def cmatrixtrrcondinf(a, n, isupper, isunit):
    pass
    __result = ctypes.c_double(0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixtrrcondinf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixtrrcondinf'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)




class x_matinvreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("r1", ctypes.c_double),
        ("rinf", ctypes.c_double)
        ]




class matinvreport(object):
    def __init__(self):
        self.r1 = 0
        self.rinf = 0


def x_matinvreport_zero_fields(x):
    x.r1 = 0
    x.rinf = 0
    return




def x_matinvreport_clear(x):
    x_matinvreport_zero_fields(x)
    return




def x_from_matinvreport(x,v):
    x.r1 = float(v.r1)
    x.rinf = float(v.rinf)
    return




def matinvreport_from_x(x):
    r = matinvreport()
    r.r1 = x.r1
    r.rinf = x.rinf
    return r


_lib_alglib.alglib_rmatrixluinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixluinverse.restype = ctypes.c_int32
def rmatrixluinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,pivots,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,pivots = functionargs
        if safe_cols("'rmatrixluinverse': incorrect parameters",a)!=safe_rows("'rmatrixluinverse': incorrect parameters",a) or safe_cols("'rmatrixluinverse': incorrect parameters",a)!=safe_len("'rmatrixluinverse': incorrect parameters",pivots):
            raise RuntimeError("Error while calling 'rmatrixluinverse': looks like one of arguments has wrong size")
        n = safe_cols("'rmatrixluinverse': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'rmatrixluinverse': function must have 2 or 3 parameters")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(pivots):
        raise ValueError("'pivots' parameter can't be cast to int_vector")
    __pivots = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__pivots, pivots, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixluinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__pivots), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixluinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_rmatrixinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixinverse.restype = ctypes.c_int32
def rmatrixinverse(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'rmatrixinverse': incorrect parameters",a)!=safe_rows("'rmatrixinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'rmatrixinverse': looks like one of arguments has wrong size")
        n = safe_cols("'rmatrixinverse': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'rmatrixinverse': function must have 1 or 2 parameters")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_cmatrixluinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixluinverse.restype = ctypes.c_int32
def cmatrixluinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,pivots,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,pivots = functionargs
        if safe_cols("'cmatrixluinverse': incorrect parameters",a)!=safe_rows("'cmatrixluinverse': incorrect parameters",a) or safe_cols("'cmatrixluinverse': incorrect parameters",a)!=safe_len("'cmatrixluinverse': incorrect parameters",pivots):
            raise RuntimeError("Error while calling 'cmatrixluinverse': looks like one of arguments has wrong size")
        n = safe_cols("'cmatrixluinverse': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'cmatrixluinverse': function must have 2 or 3 parameters")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(pivots):
        raise ValueError("'pivots' parameter can't be cast to int_vector")
    __pivots = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__pivots, pivots, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixluinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__pivots), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixluinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_cmatrixinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixinverse.restype = ctypes.c_int32
def cmatrixinverse(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'cmatrixinverse': incorrect parameters",a)!=safe_rows("'cmatrixinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'cmatrixinverse': looks like one of arguments has wrong size")
        n = safe_cols("'cmatrixinverse': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'cmatrixinverse': function must have 1 or 2 parameters")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_spdmatrixcholeskyinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholeskyinverse.restype = ctypes.c_int32
def spdmatrixcholeskyinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,n,isupper = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'spdmatrixcholeskyinverse': incorrect parameters",a)!=safe_rows("'spdmatrixcholeskyinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'spdmatrixcholeskyinverse': looks like one of arguments has wrong size")
        n = safe_cols("'spdmatrixcholeskyinverse': incorrect parameters",a)
        isupper = False
    else:
        raise RuntimeError("Error while calling 'spdmatrixcholeskyinverse': function must have 1 or 3 parameters")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholeskyinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholeskyinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_spdmatrixinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixinverse.restype = ctypes.c_int32
def spdmatrixinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,n,isupper = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'spdmatrixinverse': incorrect parameters",a)!=safe_rows("'spdmatrixinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'spdmatrixinverse': looks like one of arguments has wrong size")
        n = safe_cols("'spdmatrixinverse': incorrect parameters",a)
        isupper = False
    else:
        raise RuntimeError("Error while calling 'spdmatrixinverse': function must have 1 or 3 parameters")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to symmetric")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        if __friendly_form:
            if not x_is_symmetric(__a):
                raise ValueError("'a' parameter is not symmetric matrix")
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixinverse'")
        if __friendly_form:
            if not x_force_symmetric(__a):
                raise RuntimeError("Internal error while forcing symmetricity of 'a' parameter")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_hpdmatrixcholeskyinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixcholeskyinverse.restype = ctypes.c_int32
def hpdmatrixcholeskyinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,n,isupper = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'hpdmatrixcholeskyinverse': incorrect parameters",a)!=safe_rows("'hpdmatrixcholeskyinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'hpdmatrixcholeskyinverse': looks like one of arguments has wrong size")
        n = safe_cols("'hpdmatrixcholeskyinverse': incorrect parameters",a)
        isupper = False
    else:
        raise RuntimeError("Error while calling 'hpdmatrixcholeskyinverse': function must have 1 or 3 parameters")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixcholeskyinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixcholeskyinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_hpdmatrixinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixinverse.restype = ctypes.c_int32
def hpdmatrixinverse(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,n,isupper = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_cols("'hpdmatrixinverse': incorrect parameters",a)!=safe_rows("'hpdmatrixinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'hpdmatrixinverse': looks like one of arguments has wrong size")
        n = safe_cols("'hpdmatrixinverse': incorrect parameters",a)
        isupper = False
    else:
        raise RuntimeError("Error while calling 'hpdmatrixinverse': function must have 1 or 3 parameters")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to hermitian")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        if __friendly_form:
            if not x_is_hermitian(__a):
                raise ValueError("'a' parameter is not Hermitian matrix")
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixinverse'")
        if __friendly_form:
            if not x_force_hermitian(__a):
                raise RuntimeError("Internal error while forcing Hermitian properties of 'a' parameter")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_rmatrixtrinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixtrinverse.restype = ctypes.c_int32
def rmatrixtrinverse(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,n,isupper,isunit = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,isupper = functionargs
        if safe_cols("'rmatrixtrinverse': incorrect parameters",a)!=safe_rows("'rmatrixtrinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'rmatrixtrinverse': looks like one of arguments has wrong size")
        n = safe_cols("'rmatrixtrinverse': incorrect parameters",a)
        isunit = False
    else:
        raise RuntimeError("Error while calling 'rmatrixtrinverse': function must have 2 or 4 parameters")
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixtrinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixtrinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_cmatrixtrinverse.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixtrinverse.restype = ctypes.c_int32
def cmatrixtrinverse(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,n,isupper,isunit = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,isupper = functionargs
        if safe_cols("'cmatrixtrinverse': incorrect parameters",a)!=safe_rows("'cmatrixtrinverse': incorrect parameters",a):
            raise RuntimeError("Error while calling 'cmatrixtrinverse': looks like one of arguments has wrong size")
        n = safe_cols("'cmatrixtrinverse': incorrect parameters",a)
        isunit = False
    else:
        raise RuntimeError("Error while calling 'cmatrixtrinverse': function must have 2 or 4 parameters")
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isunit = ctypes.c_uint8(isunit)
    if __isunit.value!=0:
        __isunit = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_matinvreport()
    x_matinvreport_zero_fields(__rep)
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixtrinverse(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isunit), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixtrinverse'")
        __r__a = listlist_from_x(__a)
        __r__info = __info.value
        __r__rep = matinvreport_from_x(__rep)
        return (__r__a, __r__info, __r__rep)
    finally:
        x_matrix_clear(__a)
        x_matinvreport_clear(__rep)


_lib_alglib.alglib_fisherlda.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fisherlda.restype = ctypes.c_int32
def fisherlda(xy, npoints, nvars, nclasses):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __nclasses = c_ptrint_t(nclasses)
    if __nclasses.value!=nclasses:
        raise ValueError("Error while converting 'nclasses' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fisherlda(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__nclasses), ctypes.byref(__info), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fisherlda'")
        __r__info = __info.value
        __r__w = list_from_x(__w)
        return (__r__info, __r__w)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__w)


_lib_alglib.alglib_fisherldan.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fisherldan.restype = ctypes.c_int32
def fisherldan(xy, npoints, nvars, nclasses):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __nclasses = c_ptrint_t(nclasses)
    if __nclasses.value!=nclasses:
        raise ValueError("Error while converting 'nclasses' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __w = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fisherldan(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__nclasses), ctypes.byref(__info), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fisherldan'")
        __r__info = __info.value
        __r__w = listlist_from_x(__w)
        return (__r__info, __r__w)
    finally:
        x_matrix_clear(__xy)
        x_matrix_clear(__w)


_lib_alglib.alglib_gammafunction.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gammafunction.restype = ctypes.c_int32
def gammafunction(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gammafunction(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gammafunction'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_lngamma.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lngamma.restype = ctypes.c_int32
def lngamma(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __sgngam = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lngamma(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__sgngam))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lngamma'")
        __r__result = __result.value
        __r__sgngam = __sgngam.value
        return (__r__result, __r__sgngam)
    finally:
        pass


_lib_alglib.alglib_errorfunction.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_errorfunction.restype = ctypes.c_int32
def errorfunction(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_errorfunction(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'errorfunction'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_errorfunctionc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_errorfunctionc.restype = ctypes.c_int32
def errorfunctionc(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_errorfunctionc(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'errorfunctionc'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_normaldistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_normaldistribution.restype = ctypes.c_int32
def normaldistribution(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_normaldistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'normaldistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_inverf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_inverf.restype = ctypes.c_int32
def inverf(e):
    pass
    __result = ctypes.c_double(0)
    __e = ctypes.c_double(e)
    if __e.value!=e:
        raise ValueError("Error while converting 'e' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_inverf(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__e))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'inverf'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invnormaldistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invnormaldistribution.restype = ctypes.c_int32
def invnormaldistribution(y0):
    pass
    __result = ctypes.c_double(0)
    __y0 = ctypes.c_double(y0)
    if __y0.value!=y0:
        raise ValueError("Error while converting 'y0' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invnormaldistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__y0))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invnormaldistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_incompletegamma.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_incompletegamma.restype = ctypes.c_int32
def incompletegamma(a, x):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_incompletegamma(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'incompletegamma'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_incompletegammac.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_incompletegammac.restype = ctypes.c_int32
def incompletegammac(a, x):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_incompletegammac(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'incompletegammac'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invincompletegammac.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invincompletegammac.restype = ctypes.c_int32
def invincompletegammac(a, y0):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __y0 = ctypes.c_double(y0)
    if __y0.value!=y0:
        raise ValueError("Error while converting 'y0' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invincompletegammac(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__y0))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invincompletegammac'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_rmatrixbdsvd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixbdsvd.restype = ctypes.c_int32
def rmatrixbdsvd(d, e, n, isupper, isfractionalaccuracyrequired, u, nru, c, ncc, vt, ncvt):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_vector(d):
        raise ValueError("'d' parameter can't be cast to real_vector")
    __d = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(e):
        raise ValueError("'e' parameter can't be cast to real_vector")
    __e = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    __isfractionalaccuracyrequired = ctypes.c_uint8(isfractionalaccuracyrequired)
    if __isfractionalaccuracyrequired.value!=0:
        __isfractionalaccuracyrequired = ctypes.c_uint8(1)
    if not is_real_matrix(u):
        raise ValueError("'u' parameter can't be cast to real_matrix")
    __u = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nru = c_ptrint_t(nru)
    if __nru.value!=nru:
        raise ValueError("Error while converting 'nru' parameter to 'c_ptrint_t'")
    if not is_real_matrix(c):
        raise ValueError("'c' parameter can't be cast to real_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ncc = c_ptrint_t(ncc)
    if __ncc.value!=ncc:
        raise ValueError("Error while converting 'ncc' parameter to 'c_ptrint_t'")
    if not is_real_matrix(vt):
        raise ValueError("'vt' parameter can't be cast to real_matrix")
    __vt = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ncvt = c_ptrint_t(ncvt)
    if __ncvt.value!=ncvt:
        raise ValueError("Error while converting 'ncvt' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__d, d, DT_REAL, X_CREATE)
        x_from_list(__e, e, DT_REAL, X_CREATE)
        x_from_listlist(__u, u, DT_REAL, X_CREATE)
        x_from_listlist(__c, c, DT_REAL, X_CREATE)
        x_from_listlist(__vt, vt, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixbdsvd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__d), ctypes.byref(__e), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__isfractionalaccuracyrequired), ctypes.byref(__u), ctypes.byref(__nru), ctypes.byref(__c), ctypes.byref(__ncc), ctypes.byref(__vt), ctypes.byref(__ncvt))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixbdsvd'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__u = listlist_from_x(__u)
        __r__c = listlist_from_x(__c)
        __r__vt = listlist_from_x(__vt)
        return (__r__result, __r__d, __r__u, __r__c, __r__vt)
    finally:
        x_vector_clear(__d)
        x_vector_clear(__e)
        x_matrix_clear(__u)
        x_matrix_clear(__c)
        x_matrix_clear(__vt)


_lib_alglib.alglib_rmatrixsvd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixsvd.restype = ctypes.c_int32
def rmatrixsvd(a, m, n, uneeded, vtneeded, additionalmemory):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __uneeded = c_ptrint_t(uneeded)
    if __uneeded.value!=uneeded:
        raise ValueError("Error while converting 'uneeded' parameter to 'c_ptrint_t'")
    __vtneeded = c_ptrint_t(vtneeded)
    if __vtneeded.value!=vtneeded:
        raise ValueError("Error while converting 'vtneeded' parameter to 'c_ptrint_t'")
    __additionalmemory = c_ptrint_t(additionalmemory)
    if __additionalmemory.value!=additionalmemory:
        raise ValueError("Error while converting 'additionalmemory' parameter to 'c_ptrint_t'")
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __u = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __vt = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixsvd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__uneeded), ctypes.byref(__vtneeded), ctypes.byref(__additionalmemory), ctypes.byref(__w), ctypes.byref(__u), ctypes.byref(__vt))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixsvd'")
        __r__result = __result.value!=0
        __r__w = list_from_x(__w)
        __r__u = listlist_from_x(__u)
        __r__vt = listlist_from_x(__vt)
        return (__r__result, __r__w, __r__u, __r__vt)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__w)
        x_matrix_clear(__u)
        x_matrix_clear(__vt)


_lib_alglib.x_obj_free_linearmodel.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_linearmodel.restype = None


class linearmodel(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_linearmodel(self.ptr)


class x_lrreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("c", x_matrix),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("cvrmserror", ctypes.c_double),
        ("cvavgerror", ctypes.c_double),
        ("cvavgrelerror", ctypes.c_double),
        ("ncvdefects", c_ptrint_t),
        ("cvdefects", x_vector)
        ]




class lrreport(object):
    def __init__(self):
        self.c = [[]]
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.cvrmserror = 0
        self.cvavgerror = 0
        self.cvavgrelerror = 0
        self.ncvdefects = 0
        self.cvdefects = []


def x_lrreport_zero_fields(x):
    x.c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.cvrmserror = 0
    x.cvavgerror = 0
    x.cvavgrelerror = 0
    x.ncvdefects = 0
    x.cvdefects = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    return




def x_lrreport_clear(x):
    x_matrix_clear(x.c)
    x_vector_clear(x.cvdefects)
    x_lrreport_zero_fields(x)
    return




def x_from_lrreport(x,v):
    x_from_listlist(x.c, v.c, DT_REAL, X_CREATE)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.cvrmserror = float(v.cvrmserror)
    x.cvavgerror = float(v.cvavgerror)
    x.cvavgrelerror = float(v.cvavgrelerror)
    x.ncvdefects = int(v.ncvdefects)
    x_from_list(x.cvdefects, v.cvdefects, DT_INT, X_CREATE)
    return




def lrreport_from_x(x):
    r = lrreport()
    r.c = listlist_from_x(x.c)
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.cvrmserror = x.cvrmserror
    r.cvavgerror = x.cvavgerror
    r.cvavgrelerror = x.cvavgrelerror
    r.ncvdefects = x.ncvdefects
    r.cvdefects = list_from_x(x.cvdefects)
    return r


_lib_alglib.alglib_lrbuild.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrbuild.restype = ctypes.c_int32
def lrbuild(xy, npoints, nvars):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __lm = ctypes.c_void_p(0)
    __ar = x_lrreport()
    x_lrreport_zero_fields(__ar)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrbuild(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__info), ctypes.byref(__lm), ctypes.byref(__ar))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrbuild'")
        __r__info = __info.value
        __r__lm = linearmodel(__lm)
        __r__ar = lrreport_from_x(__ar)
        return (__r__info, __r__lm, __r__ar)
    finally:
        x_matrix_clear(__xy)
        x_lrreport_clear(__ar)


_lib_alglib.alglib_lrbuilds.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrbuilds.restype = ctypes.c_int32
def lrbuilds(xy, s, npoints, nvars):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(s):
        raise ValueError("'s' parameter can't be cast to real_vector")
    __s = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __lm = ctypes.c_void_p(0)
    __ar = x_lrreport()
    x_lrreport_zero_fields(__ar)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__s, s, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrbuilds(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__s), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__info), ctypes.byref(__lm), ctypes.byref(__ar))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrbuilds'")
        __r__info = __info.value
        __r__lm = linearmodel(__lm)
        __r__ar = lrreport_from_x(__ar)
        return (__r__info, __r__lm, __r__ar)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__s)
        x_lrreport_clear(__ar)


_lib_alglib.alglib_lrbuildzs.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrbuildzs.restype = ctypes.c_int32
def lrbuildzs(xy, s, npoints, nvars):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(s):
        raise ValueError("'s' parameter can't be cast to real_vector")
    __s = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __lm = ctypes.c_void_p(0)
    __ar = x_lrreport()
    x_lrreport_zero_fields(__ar)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__s, s, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrbuildzs(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__s), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__info), ctypes.byref(__lm), ctypes.byref(__ar))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrbuildzs'")
        __r__info = __info.value
        __r__lm = linearmodel(__lm)
        __r__ar = lrreport_from_x(__ar)
        return (__r__info, __r__lm, __r__ar)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__s)
        x_lrreport_clear(__ar)


_lib_alglib.alglib_lrbuildz.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrbuildz.restype = ctypes.c_int32
def lrbuildz(xy, npoints, nvars):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __lm = ctypes.c_void_p(0)
    __ar = x_lrreport()
    x_lrreport_zero_fields(__ar)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrbuildz(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__info), ctypes.byref(__lm), ctypes.byref(__ar))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrbuildz'")
        __r__info = __info.value
        __r__lm = linearmodel(__lm)
        __r__ar = lrreport_from_x(__ar)
        return (__r__info, __r__lm, __r__ar)
    finally:
        x_matrix_clear(__xy)
        x_lrreport_clear(__ar)


_lib_alglib.alglib_lrunpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrunpack.restype = ctypes.c_int32
def lrunpack(lm):
    pass
    __lm = lm.ptr
    __v = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nvars = c_ptrint_t(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrunpack(ctypes.byref(_error_msg), ctypes.byref(__lm), ctypes.byref(__v), ctypes.byref(__nvars))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrunpack'")
        __r__v = list_from_x(__v)
        __r__nvars = __nvars.value
        return (__r__v, __r__nvars)
    finally:
        x_vector_clear(__v)


_lib_alglib.alglib_lrpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrpack.restype = ctypes.c_int32
def lrpack(v, nvars):
    pass
    if not is_real_vector(v):
        raise ValueError("'v' parameter can't be cast to real_vector")
    __v = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __lm = ctypes.c_void_p(0)
    try:
        x_from_list(__v, v, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrpack(ctypes.byref(_error_msg), ctypes.byref(__v), ctypes.byref(__nvars), ctypes.byref(__lm))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrpack'")
        __r__lm = linearmodel(__lm)
        return __r__lm
    finally:
        x_vector_clear(__v)


_lib_alglib.alglib_lrprocess.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrprocess.restype = ctypes.c_int32
def lrprocess(lm, x):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrprocess(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrprocess'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_lrrmserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lrrmserror.restype = ctypes.c_int32
def lrrmserror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lrrmserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lrrmserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_lravgerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lravgerror.restype = ctypes.c_int32
def lravgerror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lravgerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lravgerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_lravgrelerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lravgrelerror.restype = ctypes.c_int32
def lravgrelerror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lravgrelerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lravgrelerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.x_obj_free_multilayerperceptron.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_multilayerperceptron.restype = None


class multilayerperceptron(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_multilayerperceptron(self.ptr)
_lib_alglib.alglib_mlpcreate0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreate0.restype = ctypes.c_int32
def mlpcreate0(nin, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreate0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreate0'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreate1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreate1.restype = ctypes.c_int32
def mlpcreate1(nin, nhid, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreate1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreate1'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreate2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreate2.restype = ctypes.c_int32
def mlpcreate2(nin, nhid1, nhid2, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreate2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreate2'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreateb0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreateb0.restype = ctypes.c_int32
def mlpcreateb0(nin, nout, b, d):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreateb0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreateb0'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreateb1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreateb1.restype = ctypes.c_int32
def mlpcreateb1(nin, nhid, nout, b, d):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreateb1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreateb1'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreateb2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreateb2.restype = ctypes.c_int32
def mlpcreateb2(nin, nhid1, nhid2, nout, b, d):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreateb2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreateb2'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreater0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreater0.restype = ctypes.c_int32
def mlpcreater0(nin, nout, a, b):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreater0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreater0'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreater1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreater1.restype = ctypes.c_int32
def mlpcreater1(nin, nhid, nout, a, b):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreater1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreater1'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreater2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreater2.restype = ctypes.c_int32
def mlpcreater2(nin, nhid1, nhid2, nout, a, b):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreater2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreater2'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreatec0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreatec0.restype = ctypes.c_int32
def mlpcreatec0(nin, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreatec0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreatec0'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreatec1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreatec1.restype = ctypes.c_int32
def mlpcreatec1(nin, nhid, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreatec1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreatec1'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlpcreatec2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpcreatec2.restype = ctypes.c_int32
def mlpcreatec2(nin, nhid1, nhid2, nout):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __network = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpcreatec2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpcreatec2'")
        __r__network = multilayerperceptron(__network)
        return __r__network
    finally:
        pass


_lib_alglib.alglib_mlprandomize.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlprandomize.restype = ctypes.c_int32
def mlprandomize(network):
    pass
    __network = network.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlprandomize(ctypes.byref(_error_msg), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlprandomize'")
        return
    finally:
        pass


_lib_alglib.alglib_mlprandomizefull.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlprandomizefull.restype = ctypes.c_int32
def mlprandomizefull(network):
    pass
    __network = network.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlprandomizefull(ctypes.byref(_error_msg), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlprandomizefull'")
        return
    finally:
        pass


_lib_alglib.alglib_mlpproperties.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpproperties.restype = ctypes.c_int32
def mlpproperties(network):
    pass
    __network = network.ptr
    __nin = c_ptrint_t(0)
    __nout = c_ptrint_t(0)
    __wcount = c_ptrint_t(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpproperties(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__wcount))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpproperties'")
        __r__nin = __nin.value
        __r__nout = __nout.value
        __r__wcount = __wcount.value
        return (__r__nin, __r__nout, __r__wcount)
    finally:
        pass


_lib_alglib.alglib_mlpissoftmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpissoftmax.restype = ctypes.c_int32
def mlpissoftmax(network):
    pass
    __result = ctypes.c_uint8(0)
    __network = network.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpissoftmax(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpissoftmax'")
        __r__result = __result.value!=0
        return __r__result
    finally:
        pass


_lib_alglib.alglib_mlpprocess.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpprocess.restype = ctypes.c_int32
def mlpprocess(network, x, y):
    pass
    __network = network.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpprocess(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpprocess'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mlpprocessi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpprocessi.restype = ctypes.c_int32
def mlpprocessi(network, x):
    pass
    __network = network.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __y = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpprocessi(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpprocessi'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mlperror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlperror.restype = ctypes.c_int32
def mlperror(network, xy, ssize):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlperror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlperror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlperrorn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlperrorn.restype = ctypes.c_int32
def mlperrorn(network, xy, ssize):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlperrorn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlperrorn'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpclserror.restype = ctypes.c_int32
def mlpclserror(network, xy, ssize):
    pass
    __result = c_ptrint_t(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlprelclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlprelclserror.restype = ctypes.c_int32
def mlprelclserror(network, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlprelclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlprelclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpavgce.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpavgce.restype = ctypes.c_int32
def mlpavgce(network, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpavgce(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpavgce'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlprmserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlprmserror.restype = ctypes.c_int32
def mlprmserror(network, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlprmserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlprmserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpavgerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpavgerror.restype = ctypes.c_int32
def mlpavgerror(network, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpavgerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpavgerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpavgrelerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpavgrelerror.restype = ctypes.c_int32
def mlpavgrelerror(network, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpavgrelerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpavgrelerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpgrad.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpgrad.restype = ctypes.c_int32
def mlpgrad(network, x, desiredy, grad):
    pass
    __network = network.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(desiredy):
        raise ValueError("'desiredy' parameter can't be cast to real_vector")
    __desiredy = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__desiredy, desiredy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpgrad(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__x), ctypes.byref(__desiredy), ctypes.byref(__e), ctypes.byref(__grad))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpgrad'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        return (__r__e, __r__grad)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__desiredy)
        x_vector_clear(__grad)


_lib_alglib.alglib_mlpgradn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpgradn.restype = ctypes.c_int32
def mlpgradn(network, x, desiredy, grad):
    pass
    __network = network.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(desiredy):
        raise ValueError("'desiredy' parameter can't be cast to real_vector")
    __desiredy = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__desiredy, desiredy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpgradn(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__x), ctypes.byref(__desiredy), ctypes.byref(__e), ctypes.byref(__grad))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpgradn'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        return (__r__e, __r__grad)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__desiredy)
        x_vector_clear(__grad)


_lib_alglib.alglib_mlpgradbatch.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpgradbatch.restype = ctypes.c_int32
def mlpgradbatch(network, xy, ssize, grad):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpgradbatch(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize), ctypes.byref(__e), ctypes.byref(__grad))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpgradbatch'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        return (__r__e, __r__grad)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__grad)


_lib_alglib.alglib_mlpgradnbatch.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpgradnbatch.restype = ctypes.c_int32
def mlpgradnbatch(network, xy, ssize, grad):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpgradnbatch(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize), ctypes.byref(__e), ctypes.byref(__grad))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpgradnbatch'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        return (__r__e, __r__grad)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__grad)


_lib_alglib.alglib_mlphessiannbatch.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlphessiannbatch.restype = ctypes.c_int32
def mlphessiannbatch(network, xy, ssize, grad, h):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(h):
        raise ValueError("'h' parameter can't be cast to real_matrix")
    __h = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        x_from_listlist(__h, h, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlphessiannbatch(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize), ctypes.byref(__e), ctypes.byref(__grad), ctypes.byref(__h))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlphessiannbatch'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        __r__h = listlist_from_x(__h)
        return (__r__e, __r__grad, __r__h)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__grad)
        x_matrix_clear(__h)


_lib_alglib.alglib_mlphessianbatch.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlphessianbatch.restype = ctypes.c_int32
def mlphessianbatch(network, xy, ssize, grad, h):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    __e = ctypes.c_double(0)
    if not is_real_vector(grad):
        raise ValueError("'grad' parameter can't be cast to real_vector")
    __grad = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(h):
        raise ValueError("'h' parameter can't be cast to real_matrix")
    __h = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        x_from_list(__grad, grad, DT_REAL, X_CREATE)
        x_from_listlist(__h, h, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlphessianbatch(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__ssize), ctypes.byref(__e), ctypes.byref(__grad), ctypes.byref(__h))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlphessianbatch'")
        __r__e = __e.value
        __r__grad = list_from_x(__grad)
        __r__h = listlist_from_x(__h)
        return (__r__e, __r__grad, __r__h)
    finally:
        x_matrix_clear(__xy)
        x_vector_clear(__grad)
        x_matrix_clear(__h)




class x_densesolverreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("r1", ctypes.c_double),
        ("rinf", ctypes.c_double)
        ]




class densesolverreport(object):
    def __init__(self):
        self.r1 = 0
        self.rinf = 0


def x_densesolverreport_zero_fields(x):
    x.r1 = 0
    x.rinf = 0
    return




def x_densesolverreport_clear(x):
    x_densesolverreport_zero_fields(x)
    return




def x_from_densesolverreport(x,v):
    x.r1 = float(v.r1)
    x.rinf = float(v.rinf)
    return




def densesolverreport_from_x(x):
    r = densesolverreport()
    r.r1 = x.r1
    r.rinf = x.rinf
    return r




class x_densesolverlsreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("r2", ctypes.c_double),
        ("cx", x_matrix),
        ("n", c_ptrint_t),
        ("k", c_ptrint_t)
        ]




class densesolverlsreport(object):
    def __init__(self):
        self.r2 = 0
        self.cx = [[]]
        self.n = 0
        self.k = 0


def x_densesolverlsreport_zero_fields(x):
    x.r2 = 0
    x.cx = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    x.n = 0
    x.k = 0
    return




def x_densesolverlsreport_clear(x):
    x_matrix_clear(x.cx)
    x_densesolverlsreport_zero_fields(x)
    return




def x_from_densesolverlsreport(x,v):
    x.r2 = float(v.r2)
    x_from_listlist(x.cx, v.cx, DT_REAL, X_CREATE)
    x.n = int(v.n)
    x.k = int(v.k)
    return




def densesolverlsreport_from_x(x):
    r = densesolverlsreport()
    r.r2 = x.r2
    r.cx = listlist_from_x(x.cx)
    r.n = x.n
    r.k = x.k
    return r


_lib_alglib.alglib_rmatrixsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixsolve.restype = ctypes.c_int32
def rmatrixsolve(a, n, b):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_rmatrixsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixsolvem.restype = ctypes.c_int32
def rmatrixsolvem(a, n, b, m, rfs):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __rfs = ctypes.c_uint8(rfs)
    if __rfs.value!=0:
        __rfs = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__rfs), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_rmatrixlusolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlusolve.restype = ctypes.c_int32
def rmatrixlusolve(lua, p, n, b):
    pass
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlusolve(ctypes.byref(_error_msg), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlusolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_rmatrixlusolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixlusolvem.restype = ctypes.c_int32
def rmatrixlusolvem(lua, p, n, b, m):
    pass
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixlusolvem(ctypes.byref(_error_msg), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixlusolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_rmatrixmixedsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixmixedsolve.restype = ctypes.c_int32
def rmatrixmixedsolve(a, lua, p, n, b):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixmixedsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixmixedsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_rmatrixmixedsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixmixedsolvem.restype = ctypes.c_int32
def rmatrixmixedsolvem(a, lua, p, n, b, m):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to real_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__lua, lua, DT_REAL, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixmixedsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixmixedsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixsolvem.restype = ctypes.c_int32
def cmatrixsolvem(a, n, b, m, rfs):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __rfs = ctypes.c_uint8(rfs)
    if __rfs.value!=0:
        __rfs = ctypes.c_uint8(1)
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__rfs), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixsolve.restype = ctypes.c_int32
def cmatrixsolve(a, n, b):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_cmatrixlusolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlusolvem.restype = ctypes.c_int32
def cmatrixlusolvem(lua, p, n, b, m):
    pass
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlusolvem(ctypes.byref(_error_msg), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlusolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixlusolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixlusolve.restype = ctypes.c_int32
def cmatrixlusolve(lua, p, n, b):
    pass
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixlusolve(ctypes.byref(_error_msg), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixlusolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_cmatrixmixedsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixmixedsolvem.restype = ctypes.c_int32
def cmatrixmixedsolvem(a, lua, p, n, b, m):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixmixedsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixmixedsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_cmatrixmixedsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixmixedsolve.restype = ctypes.c_int32
def cmatrixmixedsolve(a, lua, p, n, b):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_complex_matrix(lua):
        raise ValueError("'lua' parameter can't be cast to complex_matrix")
    __lua = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(p):
        raise ValueError("'p' parameter can't be cast to int_vector")
    __p = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__lua, lua, DT_COMPLEX, X_CREATE)
        x_from_list(__p, p, DT_INT, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixmixedsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__lua), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixmixedsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__lua)
        x_vector_clear(__p)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_spdmatrixsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixsolvem.restype = ctypes.c_int32
def spdmatrixsolvem(a, n, isupper, b, m):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_spdmatrixsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixsolve.restype = ctypes.c_int32
def spdmatrixsolve(a, n, isupper, b):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_spdmatrixcholeskysolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholeskysolvem.restype = ctypes.c_int32
def spdmatrixcholeskysolvem(cha, n, isupper, b, m):
    pass
    if not is_real_matrix(cha):
        raise ValueError("'cha' parameter can't be cast to real_matrix")
    __cha = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__cha, cha, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholeskysolvem(ctypes.byref(_error_msg), ctypes.byref(__cha), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholeskysolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__cha)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_spdmatrixcholeskysolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholeskysolve.restype = ctypes.c_int32
def spdmatrixcholeskysolve(cha, n, isupper, b):
    pass
    if not is_real_matrix(cha):
        raise ValueError("'cha' parameter can't be cast to real_matrix")
    __cha = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__cha, cha, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholeskysolve(ctypes.byref(_error_msg), ctypes.byref(__cha), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholeskysolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__cha)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_hpdmatrixsolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixsolvem.restype = ctypes.c_int32
def hpdmatrixsolvem(a, n, isupper, b, m):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixsolvem(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixsolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_hpdmatrixsolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixsolve.restype = ctypes.c_int32
def hpdmatrixsolve(a, n, isupper, b):
    pass
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixsolve(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixsolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_hpdmatrixcholeskysolvem.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixcholeskysolvem.restype = ctypes.c_int32
def hpdmatrixcholeskysolvem(cha, n, isupper, b, m):
    pass
    if not is_complex_matrix(cha):
        raise ValueError("'cha' parameter can't be cast to complex_matrix")
    __cha = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_complex_matrix(b):
        raise ValueError("'b' parameter can't be cast to complex_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__cha, cha, DT_COMPLEX, X_CREATE)
        x_from_listlist(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixcholeskysolvem(ctypes.byref(_error_msg), ctypes.byref(__cha), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixcholeskysolvem'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = listlist_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__cha)
        x_matrix_clear(__b)
        x_densesolverreport_clear(__rep)
        x_matrix_clear(__x)


_lib_alglib.alglib_hpdmatrixcholeskysolve.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hpdmatrixcholeskysolve.restype = ctypes.c_int32
def hpdmatrixcholeskysolve(cha, n, isupper, b):
    pass
    if not is_complex_matrix(cha):
        raise ValueError("'cha' parameter can't be cast to complex_matrix")
    __cha = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __info = c_ptrint_t(0)
    __rep = x_densesolverreport()
    x_densesolverreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__cha, cha, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hpdmatrixcholeskysolve(ctypes.byref(_error_msg), ctypes.byref(__cha), ctypes.byref(__n), ctypes.byref(__isupper), ctypes.byref(__b), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hpdmatrixcholeskysolve'")
        __r__info = __info.value
        __r__rep = densesolverreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__cha)
        x_vector_clear(__b)
        x_densesolverreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.alglib_rmatrixsolvels.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixsolvels.restype = ctypes.c_int32
def rmatrixsolvels(a, nrows, ncols, b, threshold):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nrows = c_ptrint_t(nrows)
    if __nrows.value!=nrows:
        raise ValueError("Error while converting 'nrows' parameter to 'c_ptrint_t'")
    __ncols = c_ptrint_t(ncols)
    if __ncols.value!=ncols:
        raise ValueError("Error while converting 'ncols' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __threshold = ctypes.c_double(threshold)
    if __threshold.value!=threshold:
        raise ValueError("Error while converting 'threshold' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __rep = x_densesolverlsreport()
    x_densesolverlsreport_zero_fields(__rep)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixsolvels(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__nrows), ctypes.byref(__ncols), ctypes.byref(__b), ctypes.byref(__threshold), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixsolvels'")
        __r__info = __info.value
        __r__rep = densesolverlsreport_from_x(__rep)
        __r__x = list_from_x(__x)
        return (__r__info, __r__rep, __r__x)
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__b)
        x_densesolverlsreport_clear(__rep)
        x_vector_clear(__x)


_lib_alglib.x_obj_free_logitmodel.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_logitmodel.restype = None


class logitmodel(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_logitmodel(self.ptr)


class x_mnlreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("ngrad", c_ptrint_t),
        ("nhess", c_ptrint_t)
        ]




class mnlreport(object):
    def __init__(self):
        self.ngrad = 0
        self.nhess = 0


def x_mnlreport_zero_fields(x):
    x.ngrad = 0
    x.nhess = 0
    return




def x_mnlreport_clear(x):
    x_mnlreport_zero_fields(x)
    return




def x_from_mnlreport(x,v):
    x.ngrad = int(v.ngrad)
    x.nhess = int(v.nhess)
    return




def mnlreport_from_x(x):
    r = mnlreport()
    r.ngrad = x.ngrad
    r.nhess = x.nhess
    return r


_lib_alglib.alglib_mnltrainh.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnltrainh.restype = ctypes.c_int32
def mnltrainh(xy, npoints, nvars, nclasses):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __nclasses = c_ptrint_t(nclasses)
    if __nclasses.value!=nclasses:
        raise ValueError("Error while converting 'nclasses' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __lm = ctypes.c_void_p(0)
    __rep = x_mnlreport()
    x_mnlreport_zero_fields(__rep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnltrainh(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__nclasses), ctypes.byref(__info), ctypes.byref(__lm), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnltrainh'")
        __r__info = __info.value
        __r__lm = logitmodel(__lm)
        __r__rep = mnlreport_from_x(__rep)
        return (__r__info, __r__lm, __r__rep)
    finally:
        x_matrix_clear(__xy)
        x_mnlreport_clear(__rep)


_lib_alglib.alglib_mnlprocess.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlprocess.restype = ctypes.c_int32
def mnlprocess(lm, x, y):
    pass
    __lm = lm.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlprocess(ctypes.byref(_error_msg), ctypes.byref(__lm), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlprocess'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mnlprocessi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlprocessi.restype = ctypes.c_int32
def mnlprocessi(lm, x):
    pass
    __lm = lm.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __y = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlprocessi(ctypes.byref(_error_msg), ctypes.byref(__lm), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlprocessi'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mnlunpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlunpack.restype = ctypes.c_int32
def mnlunpack(lm):
    pass
    __lm = lm.ptr
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nvars = c_ptrint_t(0)
    __nclasses = c_ptrint_t(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlunpack(ctypes.byref(_error_msg), ctypes.byref(__lm), ctypes.byref(__a), ctypes.byref(__nvars), ctypes.byref(__nclasses))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlunpack'")
        __r__a = listlist_from_x(__a)
        __r__nvars = __nvars.value
        __r__nclasses = __nclasses.value
        return (__r__a, __r__nvars, __r__nclasses)
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_mnlpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlpack.restype = ctypes.c_int32
def mnlpack(a, nvars, nclasses):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __nclasses = c_ptrint_t(nclasses)
    if __nclasses.value!=nclasses:
        raise ValueError("Error while converting 'nclasses' parameter to 'c_ptrint_t'")
    __lm = ctypes.c_void_p(0)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlpack(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__nvars), ctypes.byref(__nclasses), ctypes.byref(__lm))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlpack'")
        __r__lm = logitmodel(__lm)
        return __r__lm
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_mnlavgce.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlavgce.restype = ctypes.c_int32
def mnlavgce(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlavgce(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlavgce'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mnlrelclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlrelclserror.restype = ctypes.c_int32
def mnlrelclserror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlrelclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlrelclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mnlrmserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlrmserror.restype = ctypes.c_int32
def mnlrmserror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlrmserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlrmserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mnlavgerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlavgerror.restype = ctypes.c_int32
def mnlavgerror(lm, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlavgerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlavgerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mnlavgrelerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlavgrelerror.restype = ctypes.c_int32
def mnlavgrelerror(lm, xy, ssize):
    pass
    __result = ctypes.c_double(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ssize = c_ptrint_t(ssize)
    if __ssize.value!=ssize:
        raise ValueError("Error while converting 'ssize' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlavgrelerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__ssize))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlavgrelerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mnlclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mnlclserror.restype = ctypes.c_int32
def mnlclserror(lm, xy, npoints):
    pass
    __result = c_ptrint_t(0)
    __lm = lm.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mnlclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__lm), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mnlclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.x_obj_free_minlbfgsstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_minlbfgsstate.restype = None
_lib_alglib.x_minlbfgsstate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_get_needfg.restype = None
_lib_alglib.x_minlbfgsstate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_set_needfg.restype = None
_lib_alglib.x_minlbfgsstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_get_xupdated.restype = None
_lib_alglib.x_minlbfgsstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_set_xupdated.restype = None
_lib_alglib.x_minlbfgsstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_get_f.restype = None
_lib_alglib.x_minlbfgsstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_set_f.restype = None
_lib_alglib.x_minlbfgsstate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_get_g.restype = None
_lib_alglib.x_minlbfgsstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlbfgsstate_get_x.restype = None


class minlbfgsstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_minlbfgsstate(self.ptr)


class x_minlbfgsreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("iterationscount", c_ptrint_t),
        ("nfev", c_ptrint_t),
        ("terminationtype", c_ptrint_t)
        ]




class minlbfgsreport(object):
    def __init__(self):
        self.iterationscount = 0
        self.nfev = 0
        self.terminationtype = 0


def x_minlbfgsreport_zero_fields(x):
    x.iterationscount = 0
    x.nfev = 0
    x.terminationtype = 0
    return




def x_minlbfgsreport_clear(x):
    x_minlbfgsreport_zero_fields(x)
    return




def x_from_minlbfgsreport(x,v):
    x.iterationscount = int(v.iterationscount)
    x.nfev = int(v.nfev)
    x.terminationtype = int(v.terminationtype)
    return




def minlbfgsreport_from_x(x):
    r = minlbfgsreport()
    r.iterationscount = x.iterationscount
    r.nfev = x.nfev
    r.terminationtype = x.terminationtype
    return r


_lib_alglib.alglib_minlbfgscreate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgscreate.restype = ctypes.c_int32
def minlbfgscreate(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'minlbfgscreate': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlbfgscreate': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgscreate(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgscreate'")
        __r__state = minlbfgsstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlbfgssetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgssetcond.restype = ctypes.c_int32
def minlbfgssetcond(state, epsg, epsf, epsx, maxits):
    pass
    __state = state.ptr
    __epsg = ctypes.c_double(epsg)
    if __epsg.value!=epsg:
        raise ValueError("Error while converting 'epsg' parameter to 'ctypes.c_double'")
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgssetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsg), ctypes.byref(__epsf), ctypes.byref(__epsx), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgssetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_minlbfgssetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgssetxrep.restype = ctypes.c_int32
def minlbfgssetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgssetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgssetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_minlbfgssetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgssetstpmax.restype = ctypes.c_int32
def minlbfgssetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgssetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgssetstpmax'")
        return
    finally:
        pass


_lib_alglib.alglib_minlbfgssetdefaultpreconditioner.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgssetdefaultpreconditioner.restype = ctypes.c_int32
def minlbfgssetdefaultpreconditioner(state):
    pass
    __state = state.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgssetdefaultpreconditioner(ctypes.byref(_error_msg), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgssetdefaultpreconditioner'")
        return
    finally:
        pass


_lib_alglib.alglib_minlbfgssetcholeskypreconditioner.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgssetcholeskypreconditioner.restype = ctypes.c_int32
def minlbfgssetcholeskypreconditioner(state, p, isupper):
    pass
    __state = state.ptr
    if not is_real_matrix(p):
        raise ValueError("'p' parameter can't be cast to real_matrix")
    __p = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__p, p, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgssetcholeskypreconditioner(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__p), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgssetcholeskypreconditioner'")
        return
    finally:
        x_matrix_clear(__p)




def minlbfgsoptimize_g(state, grad, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlbfgsstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_minlbfgsstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    while True:
        retval = _lib_alglib.alglib_minlbfgsiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgsiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlbfgsstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_minlbfgsstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlbfgsstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlbfgsstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlbfgsoptimize' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_minlbfgsresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgsresults.restype = ctypes.c_int32
def minlbfgsresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minlbfgsreport()
    x_minlbfgsreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgsresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgsresults'")
        __r__x = list_from_x(__x)
        __r__rep = minlbfgsreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minlbfgsreport_clear(__rep)


_lib_alglib.alglib_minlbfgsresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgsresultsbuf.restype = ctypes.c_int32
def minlbfgsresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minlbfgsreport()
    x_minlbfgsreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_minlbfgsreport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgsresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgsresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = minlbfgsreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minlbfgsreport_clear(__rep)


_lib_alglib.alglib_minlbfgsrestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlbfgsrestartfrom.restype = ctypes.c_int32
def minlbfgsrestartfrom(state, x):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlbfgsrestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlbfgsrestartfrom'")
        return
    finally:
        x_vector_clear(__x)




class x_mlpreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("ngrad", c_ptrint_t),
        ("nhess", c_ptrint_t),
        ("ncholesky", c_ptrint_t)
        ]




class mlpreport(object):
    def __init__(self):
        self.ngrad = 0
        self.nhess = 0
        self.ncholesky = 0


def x_mlpreport_zero_fields(x):
    x.ngrad = 0
    x.nhess = 0
    x.ncholesky = 0
    return




def x_mlpreport_clear(x):
    x_mlpreport_zero_fields(x)
    return




def x_from_mlpreport(x,v):
    x.ngrad = int(v.ngrad)
    x.nhess = int(v.nhess)
    x.ncholesky = int(v.ncholesky)
    return




def mlpreport_from_x(x):
    r = mlpreport()
    r.ngrad = x.ngrad
    r.nhess = x.nhess
    r.ncholesky = x.ncholesky
    return r




class x_mlpcvreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("relclserror", ctypes.c_double),
        ("avgce", ctypes.c_double),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double)
        ]




class mlpcvreport(object):
    def __init__(self):
        self.relclserror = 0
        self.avgce = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0


def x_mlpcvreport_zero_fields(x):
    x.relclserror = 0
    x.avgce = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    return




def x_mlpcvreport_clear(x):
    x_mlpcvreport_zero_fields(x)
    return




def x_from_mlpcvreport(x,v):
    x.relclserror = float(v.relclserror)
    x.avgce = float(v.avgce)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    return




def mlpcvreport_from_x(x):
    r = mlpcvreport()
    r.relclserror = x.relclserror
    r.avgce = x.avgce
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    return r


_lib_alglib.alglib_mlptrainlm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlptrainlm.restype = ctypes.c_int32
def mlptrainlm(network, xy, npoints, decay, restarts):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlptrainlm(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlptrainlm'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        return (__r__info, __r__rep)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)


_lib_alglib.alglib_mlptrainlbfgs.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlptrainlbfgs.restype = ctypes.c_int32
def mlptrainlbfgs(network, xy, npoints, decay, restarts, wstep, maxits):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __wstep = ctypes.c_double(wstep)
    if __wstep.value!=wstep:
        raise ValueError("Error while converting 'wstep' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlptrainlbfgs(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__wstep), ctypes.byref(__maxits), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlptrainlbfgs'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        return (__r__info, __r__rep)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)


_lib_alglib.alglib_mlptraines.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlptraines.restype = ctypes.c_int32
def mlptraines(network, trnxy, trnsize, valxy, valsize, decay, restarts):
    pass
    __network = network.ptr
    if not is_real_matrix(trnxy):
        raise ValueError("'trnxy' parameter can't be cast to real_matrix")
    __trnxy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __trnsize = c_ptrint_t(trnsize)
    if __trnsize.value!=trnsize:
        raise ValueError("Error while converting 'trnsize' parameter to 'c_ptrint_t'")
    if not is_real_matrix(valxy):
        raise ValueError("'valxy' parameter can't be cast to real_matrix")
    __valxy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __valsize = c_ptrint_t(valsize)
    if __valsize.value!=valsize:
        raise ValueError("Error while converting 'valsize' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    try:
        x_from_listlist(__trnxy, trnxy, DT_REAL, X_CREATE)
        x_from_listlist(__valxy, valxy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlptraines(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__trnxy), ctypes.byref(__trnsize), ctypes.byref(__valxy), ctypes.byref(__valsize), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlptraines'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        return (__r__info, __r__rep)
    finally:
        x_matrix_clear(__trnxy)
        x_matrix_clear(__valxy)
        x_mlpreport_clear(__rep)


_lib_alglib.alglib_mlpkfoldcvlbfgs.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpkfoldcvlbfgs.restype = ctypes.c_int32
def mlpkfoldcvlbfgs(network, xy, npoints, decay, restarts, wstep, maxits, foldscount):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __wstep = ctypes.c_double(wstep)
    if __wstep.value!=wstep:
        raise ValueError("Error while converting 'wstep' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    __foldscount = c_ptrint_t(foldscount)
    if __foldscount.value!=foldscount:
        raise ValueError("Error while converting 'foldscount' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    __cvrep = x_mlpcvreport()
    x_mlpcvreport_zero_fields(__cvrep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpkfoldcvlbfgs(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__wstep), ctypes.byref(__maxits), ctypes.byref(__foldscount), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__cvrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpkfoldcvlbfgs'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        __r__cvrep = mlpcvreport_from_x(__cvrep)
        return (__r__info, __r__rep, __r__cvrep)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)
        x_mlpcvreport_clear(__cvrep)


_lib_alglib.alglib_mlpkfoldcvlm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpkfoldcvlm.restype = ctypes.c_int32
def mlpkfoldcvlm(network, xy, npoints, decay, restarts, foldscount):
    pass
    __network = network.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __foldscount = c_ptrint_t(foldscount)
    if __foldscount.value!=foldscount:
        raise ValueError("Error while converting 'foldscount' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    __cvrep = x_mlpcvreport()
    x_mlpcvreport_zero_fields(__cvrep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpkfoldcvlm(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__foldscount), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__cvrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpkfoldcvlm'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        __r__cvrep = mlpcvreport_from_x(__cvrep)
        return (__r__info, __r__rep, __r__cvrep)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)
        x_mlpcvreport_clear(__cvrep)


_lib_alglib.x_obj_free_mlpensemble.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_mlpensemble.restype = None


class mlpensemble(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_mlpensemble(self.ptr)
_lib_alglib.alglib_mlpecreate0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreate0.restype = ctypes.c_int32
def mlpecreate0(nin, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreate0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreate0'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreate1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreate1.restype = ctypes.c_int32
def mlpecreate1(nin, nhid, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreate1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreate1'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreate2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreate2.restype = ctypes.c_int32
def mlpecreate2(nin, nhid1, nhid2, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreate2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreate2'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreateb0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreateb0.restype = ctypes.c_int32
def mlpecreateb0(nin, nout, b, d, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreateb0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreateb0'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreateb1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreateb1.restype = ctypes.c_int32
def mlpecreateb1(nin, nhid, nout, b, d, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreateb1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreateb1'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreateb2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreateb2.restype = ctypes.c_int32
def mlpecreateb2(nin, nhid1, nhid2, nout, b, d, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __d = ctypes.c_double(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreateb2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__b), ctypes.byref(__d), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreateb2'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreater0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreater0.restype = ctypes.c_int32
def mlpecreater0(nin, nout, a, b, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreater0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreater0'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreater1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreater1.restype = ctypes.c_int32
def mlpecreater1(nin, nhid, nout, a, b, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreater1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreater1'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreater2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreater2.restype = ctypes.c_int32
def mlpecreater2(nin, nhid1, nhid2, nout, a, b, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreater2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreater2'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreatec0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreatec0.restype = ctypes.c_int32
def mlpecreatec0(nin, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreatec0(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreatec0'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreatec1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreatec1.restype = ctypes.c_int32
def mlpecreatec1(nin, nhid, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid = c_ptrint_t(nhid)
    if __nhid.value!=nhid:
        raise ValueError("Error while converting 'nhid' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreatec1(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreatec1'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreatec2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreatec2.restype = ctypes.c_int32
def mlpecreatec2(nin, nhid1, nhid2, nout, ensemblesize):
    pass
    __nin = c_ptrint_t(nin)
    if __nin.value!=nin:
        raise ValueError("Error while converting 'nin' parameter to 'c_ptrint_t'")
    __nhid1 = c_ptrint_t(nhid1)
    if __nhid1.value!=nhid1:
        raise ValueError("Error while converting 'nhid1' parameter to 'c_ptrint_t'")
    __nhid2 = c_ptrint_t(nhid2)
    if __nhid2.value!=nhid2:
        raise ValueError("Error while converting 'nhid2' parameter to 'c_ptrint_t'")
    __nout = c_ptrint_t(nout)
    if __nout.value!=nout:
        raise ValueError("Error while converting 'nout' parameter to 'c_ptrint_t'")
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreatec2(ctypes.byref(_error_msg), ctypes.byref(__nin), ctypes.byref(__nhid1), ctypes.byref(__nhid2), ctypes.byref(__nout), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreatec2'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlpecreatefromnetwork.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpecreatefromnetwork.restype = ctypes.c_int32
def mlpecreatefromnetwork(network, ensemblesize):
    pass
    __network = network.ptr
    __ensemblesize = c_ptrint_t(ensemblesize)
    if __ensemblesize.value!=ensemblesize:
        raise ValueError("Error while converting 'ensemblesize' parameter to 'c_ptrint_t'")
    __ensemble = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpecreatefromnetwork(ctypes.byref(_error_msg), ctypes.byref(__network), ctypes.byref(__ensemblesize), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpecreatefromnetwork'")
        __r__ensemble = mlpensemble(__ensemble)
        return __r__ensemble
    finally:
        pass


_lib_alglib.alglib_mlperandomize.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlperandomize.restype = ctypes.c_int32
def mlperandomize(ensemble):
    pass
    __ensemble = ensemble.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlperandomize(ctypes.byref(_error_msg), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlperandomize'")
        return
    finally:
        pass


_lib_alglib.alglib_mlpeproperties.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeproperties.restype = ctypes.c_int32
def mlpeproperties(ensemble):
    pass
    __ensemble = ensemble.ptr
    __nin = c_ptrint_t(0)
    __nout = c_ptrint_t(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeproperties(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__nin), ctypes.byref(__nout))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeproperties'")
        __r__nin = __nin.value
        __r__nout = __nout.value
        return (__r__nin, __r__nout)
    finally:
        pass


_lib_alglib.alglib_mlpeissoftmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeissoftmax.restype = ctypes.c_int32
def mlpeissoftmax(ensemble):
    pass
    __result = ctypes.c_uint8(0)
    __ensemble = ensemble.ptr
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeissoftmax(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeissoftmax'")
        __r__result = __result.value!=0
        return __r__result
    finally:
        pass


_lib_alglib.alglib_mlpeprocess.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeprocess.restype = ctypes.c_int32
def mlpeprocess(ensemble, x, y):
    pass
    __ensemble = ensemble.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeprocess(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeprocess'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mlpeprocessi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeprocessi.restype = ctypes.c_int32
def mlpeprocessi(ensemble, x):
    pass
    __ensemble = ensemble.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __y = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeprocessi(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeprocessi'")
        __r__y = list_from_x(__y)
        return __r__y
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_mlperelclserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlperelclserror.restype = ctypes.c_int32
def mlperelclserror(ensemble, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlperelclserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlperelclserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpeavgce.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeavgce.restype = ctypes.c_int32
def mlpeavgce(ensemble, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeavgce(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeavgce'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpermserror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpermserror.restype = ctypes.c_int32
def mlpermserror(ensemble, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpermserror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpermserror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpeavgerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeavgerror.restype = ctypes.c_int32
def mlpeavgerror(ensemble, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeavgerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeavgerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpeavgrelerror.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpeavgrelerror.restype = ctypes.c_int32
def mlpeavgrelerror(ensemble, xy, npoints):
    pass
    __result = ctypes.c_double(0)
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpeavgrelerror(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpeavgrelerror'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_mlpebagginglm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpebagginglm.restype = ctypes.c_int32
def mlpebagginglm(ensemble, xy, npoints, decay, restarts):
    pass
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    __ooberrors = x_mlpcvreport()
    x_mlpcvreport_zero_fields(__ooberrors)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpebagginglm(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__ooberrors))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpebagginglm'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        __r__ooberrors = mlpcvreport_from_x(__ooberrors)
        return (__r__info, __r__rep, __r__ooberrors)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)
        x_mlpcvreport_clear(__ooberrors)


_lib_alglib.alglib_mlpebagginglbfgs.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpebagginglbfgs.restype = ctypes.c_int32
def mlpebagginglbfgs(ensemble, xy, npoints, decay, restarts, wstep, maxits):
    pass
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __wstep = ctypes.c_double(wstep)
    if __wstep.value!=wstep:
        raise ValueError("Error while converting 'wstep' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    __ooberrors = x_mlpcvreport()
    x_mlpcvreport_zero_fields(__ooberrors)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpebagginglbfgs(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__wstep), ctypes.byref(__maxits), ctypes.byref(__info), ctypes.byref(__rep), ctypes.byref(__ooberrors))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpebagginglbfgs'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        __r__ooberrors = mlpcvreport_from_x(__ooberrors)
        return (__r__info, __r__rep, __r__ooberrors)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)
        x_mlpcvreport_clear(__ooberrors)


_lib_alglib.alglib_mlpetraines.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mlpetraines.restype = ctypes.c_int32
def mlpetraines(ensemble, xy, npoints, decay, restarts):
    pass
    __ensemble = ensemble.ptr
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __decay = ctypes.c_double(decay)
    if __decay.value!=decay:
        raise ValueError("Error while converting 'decay' parameter to 'ctypes.c_double'")
    __restarts = c_ptrint_t(restarts)
    if __restarts.value!=restarts:
        raise ValueError("Error while converting 'restarts' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __rep = x_mlpreport()
    x_mlpreport_zero_fields(__rep)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mlpetraines(ctypes.byref(_error_msg), ctypes.byref(__ensemble), ctypes.byref(__xy), ctypes.byref(__npoints), ctypes.byref(__decay), ctypes.byref(__restarts), ctypes.byref(__info), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mlpetraines'")
        __r__info = __info.value
        __r__rep = mlpreport_from_x(__rep)
        return (__r__info, __r__rep)
    finally:
        x_matrix_clear(__xy)
        x_mlpreport_clear(__rep)


_lib_alglib.alglib_pcabuildbasis.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pcabuildbasis.restype = ctypes.c_int32
def pcabuildbasis(x, npoints, nvars):
    pass
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __npoints = c_ptrint_t(npoints)
    if __npoints.value!=npoints:
        raise ValueError("Error while converting 'npoints' parameter to 'c_ptrint_t'")
    __nvars = c_ptrint_t(nvars)
    if __nvars.value!=nvars:
        raise ValueError("Error while converting 'nvars' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __s2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __v = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pcabuildbasis(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__npoints), ctypes.byref(__nvars), ctypes.byref(__info), ctypes.byref(__s2), ctypes.byref(__v))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pcabuildbasis'")
        __r__info = __info.value
        __r__s2 = list_from_x(__s2)
        __r__v = listlist_from_x(__v)
        return (__r__info, __r__s2, __r__v)
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__s2)
        x_matrix_clear(__v)


_lib_alglib.x_obj_free_odesolverstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_odesolverstate.restype = None
_lib_alglib.x_odesolverstate_get_needdy.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_get_needdy.restype = None
_lib_alglib.x_odesolverstate_set_needdy.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_set_needdy.restype = None
_lib_alglib.x_odesolverstate_get_y.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_get_y.restype = None
_lib_alglib.x_odesolverstate_get_dy.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_get_dy.restype = None
_lib_alglib.x_odesolverstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_get_x.restype = None
_lib_alglib.x_odesolverstate_set_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_odesolverstate_set_x.restype = None


class odesolverstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_odesolverstate(self.ptr)


class x_odesolverreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("nfev", c_ptrint_t),
        ("terminationtype", c_ptrint_t)
        ]




class odesolverreport(object):
    def __init__(self):
        self.nfev = 0
        self.terminationtype = 0


def x_odesolverreport_zero_fields(x):
    x.nfev = 0
    x.terminationtype = 0
    return




def x_odesolverreport_clear(x):
    x_odesolverreport_zero_fields(x)
    return




def x_from_odesolverreport(x,v):
    x.nfev = int(v.nfev)
    x.terminationtype = int(v.terminationtype)
    return




def odesolverreport_from_x(x):
    r = odesolverreport()
    r.nfev = x.nfev
    r.terminationtype = x.terminationtype
    return r


_lib_alglib.alglib_odesolverrkck.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_odesolverrkck.restype = ctypes.c_int32
def odesolverrkck(*functionargs):
    if len(functionargs)==6:
        __friendly_form = False
        y,n,x,m,eps,h = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        y,x,eps,h = functionargs
        n = safe_len("'odesolverrkck': incorrect parameters",y)
        m = safe_len("'odesolverrkck': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'odesolverrkck': function must have 4 or 6 parameters")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __eps = ctypes.c_double(eps)
    if __eps.value!=eps:
        raise ValueError("Error while converting 'eps' parameter to 'ctypes.c_double'")
    __h = ctypes.c_double(h)
    if __h.value!=h:
        raise ValueError("Error while converting 'h' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_odesolverrkck(ctypes.byref(_error_msg), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__m), ctypes.byref(__eps), ctypes.byref(__h), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'odesolverrkck'")
        __r__state = odesolverstate(__state)
        return __r__state
    finally:
        x_vector_clear(__y)
        x_vector_clear(__x)




def odesolversolve(state, dy, param = None):
    # initialize temporaries
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_flag = ctypes.c_uint8()
    
    # initialize reverse communication variables
    _xc_y = x_vector()
    _lib_alglib.x_odesolverstate_get_y(state.ptr, ctypes.byref(_xc_y))
    _py_y = create_real_vector(_xc_y.cnt)
    _xc_x = ctypes.c_double()
    _xc_dy = x_vector()
    _lib_alglib.x_odesolverstate_get_dy(state.ptr, ctypes.byref(_xc_dy))
    _py_dy = create_real_vector(_xc_dy.cnt)
    
    # algorithm iterations
    while True:
        retval = _lib_alglib.alglib_odesolveriteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'odesolveriteration'")
        if not _xc_result:
            break
        _lib_alglib.x_odesolverstate_get_needdy(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_y, _py_y)
            _lib_alglib.x_odesolverstate_get_x(state.ptr, ctypes.byref(_xc_x))
            dy(_py_y, _xc_x.value, _py_dy, param)
            x_from_list(_xc_dy, _py_dy, DT_REAL, X_REWRITE)
            continue
        raise RuntimeError("ALGLIB: unexpected error in 'odesolversolve'")
    return


_lib_alglib.alglib_odesolverresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_odesolverresults.restype = ctypes.c_int32
def odesolverresults(state):
    pass
    __state = state.ptr
    __m = c_ptrint_t(0)
    __xtbl = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __ytbl = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_odesolverreport()
    x_odesolverreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_odesolverresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__m), ctypes.byref(__xtbl), ctypes.byref(__ytbl), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'odesolverresults'")
        __r__m = __m.value
        __r__xtbl = list_from_x(__xtbl)
        __r__ytbl = listlist_from_x(__ytbl)
        __r__rep = odesolverreport_from_x(__rep)
        return (__r__m, __r__xtbl, __r__ytbl, __r__rep)
    finally:
        x_vector_clear(__xtbl)
        x_matrix_clear(__ytbl)
        x_odesolverreport_clear(__rep)


_lib_alglib.alglib_fftc1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fftc1d.restype = ctypes.c_int32
def fftc1d(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        n = safe_len("'fftc1d': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'fftc1d': function must have 1 or 2 parameters")
    if not is_complex_vector(a):
        raise ValueError("'a' parameter can't be cast to complex_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fftc1d(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fftc1d'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_fftc1dinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fftc1dinv.restype = ctypes.c_int32
def fftc1dinv(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        n = safe_len("'fftc1dinv': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'fftc1dinv': function must have 1 or 2 parameters")
    if not is_complex_vector(a):
        raise ValueError("'a' parameter can't be cast to complex_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fftc1dinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fftc1dinv'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_fftr1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fftr1d.restype = ctypes.c_int32
def fftr1d(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        n = safe_len("'fftr1d': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'fftr1d': function must have 1 or 2 parameters")
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __f = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fftr1d(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__f))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fftr1d'")
        __r__f = list_from_x(__f)
        return __r__f
    finally:
        x_vector_clear(__a)
        x_vector_clear(__f)


_lib_alglib.alglib_fftr1dinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fftr1dinv.restype = ctypes.c_int32
def fftr1dinv(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        f,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        f, = functionargs
        n = safe_len("'fftr1dinv': incorrect parameters",f)
    else:
        raise RuntimeError("Error while calling 'fftr1dinv': function must have 1 or 2 parameters")
    if not is_complex_vector(f):
        raise ValueError("'f' parameter can't be cast to complex_vector")
    __f = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __a = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__f, f, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fftr1dinv(ctypes.byref(_error_msg), ctypes.byref(__f), ctypes.byref(__n), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fftr1dinv'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__f)
        x_vector_clear(__a)


_lib_alglib.alglib_convc1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convc1d.restype = ctypes.c_int32
def convc1d(a, m, b, n):
    pass
    if not is_complex_vector(a):
        raise ValueError("'a' parameter can't be cast to complex_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convc1d(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convc1d'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_convc1dinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convc1dinv.restype = ctypes.c_int32
def convc1dinv(a, m, b, n):
    pass
    if not is_complex_vector(a):
        raise ValueError("'a' parameter can't be cast to complex_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convc1dinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convc1dinv'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_convc1dcircular.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convc1dcircular.restype = ctypes.c_int32
def convc1dcircular(s, m, r, n):
    pass
    if not is_complex_vector(s):
        raise ValueError("'s' parameter can't be cast to complex_vector")
    __s = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_complex_vector(r):
        raise ValueError("'r' parameter can't be cast to complex_vector")
    __r = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__s, s, DT_COMPLEX, X_CREATE)
        x_from_list(__r, r, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convc1dcircular(ctypes.byref(_error_msg), ctypes.byref(__s), ctypes.byref(__m), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convc1dcircular'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__s)
        x_vector_clear(__r)
        x_vector_clear(__c)


_lib_alglib.alglib_convc1dcircularinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convc1dcircularinv.restype = ctypes.c_int32
def convc1dcircularinv(a, m, b, n):
    pass
    if not is_complex_vector(a):
        raise ValueError("'a' parameter can't be cast to complex_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_complex_vector(b):
        raise ValueError("'b' parameter can't be cast to complex_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__b, b, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convc1dcircularinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convc1dcircularinv'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_convr1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convr1d.restype = ctypes.c_int32
def convr1d(a, m, b, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convr1d(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convr1d'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_convr1dinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convr1dinv.restype = ctypes.c_int32
def convr1dinv(a, m, b, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convr1dinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convr1dinv'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_convr1dcircular.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convr1dcircular.restype = ctypes.c_int32
def convr1dcircular(s, m, r, n):
    pass
    if not is_real_vector(s):
        raise ValueError("'s' parameter can't be cast to real_vector")
    __s = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(r):
        raise ValueError("'r' parameter can't be cast to real_vector")
    __r = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__s, s, DT_REAL, X_CREATE)
        x_from_list(__r, r, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convr1dcircular(ctypes.byref(_error_msg), ctypes.byref(__s), ctypes.byref(__m), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convr1dcircular'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__s)
        x_vector_clear(__r)
        x_vector_clear(__c)


_lib_alglib.alglib_convr1dcircularinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_convr1dcircularinv.restype = ctypes.c_int32
def convr1dcircularinv(a, m, b, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(b):
        raise ValueError("'b' parameter can't be cast to real_vector")
    __b = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        x_from_list(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_convr1dcircularinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__m), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'convr1dcircularinv'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)
        x_vector_clear(__r)


_lib_alglib.alglib_corrc1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_corrc1d.restype = ctypes.c_int32
def corrc1d(signal, n, pattern, m):
    pass
    if not is_complex_vector(signal):
        raise ValueError("'signal' parameter can't be cast to complex_vector")
    __signal = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_complex_vector(pattern):
        raise ValueError("'pattern' parameter can't be cast to complex_vector")
    __pattern = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__signal, signal, DT_COMPLEX, X_CREATE)
        x_from_list(__pattern, pattern, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_corrc1d(ctypes.byref(_error_msg), ctypes.byref(__signal), ctypes.byref(__n), ctypes.byref(__pattern), ctypes.byref(__m), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'corrc1d'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__signal)
        x_vector_clear(__pattern)
        x_vector_clear(__r)


_lib_alglib.alglib_corrc1dcircular.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_corrc1dcircular.restype = ctypes.c_int32
def corrc1dcircular(signal, m, pattern, n):
    pass
    if not is_complex_vector(signal):
        raise ValueError("'signal' parameter can't be cast to complex_vector")
    __signal = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_complex_vector(pattern):
        raise ValueError("'pattern' parameter can't be cast to complex_vector")
    __pattern = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_COMPLEX,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__signal, signal, DT_COMPLEX, X_CREATE)
        x_from_list(__pattern, pattern, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_corrc1dcircular(ctypes.byref(_error_msg), ctypes.byref(__signal), ctypes.byref(__m), ctypes.byref(__pattern), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'corrc1dcircular'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__signal)
        x_vector_clear(__pattern)
        x_vector_clear(__c)


_lib_alglib.alglib_corrr1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_corrr1d.restype = ctypes.c_int32
def corrr1d(signal, n, pattern, m):
    pass
    if not is_real_vector(signal):
        raise ValueError("'signal' parameter can't be cast to real_vector")
    __signal = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(pattern):
        raise ValueError("'pattern' parameter can't be cast to real_vector")
    __pattern = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __r = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__signal, signal, DT_REAL, X_CREATE)
        x_from_list(__pattern, pattern, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_corrr1d(ctypes.byref(_error_msg), ctypes.byref(__signal), ctypes.byref(__n), ctypes.byref(__pattern), ctypes.byref(__m), ctypes.byref(__r))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'corrr1d'")
        __r__r = list_from_x(__r)
        return __r__r
    finally:
        x_vector_clear(__signal)
        x_vector_clear(__pattern)
        x_vector_clear(__r)


_lib_alglib.alglib_corrr1dcircular.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_corrr1dcircular.restype = ctypes.c_int32
def corrr1dcircular(signal, m, pattern, n):
    pass
    if not is_real_vector(signal):
        raise ValueError("'signal' parameter can't be cast to real_vector")
    __signal = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(pattern):
        raise ValueError("'pattern' parameter can't be cast to real_vector")
    __pattern = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__signal, signal, DT_REAL, X_CREATE)
        x_from_list(__pattern, pattern, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_corrr1dcircular(ctypes.byref(_error_msg), ctypes.byref(__signal), ctypes.byref(__m), ctypes.byref(__pattern), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'corrr1dcircular'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__signal)
        x_vector_clear(__pattern)
        x_vector_clear(__c)


_lib_alglib.alglib_fhtr1d.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fhtr1d.restype = ctypes.c_int32
def fhtr1d(a, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fhtr1d(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fhtr1d'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_fhtr1dinv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fhtr1dinv.restype = ctypes.c_int32
def fhtr1dinv(a, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fhtr1dinv(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fhtr1dinv'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_gqgeneraterec.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgeneraterec.restype = ctypes.c_int32
def gqgeneraterec(alpha, beta, mu0, n):
    pass
    if not is_real_vector(alpha):
        raise ValueError("'alpha' parameter can't be cast to real_vector")
    __alpha = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(beta):
        raise ValueError("'beta' parameter can't be cast to real_vector")
    __beta = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __mu0 = ctypes.c_double(mu0)
    if __mu0.value!=mu0:
        raise ValueError("Error while converting 'mu0' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__alpha, alpha, DT_REAL, X_CREATE)
        x_from_list(__beta, beta, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgeneraterec(ctypes.byref(_error_msg), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__mu0), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgeneraterec'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__alpha)
        x_vector_clear(__beta)
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategausslobattorec.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategausslobattorec.restype = ctypes.c_int32
def gqgenerategausslobattorec(alpha, beta, mu0, a, b, n):
    pass
    if not is_real_vector(alpha):
        raise ValueError("'alpha' parameter can't be cast to real_vector")
    __alpha = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(beta):
        raise ValueError("'beta' parameter can't be cast to real_vector")
    __beta = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __mu0 = ctypes.c_double(mu0)
    if __mu0.value!=mu0:
        raise ValueError("Error while converting 'mu0' parameter to 'ctypes.c_double'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__alpha, alpha, DT_REAL, X_CREATE)
        x_from_list(__beta, beta, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategausslobattorec(ctypes.byref(_error_msg), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__mu0), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategausslobattorec'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__alpha)
        x_vector_clear(__beta)
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategaussradaurec.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategaussradaurec.restype = ctypes.c_int32
def gqgenerategaussradaurec(alpha, beta, mu0, a, n):
    pass
    if not is_real_vector(alpha):
        raise ValueError("'alpha' parameter can't be cast to real_vector")
    __alpha = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(beta):
        raise ValueError("'beta' parameter can't be cast to real_vector")
    __beta = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __mu0 = ctypes.c_double(mu0)
    if __mu0.value!=mu0:
        raise ValueError("Error while converting 'mu0' parameter to 'ctypes.c_double'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__alpha, alpha, DT_REAL, X_CREATE)
        x_from_list(__beta, beta, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategaussradaurec(ctypes.byref(_error_msg), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__mu0), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategaussradaurec'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__alpha)
        x_vector_clear(__beta)
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategausslegendre.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategausslegendre.restype = ctypes.c_int32
def gqgenerategausslegendre(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategausslegendre(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategausslegendre'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategaussjacobi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategaussjacobi.restype = ctypes.c_int32
def gqgenerategaussjacobi(n, alpha, beta):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategaussjacobi(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategaussjacobi'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategausslaguerre.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategausslaguerre.restype = ctypes.c_int32
def gqgenerategausslaguerre(n, alpha):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategausslaguerre(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__alpha), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategausslaguerre'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gqgenerategausshermite.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gqgenerategausshermite.restype = ctypes.c_int32
def gqgenerategausshermite(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gqgenerategausshermite(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gqgenerategausshermite'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__w = list_from_x(__w)
        return (__r__info, __r__x, __r__w)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__w)


_lib_alglib.alglib_gkqgeneraterec.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gkqgeneraterec.restype = ctypes.c_int32
def gkqgeneraterec(alpha, beta, mu0, n):
    pass
    if not is_real_vector(alpha):
        raise ValueError("'alpha' parameter can't be cast to real_vector")
    __alpha = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(beta):
        raise ValueError("'beta' parameter can't be cast to real_vector")
    __beta = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __mu0 = ctypes.c_double(mu0)
    if __mu0.value!=mu0:
        raise ValueError("Error while converting 'mu0' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wkronrod = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wgauss = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__alpha, alpha, DT_REAL, X_CREATE)
        x_from_list(__beta, beta, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gkqgeneraterec(ctypes.byref(_error_msg), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__mu0), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__wkronrod), ctypes.byref(__wgauss))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gkqgeneraterec'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__wkronrod = list_from_x(__wkronrod)
        __r__wgauss = list_from_x(__wgauss)
        return (__r__info, __r__x, __r__wkronrod, __r__wgauss)
    finally:
        x_vector_clear(__alpha)
        x_vector_clear(__beta)
        x_vector_clear(__x)
        x_vector_clear(__wkronrod)
        x_vector_clear(__wgauss)


_lib_alglib.alglib_gkqgenerategausslegendre.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gkqgenerategausslegendre.restype = ctypes.c_int32
def gkqgenerategausslegendre(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wkronrod = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wgauss = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gkqgenerategausslegendre(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__wkronrod), ctypes.byref(__wgauss))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gkqgenerategausslegendre'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__wkronrod = list_from_x(__wkronrod)
        __r__wgauss = list_from_x(__wgauss)
        return (__r__info, __r__x, __r__wkronrod, __r__wgauss)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__wkronrod)
        x_vector_clear(__wgauss)


_lib_alglib.alglib_gkqgenerategaussjacobi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gkqgenerategaussjacobi.restype = ctypes.c_int32
def gkqgenerategaussjacobi(n, alpha, beta):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wkronrod = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wgauss = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gkqgenerategaussjacobi(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__wkronrod), ctypes.byref(__wgauss))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gkqgenerategaussjacobi'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__wkronrod = list_from_x(__wkronrod)
        __r__wgauss = list_from_x(__wgauss)
        return (__r__info, __r__x, __r__wkronrod, __r__wgauss)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__wkronrod)
        x_vector_clear(__wgauss)


_lib_alglib.alglib_gkqlegendrecalc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gkqlegendrecalc.restype = ctypes.c_int32
def gkqlegendrecalc(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wkronrod = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wgauss = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gkqlegendrecalc(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__info), ctypes.byref(__x), ctypes.byref(__wkronrod), ctypes.byref(__wgauss))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gkqlegendrecalc'")
        __r__info = __info.value
        __r__x = list_from_x(__x)
        __r__wkronrod = list_from_x(__wkronrod)
        __r__wgauss = list_from_x(__wgauss)
        return (__r__info, __r__x, __r__wkronrod, __r__wgauss)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__wkronrod)
        x_vector_clear(__wgauss)


_lib_alglib.alglib_gkqlegendretbl.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_gkqlegendretbl.restype = ctypes.c_int32
def gkqlegendretbl(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wkronrod = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __wgauss = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __eps = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_gkqlegendretbl(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__wkronrod), ctypes.byref(__wgauss), ctypes.byref(__eps))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'gkqlegendretbl'")
        __r__x = list_from_x(__x)
        __r__wkronrod = list_from_x(__wkronrod)
        __r__wgauss = list_from_x(__wgauss)
        __r__eps = __eps.value
        return (__r__x, __r__wkronrod, __r__wgauss, __r__eps)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__wkronrod)
        x_vector_clear(__wgauss)




class x_autogkreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("terminationtype", c_ptrint_t),
        ("nfev", c_ptrint_t),
        ("nintervals", c_ptrint_t)
        ]




class autogkreport(object):
    def __init__(self):
        self.terminationtype = 0
        self.nfev = 0
        self.nintervals = 0


def x_autogkreport_zero_fields(x):
    x.terminationtype = 0
    x.nfev = 0
    x.nintervals = 0
    return




def x_autogkreport_clear(x):
    x_autogkreport_zero_fields(x)
    return




def x_from_autogkreport(x,v):
    x.terminationtype = int(v.terminationtype)
    x.nfev = int(v.nfev)
    x.nintervals = int(v.nintervals)
    return




def autogkreport_from_x(x):
    r = autogkreport()
    r.terminationtype = x.terminationtype
    r.nfev = x.nfev
    r.nintervals = x.nintervals
    return r


_lib_alglib.x_obj_free_autogkstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_autogkstate.restype = None
_lib_alglib.x_autogkstate_get_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_get_needf.restype = None
_lib_alglib.x_autogkstate_set_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_set_needf.restype = None
_lib_alglib.x_autogkstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_get_x.restype = None
_lib_alglib.x_autogkstate_set_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_set_x.restype = None
_lib_alglib.x_autogkstate_get_xminusa.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_get_xminusa.restype = None
_lib_alglib.x_autogkstate_set_xminusa.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_set_xminusa.restype = None
_lib_alglib.x_autogkstate_get_bminusx.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_get_bminusx.restype = None
_lib_alglib.x_autogkstate_set_bminusx.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_set_bminusx.restype = None
_lib_alglib.x_autogkstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_get_f.restype = None
_lib_alglib.x_autogkstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_autogkstate_set_f.restype = None


class autogkstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_autogkstate(self.ptr)
_lib_alglib.alglib_autogksmooth.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_autogksmooth.restype = ctypes.c_int32
def autogksmooth(a, b):
    pass
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_autogksmooth(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'autogksmooth'")
        __r__state = autogkstate(__state)
        return __r__state
    finally:
        pass


_lib_alglib.alglib_autogksmoothw.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_autogksmoothw.restype = ctypes.c_int32
def autogksmoothw(a, b, xwidth):
    pass
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __xwidth = ctypes.c_double(xwidth)
    if __xwidth.value!=xwidth:
        raise ValueError("Error while converting 'xwidth' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_autogksmoothw(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__xwidth), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'autogksmoothw'")
        __r__state = autogkstate(__state)
        return __r__state
    finally:
        pass


_lib_alglib.alglib_autogksingular.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_autogksingular.restype = ctypes.c_int32
def autogksingular(a, b, alpha, beta):
    pass
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __alpha = ctypes.c_double(alpha)
    if __alpha.value!=alpha:
        raise ValueError("Error while converting 'alpha' parameter to 'ctypes.c_double'")
    __beta = ctypes.c_double(beta)
    if __beta.value!=beta:
        raise ValueError("Error while converting 'beta' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_autogksingular(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__alpha), ctypes.byref(__beta), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'autogksingular'")
        __r__state = autogkstate(__state)
        return __r__state
    finally:
        pass




def autogkintegrate(state, func, param = None):
    # initialize temporaries
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_flag = ctypes.c_uint8()
    
    # initialize reverse communication variables
    _xc_x = ctypes.c_double()
    _xc_xminusa = ctypes.c_double()
    _xc_bminusx = ctypes.c_double()
    _xc_f = ctypes.c_double()
    
    # algorithm iterations
    while True:
        retval = _lib_alglib.alglib_autogkiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'autogkiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_autogkstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            _lib_alglib.x_autogkstate_get_x(state.ptr, ctypes.byref(_xc_x))
            _lib_alglib.x_autogkstate_get_xminusa(state.ptr, ctypes.byref(_xc_xminusa))
            _lib_alglib.x_autogkstate_get_bminusx(state.ptr, ctypes.byref(_xc_bminusx))
            _xc_f.value = func(_xc_x.value, _xc_xminusa.value, _xc_bminusx.value, param)
            _lib_alglib.x_autogkstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        raise RuntimeError("ALGLIB: unexpected error in 'autogkintegrate'")
    return


_lib_alglib.alglib_autogkresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_autogkresults.restype = ctypes.c_int32
def autogkresults(state):
    pass
    __state = state.ptr
    __v = ctypes.c_double(0)
    __rep = x_autogkreport()
    x_autogkreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_autogkresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__v), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'autogkresults'")
        __r__v = __v.value
        __r__rep = autogkreport_from_x(__rep)
        return (__r__v, __r__rep)
    finally:
        x_autogkreport_clear(__rep)


_lib_alglib.x_obj_free_idwinterpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_idwinterpolant.restype = None


class idwinterpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_idwinterpolant(self.ptr)
_lib_alglib.alglib_idwcalc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_idwcalc.restype = ctypes.c_int32
def idwcalc(z, x):
    pass
    __result = ctypes.c_double(0)
    __z = z.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_idwcalc(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__z), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'idwcalc'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_idwbuildmodifiedshepard.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_idwbuildmodifiedshepard.restype = ctypes.c_int32
def idwbuildmodifiedshepard(xy, n, nx, d, nq, nw):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nx = c_ptrint_t(nx)
    if __nx.value!=nx:
        raise ValueError("Error while converting 'nx' parameter to 'c_ptrint_t'")
    __d = c_ptrint_t(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'c_ptrint_t'")
    __nq = c_ptrint_t(nq)
    if __nq.value!=nq:
        raise ValueError("Error while converting 'nq' parameter to 'c_ptrint_t'")
    __nw = c_ptrint_t(nw)
    if __nw.value!=nw:
        raise ValueError("Error while converting 'nw' parameter to 'c_ptrint_t'")
    __z = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_idwbuildmodifiedshepard(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__nx), ctypes.byref(__d), ctypes.byref(__nq), ctypes.byref(__nw), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'idwbuildmodifiedshepard'")
        __r__z = idwinterpolant(__z)
        return __r__z
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_idwbuildmodifiedshepardr.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_idwbuildmodifiedshepardr.restype = ctypes.c_int32
def idwbuildmodifiedshepardr(xy, n, nx, r):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nx = c_ptrint_t(nx)
    if __nx.value!=nx:
        raise ValueError("Error while converting 'nx' parameter to 'c_ptrint_t'")
    __r = ctypes.c_double(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'ctypes.c_double'")
    __z = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_idwbuildmodifiedshepardr(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__nx), ctypes.byref(__r), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'idwbuildmodifiedshepardr'")
        __r__z = idwinterpolant(__z)
        return __r__z
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_idwbuildnoisy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_idwbuildnoisy.restype = ctypes.c_int32
def idwbuildnoisy(xy, n, nx, d, nq, nw):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __nx = c_ptrint_t(nx)
    if __nx.value!=nx:
        raise ValueError("Error while converting 'nx' parameter to 'c_ptrint_t'")
    __d = c_ptrint_t(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'c_ptrint_t'")
    __nq = c_ptrint_t(nq)
    if __nq.value!=nq:
        raise ValueError("Error while converting 'nq' parameter to 'c_ptrint_t'")
    __nw = c_ptrint_t(nw)
    if __nw.value!=nw:
        raise ValueError("Error while converting 'nw' parameter to 'c_ptrint_t'")
    __z = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_idwbuildnoisy(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__nx), ctypes.byref(__d), ctypes.byref(__nq), ctypes.byref(__nw), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'idwbuildnoisy'")
        __r__z = idwinterpolant(__z)
        return __r__z
    finally:
        x_matrix_clear(__xy)


_lib_alglib.x_obj_free_barycentricinterpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_barycentricinterpolant.restype = None


class barycentricinterpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_barycentricinterpolant(self.ptr)
_lib_alglib.alglib_barycentriccalc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentriccalc.restype = ctypes.c_int32
def barycentriccalc(b, t):
    pass
    __result = ctypes.c_double(0)
    __b = b.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentriccalc(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__b), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentriccalc'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_barycentricdiff1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricdiff1.restype = ctypes.c_int32
def barycentricdiff1(b, t):
    pass
    __b = b.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __f = ctypes.c_double(0)
    __df = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricdiff1(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__t), ctypes.byref(__f), ctypes.byref(__df))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricdiff1'")
        __r__f = __f.value
        __r__df = __df.value
        return (__r__f, __r__df)
    finally:
        pass


_lib_alglib.alglib_barycentricdiff2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricdiff2.restype = ctypes.c_int32
def barycentricdiff2(b, t):
    pass
    __b = b.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __f = ctypes.c_double(0)
    __df = ctypes.c_double(0)
    __d2f = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricdiff2(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__t), ctypes.byref(__f), ctypes.byref(__df), ctypes.byref(__d2f))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricdiff2'")
        __r__f = __f.value
        __r__df = __df.value
        __r__d2f = __d2f.value
        return (__r__f, __r__df, __r__d2f)
    finally:
        pass


_lib_alglib.alglib_barycentriclintransx.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentriclintransx.restype = ctypes.c_int32
def barycentriclintransx(b, ca, cb):
    pass
    __b = b.ptr
    __ca = ctypes.c_double(ca)
    if __ca.value!=ca:
        raise ValueError("Error while converting 'ca' parameter to 'ctypes.c_double'")
    __cb = ctypes.c_double(cb)
    if __cb.value!=cb:
        raise ValueError("Error while converting 'cb' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentriclintransx(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__ca), ctypes.byref(__cb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentriclintransx'")
        return
    finally:
        pass


_lib_alglib.alglib_barycentriclintransy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentriclintransy.restype = ctypes.c_int32
def barycentriclintransy(b, ca, cb):
    pass
    __b = b.ptr
    __ca = ctypes.c_double(ca)
    if __ca.value!=ca:
        raise ValueError("Error while converting 'ca' parameter to 'ctypes.c_double'")
    __cb = ctypes.c_double(cb)
    if __cb.value!=cb:
        raise ValueError("Error while converting 'cb' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentriclintransy(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__ca), ctypes.byref(__cb))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentriclintransy'")
        return
    finally:
        pass


_lib_alglib.alglib_barycentricunpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricunpack.restype = ctypes.c_int32
def barycentricunpack(b):
    pass
    __b = b.ptr
    __n = c_ptrint_t(0)
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __y = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __w = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricunpack(ctypes.byref(_error_msg), ctypes.byref(__b), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricunpack'")
        __r__n = __n.value
        __r__x = list_from_x(__x)
        __r__y = list_from_x(__y)
        __r__w = list_from_x(__w)
        return (__r__n, __r__x, __r__y, __r__w)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)


_lib_alglib.alglib_barycentricbuildxyw.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricbuildxyw.restype = ctypes.c_int32
def barycentricbuildxyw(x, y, w, n):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __b = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricbuildxyw(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricbuildxyw'")
        __r__b = barycentricinterpolant(__b)
        return __r__b
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)


_lib_alglib.alglib_barycentricbuildfloaterhormann.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricbuildfloaterhormann.restype = ctypes.c_int32
def barycentricbuildfloaterhormann(x, y, n, d):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __d = c_ptrint_t(d)
    if __d.value!=d:
        raise ValueError("Error while converting 'd' parameter to 'c_ptrint_t'")
    __b = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricbuildfloaterhormann(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__d), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricbuildfloaterhormann'")
        __r__b = barycentricinterpolant(__b)
        return __r__b
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_polynomialbar2cheb.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbar2cheb.restype = ctypes.c_int32
def polynomialbar2cheb(p, a, b):
    pass
    __p = p.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __t = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbar2cheb(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbar2cheb'")
        __r__t = list_from_x(__t)
        return __r__t
    finally:
        x_vector_clear(__t)


_lib_alglib.alglib_polynomialcheb2bar.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialcheb2bar.restype = ctypes.c_int32
def polynomialcheb2bar(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        t,n,a,b = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        t,a,b = functionargs
        n = safe_len("'polynomialcheb2bar': incorrect parameters",t)
    else:
        raise RuntimeError("Error while calling 'polynomialcheb2bar': function must have 3 or 4 parameters")
    if not is_real_vector(t):
        raise ValueError("'t' parameter can't be cast to real_vector")
    __t = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__t, t, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialcheb2bar(ctypes.byref(_error_msg), ctypes.byref(__t), ctypes.byref(__n), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialcheb2bar'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__t)


_lib_alglib.alglib_polynomialbar2pow.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbar2pow.restype = ctypes.c_int32
def polynomialbar2pow(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        p,c,s = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        p, = functionargs
        c = 0
        s = 1
    else:
        raise RuntimeError("Error while calling 'polynomialbar2pow': function must have 1 or 3 parameters")
    __p = p.ptr
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __s = ctypes.c_double(s)
    if __s.value!=s:
        raise ValueError("Error while converting 's' parameter to 'ctypes.c_double'")
    __a = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbar2pow(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__c), ctypes.byref(__s), ctypes.byref(__a))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbar2pow'")
        __r__a = list_from_x(__a)
        return __r__a
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_polynomialpow2bar.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialpow2bar.restype = ctypes.c_int32
def polynomialpow2bar(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,n,c,s = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        n = safe_len("'polynomialpow2bar': incorrect parameters",a)
        c = 0
        s = 1
    else:
        raise RuntimeError("Error while calling 'polynomialpow2bar': function must have 1 or 4 parameters")
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __s = ctypes.c_double(s)
    if __s.value!=s:
        raise ValueError("Error while converting 's' parameter to 'ctypes.c_double'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialpow2bar(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__c), ctypes.byref(__s), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialpow2bar'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__a)


_lib_alglib.alglib_polynomialbuild.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbuild.restype = ctypes.c_int32
def polynomialbuild(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'polynomialbuild': incorrect parameters",x)!=safe_len("'polynomialbuild': incorrect parameters",y):
            raise RuntimeError("Error while calling 'polynomialbuild': looks like one of arguments has wrong size")
        n = safe_len("'polynomialbuild': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'polynomialbuild': function must have 2 or 3 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbuild(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbuild'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_polynomialbuildeqdist.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbuildeqdist.restype = ctypes.c_int32
def polynomialbuildeqdist(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,b,y,n = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        a,b,y = functionargs
        n = safe_len("'polynomialbuildeqdist': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'polynomialbuildeqdist': function must have 3 or 4 parameters")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbuildeqdist(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbuildeqdist'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__y)


_lib_alglib.alglib_polynomialbuildcheb1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbuildcheb1.restype = ctypes.c_int32
def polynomialbuildcheb1(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,b,y,n = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        a,b,y = functionargs
        n = safe_len("'polynomialbuildcheb1': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'polynomialbuildcheb1': function must have 3 or 4 parameters")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbuildcheb1(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbuildcheb1'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__y)


_lib_alglib.alglib_polynomialbuildcheb2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialbuildcheb2.restype = ctypes.c_int32
def polynomialbuildcheb2(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        a,b,y,n = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        a,b,y = functionargs
        n = safe_len("'polynomialbuildcheb2': incorrect parameters",y)
    else:
        raise RuntimeError("Error while calling 'polynomialbuildcheb2': function must have 3 or 4 parameters")
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialbuildcheb2(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialbuildcheb2'")
        __r__p = barycentricinterpolant(__p)
        return __r__p
    finally:
        x_vector_clear(__y)


_lib_alglib.alglib_polynomialcalceqdist.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialcalceqdist.restype = ctypes.c_int32
def polynomialcalceqdist(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        a,b,f,n,t = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        a,b,f,t = functionargs
        n = safe_len("'polynomialcalceqdist': incorrect parameters",f)
    else:
        raise RuntimeError("Error while calling 'polynomialcalceqdist': function must have 4 or 5 parameters")
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(f):
        raise ValueError("'f' parameter can't be cast to real_vector")
    __f = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__f, f, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialcalceqdist(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__f), ctypes.byref(__n), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialcalceqdist'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__f)


_lib_alglib.alglib_polynomialcalccheb1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialcalccheb1.restype = ctypes.c_int32
def polynomialcalccheb1(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        a,b,f,n,t = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        a,b,f,t = functionargs
        n = safe_len("'polynomialcalccheb1': incorrect parameters",f)
    else:
        raise RuntimeError("Error while calling 'polynomialcalccheb1': function must have 4 or 5 parameters")
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(f):
        raise ValueError("'f' parameter can't be cast to real_vector")
    __f = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__f, f, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialcalccheb1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__f), ctypes.byref(__n), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialcalccheb1'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__f)


_lib_alglib.alglib_polynomialcalccheb2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialcalccheb2.restype = ctypes.c_int32
def polynomialcalccheb2(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        a,b,f,n,t = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        a,b,f,t = functionargs
        n = safe_len("'polynomialcalccheb2': incorrect parameters",f)
    else:
        raise RuntimeError("Error while calling 'polynomialcalccheb2': function must have 4 or 5 parameters")
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    if not is_real_vector(f):
        raise ValueError("'f' parameter can't be cast to real_vector")
    __f = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__f, f, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialcalccheb2(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__f), ctypes.byref(__n), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialcalccheb2'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__f)


_lib_alglib.x_obj_free_spline1dinterpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_spline1dinterpolant.restype = None


class spline1dinterpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_spline1dinterpolant(self.ptr)
_lib_alglib.alglib_spline1dbuildlinear.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dbuildlinear.restype = ctypes.c_int32
def spline1dbuildlinear(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dbuildlinear': incorrect parameters",x)!=safe_len("'spline1dbuildlinear': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dbuildlinear': looks like one of arguments has wrong size")
        n = safe_len("'spline1dbuildlinear': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dbuildlinear': function must have 2 or 3 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dbuildlinear(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dbuildlinear'")
        __r__c = spline1dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spline1dbuildcubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dbuildcubic.restype = ctypes.c_int32
def spline1dbuildcubic(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dbuildcubic': incorrect parameters",x)!=safe_len("'spline1dbuildcubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dbuildcubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dbuildcubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
    else:
        raise RuntimeError("Error while calling 'spline1dbuildcubic': function must have 2 or 7 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dbuildcubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dbuildcubic'")
        __r__c = spline1dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spline1dgriddiffcubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dgriddiffcubic.restype = ctypes.c_int32
def spline1dgriddiffcubic(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dgriddiffcubic': incorrect parameters",x)!=safe_len("'spline1dgriddiffcubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dgriddiffcubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dgriddiffcubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
    else:
        raise RuntimeError("Error while calling 'spline1dgriddiffcubic': function must have 2 or 7 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dgriddiffcubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__d))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dgriddiffcubic'")
        __r__d = list_from_x(__d)
        return __r__d
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__d)


_lib_alglib.alglib_spline1dgriddiff2cubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dgriddiff2cubic.restype = ctypes.c_int32
def spline1dgriddiff2cubic(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dgriddiff2cubic': incorrect parameters",x)!=safe_len("'spline1dgriddiff2cubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dgriddiff2cubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dgriddiff2cubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
    else:
        raise RuntimeError("Error while calling 'spline1dgriddiff2cubic': function must have 2 or 7 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    __d1 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __d2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dgriddiff2cubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__d1), ctypes.byref(__d2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dgriddiff2cubic'")
        __r__d1 = list_from_x(__d1)
        __r__d2 = list_from_x(__d2)
        return (__r__d1, __r__d2)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__d1)
        x_vector_clear(__d2)


_lib_alglib.alglib_spline1dconvcubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dconvcubic.restype = ctypes.c_int32
def spline1dconvcubic(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr,x2,n2 = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,x2 = functionargs
        if safe_len("'spline1dconvcubic': incorrect parameters",x)!=safe_len("'spline1dconvcubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dconvcubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dconvcubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
        n2 = safe_len("'spline1dconvcubic': incorrect parameters",x2)
    else:
        raise RuntimeError("Error while calling 'spline1dconvcubic': function must have 3 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    if not is_real_vector(x2):
        raise ValueError("'x2' parameter can't be cast to real_vector")
    __x2 = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n2 = c_ptrint_t(n2)
    if __n2.value!=n2:
        raise ValueError("Error while converting 'n2' parameter to 'c_ptrint_t'")
    __y2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__x2, x2, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dconvcubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__x2), ctypes.byref(__n2), ctypes.byref(__y2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dconvcubic'")
        __r__y2 = list_from_x(__y2)
        return __r__y2
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__x2)
        x_vector_clear(__y2)


_lib_alglib.alglib_spline1dconvdiffcubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dconvdiffcubic.restype = ctypes.c_int32
def spline1dconvdiffcubic(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr,x2,n2 = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,x2 = functionargs
        if safe_len("'spline1dconvdiffcubic': incorrect parameters",x)!=safe_len("'spline1dconvdiffcubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dconvdiffcubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dconvdiffcubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
        n2 = safe_len("'spline1dconvdiffcubic': incorrect parameters",x2)
    else:
        raise RuntimeError("Error while calling 'spline1dconvdiffcubic': function must have 3 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    if not is_real_vector(x2):
        raise ValueError("'x2' parameter can't be cast to real_vector")
    __x2 = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n2 = c_ptrint_t(n2)
    if __n2.value!=n2:
        raise ValueError("Error while converting 'n2' parameter to 'c_ptrint_t'")
    __y2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __d2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__x2, x2, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dconvdiffcubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__x2), ctypes.byref(__n2), ctypes.byref(__y2), ctypes.byref(__d2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dconvdiffcubic'")
        __r__y2 = list_from_x(__y2)
        __r__d2 = list_from_x(__d2)
        return (__r__y2, __r__d2)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__x2)
        x_vector_clear(__y2)
        x_vector_clear(__d2)


_lib_alglib.alglib_spline1dconvdiff2cubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dconvdiff2cubic.restype = ctypes.c_int32
def spline1dconvdiff2cubic(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,n,boundltype,boundl,boundrtype,boundr,x2,n2 = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,x2 = functionargs
        if safe_len("'spline1dconvdiff2cubic': incorrect parameters",x)!=safe_len("'spline1dconvdiff2cubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dconvdiff2cubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dconvdiff2cubic': incorrect parameters",x)
        boundltype = 0
        boundl = 0
        boundrtype = 0
        boundr = 0
        n2 = safe_len("'spline1dconvdiff2cubic': incorrect parameters",x2)
    else:
        raise RuntimeError("Error while calling 'spline1dconvdiff2cubic': function must have 3 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundltype = c_ptrint_t(boundltype)
    if __boundltype.value!=boundltype:
        raise ValueError("Error while converting 'boundltype' parameter to 'c_ptrint_t'")
    __boundl = ctypes.c_double(boundl)
    if __boundl.value!=boundl:
        raise ValueError("Error while converting 'boundl' parameter to 'ctypes.c_double'")
    __boundrtype = c_ptrint_t(boundrtype)
    if __boundrtype.value!=boundrtype:
        raise ValueError("Error while converting 'boundrtype' parameter to 'c_ptrint_t'")
    __boundr = ctypes.c_double(boundr)
    if __boundr.value!=boundr:
        raise ValueError("Error while converting 'boundr' parameter to 'ctypes.c_double'")
    if not is_real_vector(x2):
        raise ValueError("'x2' parameter can't be cast to real_vector")
    __x2 = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n2 = c_ptrint_t(n2)
    if __n2.value!=n2:
        raise ValueError("Error while converting 'n2' parameter to 'c_ptrint_t'")
    __y2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __d2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __dd2 = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__x2, x2, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dconvdiff2cubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundltype), ctypes.byref(__boundl), ctypes.byref(__boundrtype), ctypes.byref(__boundr), ctypes.byref(__x2), ctypes.byref(__n2), ctypes.byref(__y2), ctypes.byref(__d2), ctypes.byref(__dd2))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dconvdiff2cubic'")
        __r__y2 = list_from_x(__y2)
        __r__d2 = list_from_x(__d2)
        __r__dd2 = list_from_x(__dd2)
        return (__r__y2, __r__d2, __r__dd2)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__x2)
        x_vector_clear(__y2)
        x_vector_clear(__d2)
        x_vector_clear(__dd2)


_lib_alglib.alglib_spline1dbuildcatmullrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dbuildcatmullrom.restype = ctypes.c_int32
def spline1dbuildcatmullrom(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        x,y,n,boundtype,tension = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dbuildcatmullrom': incorrect parameters",x)!=safe_len("'spline1dbuildcatmullrom': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dbuildcatmullrom': looks like one of arguments has wrong size")
        n = safe_len("'spline1dbuildcatmullrom': incorrect parameters",x)
        boundtype = 0
        tension = 0
    else:
        raise RuntimeError("Error while calling 'spline1dbuildcatmullrom': function must have 2 or 5 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __boundtype = c_ptrint_t(boundtype)
    if __boundtype.value!=boundtype:
        raise ValueError("Error while converting 'boundtype' parameter to 'c_ptrint_t'")
    __tension = ctypes.c_double(tension)
    if __tension.value!=tension:
        raise ValueError("Error while converting 'tension' parameter to 'ctypes.c_double'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dbuildcatmullrom(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__boundtype), ctypes.byref(__tension), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dbuildcatmullrom'")
        __r__c = spline1dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spline1dbuildhermite.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dbuildhermite.restype = ctypes.c_int32
def spline1dbuildhermite(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        x,y,d,n = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,d = functionargs
        if safe_len("'spline1dbuildhermite': incorrect parameters",x)!=safe_len("'spline1dbuildhermite': incorrect parameters",y) or safe_len("'spline1dbuildhermite': incorrect parameters",x)!=safe_len("'spline1dbuildhermite': incorrect parameters",d):
            raise RuntimeError("Error while calling 'spline1dbuildhermite': looks like one of arguments has wrong size")
        n = safe_len("'spline1dbuildhermite': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dbuildhermite': function must have 3 or 4 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(d):
        raise ValueError("'d' parameter can't be cast to real_vector")
    __d = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__d, d, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dbuildhermite(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__d), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dbuildhermite'")
        __r__c = spline1dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__d)


_lib_alglib.alglib_spline1dbuildakima.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dbuildakima.restype = ctypes.c_int32
def spline1dbuildakima(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        x,y,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        x,y = functionargs
        if safe_len("'spline1dbuildakima': incorrect parameters",x)!=safe_len("'spline1dbuildakima': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dbuildakima': looks like one of arguments has wrong size")
        n = safe_len("'spline1dbuildakima': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dbuildakima': function must have 2 or 3 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dbuildakima(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dbuildakima'")
        __r__c = spline1dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_spline1dcalc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dcalc.restype = ctypes.c_int32
def spline1dcalc(c, x):
    pass
    __result = ctypes.c_double(0)
    __c = c.ptr
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dcalc(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dcalc'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_spline1ddiff.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1ddiff.restype = ctypes.c_int32
def spline1ddiff(c, x):
    pass
    __c = c.ptr
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __s = ctypes.c_double(0)
    __ds = ctypes.c_double(0)
    __d2s = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1ddiff(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__x), ctypes.byref(__s), ctypes.byref(__ds), ctypes.byref(__d2s))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1ddiff'")
        __r__s = __s.value
        __r__ds = __ds.value
        __r__d2s = __d2s.value
        return (__r__s, __r__ds, __r__d2s)
    finally:
        pass


_lib_alglib.alglib_spline1dunpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dunpack.restype = ctypes.c_int32
def spline1dunpack(c):
    pass
    __c = c.ptr
    __n = c_ptrint_t(0)
    __tbl = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dunpack(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__tbl))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dunpack'")
        __r__n = __n.value
        __r__tbl = listlist_from_x(__tbl)
        return (__r__n, __r__tbl)
    finally:
        x_matrix_clear(__tbl)


_lib_alglib.alglib_spline1dlintransx.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dlintransx.restype = ctypes.c_int32
def spline1dlintransx(c, a, b):
    pass
    __c = c.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dlintransx(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dlintransx'")
        return
    finally:
        pass


_lib_alglib.alglib_spline1dlintransy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dlintransy.restype = ctypes.c_int32
def spline1dlintransy(c, a, b):
    pass
    __c = c.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dlintransy(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dlintransy'")
        return
    finally:
        pass


_lib_alglib.alglib_spline1dintegrate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dintegrate.restype = ctypes.c_int32
def spline1dintegrate(c, x):
    pass
    __result = ctypes.c_double(0)
    __c = c.ptr
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dintegrate(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dintegrate'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.x_obj_free_minlmstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_minlmstate.restype = None
_lib_alglib.x_minlmstate_get_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_needf.restype = None
_lib_alglib.x_minlmstate_set_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_needf.restype = None
_lib_alglib.x_minlmstate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_needfg.restype = None
_lib_alglib.x_minlmstate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_needfg.restype = None
_lib_alglib.x_minlmstate_get_needfgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_needfgh.restype = None
_lib_alglib.x_minlmstate_set_needfgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_needfgh.restype = None
_lib_alglib.x_minlmstate_get_needfi.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_needfi.restype = None
_lib_alglib.x_minlmstate_set_needfi.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_needfi.restype = None
_lib_alglib.x_minlmstate_get_needfij.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_needfij.restype = None
_lib_alglib.x_minlmstate_set_needfij.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_needfij.restype = None
_lib_alglib.x_minlmstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_xupdated.restype = None
_lib_alglib.x_minlmstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_xupdated.restype = None
_lib_alglib.x_minlmstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_f.restype = None
_lib_alglib.x_minlmstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_set_f.restype = None
_lib_alglib.x_minlmstate_get_fi.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_fi.restype = None
_lib_alglib.x_minlmstate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_g.restype = None
_lib_alglib.x_minlmstate_get_h.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_h.restype = None
_lib_alglib.x_minlmstate_get_j.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_j.restype = None
_lib_alglib.x_minlmstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minlmstate_get_x.restype = None


class minlmstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_minlmstate(self.ptr)


class x_minlmreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("iterationscount", c_ptrint_t),
        ("terminationtype", c_ptrint_t),
        ("nfunc", c_ptrint_t),
        ("njac", c_ptrint_t),
        ("ngrad", c_ptrint_t),
        ("nhess", c_ptrint_t),
        ("ncholesky", c_ptrint_t)
        ]




class minlmreport(object):
    def __init__(self):
        self.iterationscount = 0
        self.terminationtype = 0
        self.nfunc = 0
        self.njac = 0
        self.ngrad = 0
        self.nhess = 0
        self.ncholesky = 0


def x_minlmreport_zero_fields(x):
    x.iterationscount = 0
    x.terminationtype = 0
    x.nfunc = 0
    x.njac = 0
    x.ngrad = 0
    x.nhess = 0
    x.ncholesky = 0
    return




def x_minlmreport_clear(x):
    x_minlmreport_zero_fields(x)
    return




def x_from_minlmreport(x,v):
    x.iterationscount = int(v.iterationscount)
    x.terminationtype = int(v.terminationtype)
    x.nfunc = int(v.nfunc)
    x.njac = int(v.njac)
    x.ngrad = int(v.ngrad)
    x.nhess = int(v.nhess)
    x.ncholesky = int(v.ncholesky)
    return




def minlmreport_from_x(x):
    r = minlmreport()
    r.iterationscount = x.iterationscount
    r.terminationtype = x.terminationtype
    r.nfunc = x.nfunc
    r.njac = x.njac
    r.ngrad = x.ngrad
    r.nhess = x.nhess
    r.ncholesky = x.ncholesky
    return r


_lib_alglib.alglib_minlmcreatevj.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatevj.restype = ctypes.c_int32
def minlmcreatevj(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'minlmcreatevj': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatevj': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatevj(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatevj'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmcreatev.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatev.restype = ctypes.c_int32
def minlmcreatev(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        n,m,x,diffstep = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        m,x,diffstep = functionargs
        n = safe_len("'minlmcreatev': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatev': function must have 3 or 4 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __diffstep = ctypes.c_double(diffstep)
    if __diffstep.value!=diffstep:
        raise ValueError("Error while converting 'diffstep' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatev(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__diffstep), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatev'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmcreatefgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatefgh.restype = ctypes.c_int32
def minlmcreatefgh(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        n,x = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'minlmcreatefgh': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatefgh': function must have 1 or 2 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatefgh(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatefgh'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmcreatevgj.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatevgj.restype = ctypes.c_int32
def minlmcreatevgj(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'minlmcreatevgj': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatevgj': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatevgj(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatevgj'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmcreatefgj.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatefgj.restype = ctypes.c_int32
def minlmcreatefgj(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'minlmcreatefgj': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatefgj': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatefgj(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatefgj'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmcreatefj.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmcreatefj.restype = ctypes.c_int32
def minlmcreatefj(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'minlmcreatefj': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minlmcreatefj': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmcreatefj(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmcreatefj'")
        __r__state = minlmstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minlmsetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmsetcond.restype = ctypes.c_int32
def minlmsetcond(state, epsg, epsf, epsx, maxits):
    pass
    __state = state.ptr
    __epsg = ctypes.c_double(epsg)
    if __epsg.value!=epsg:
        raise ValueError("Error while converting 'epsg' parameter to 'ctypes.c_double'")
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmsetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsg), ctypes.byref(__epsf), ctypes.byref(__epsx), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmsetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_minlmsetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmsetxrep.restype = ctypes.c_int32
def minlmsetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmsetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmsetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_minlmsetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmsetstpmax.restype = ctypes.c_int32
def minlmsetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmsetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmsetstpmax'")
        return
    finally:
        pass


_lib_alglib.alglib_minlmsetacctype.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmsetacctype.restype = ctypes.c_int32
def minlmsetacctype(state, acctype):
    pass
    __state = state.ptr
    __acctype = c_ptrint_t(acctype)
    if __acctype.value!=acctype:
        raise ValueError("Error while converting 'acctype' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmsetacctype(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__acctype))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmsetacctype'")
        return
    finally:
        pass




def minlmoptimize_v(state, fvec, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlmstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_fi  = x_vector()
    _lib_alglib.x_minlmstate_get_fi(state.ptr, ctypes.byref(_xc_fi))
    _py_fi = create_real_vector(_xc_fi.cnt)
    while True:
        retval = _lib_alglib.alglib_minlmiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlmstate_get_needfi(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            fvec(_py_x, _py_fi, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlmstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)")
    return


def minlmoptimize_vj(state, fvec, jac, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlmstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_fi  = x_vector()
    _lib_alglib.x_minlmstate_get_fi(state.ptr, ctypes.byref(_xc_fi))
    _py_fi = create_real_vector(_xc_fi.cnt)
    _xc_j  = x_matrix()
    _lib_alglib.x_minlmstate_get_j(state.ptr, ctypes.byref(_xc_j))
    _py_j = create_real_matrix(_xc_j.rows,_xc_j.cols)
    while True:
        retval = _lib_alglib.alglib_minlmiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlmstate_get_needfi(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            fvec(_py_x, _py_fi, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_needfij(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            jac(_py_x, _py_fi, _py_j, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_j, _py_j, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlmstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)")
    return


def minlmoptimize_fgh(state, func, grad, hess, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlmstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_minlmstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    _xc_h  = x_matrix()
    _lib_alglib.x_minlmstate_get_h(state.ptr, ctypes.byref(_xc_h))
    _py_h = create_real_matrix(_xc_h.rows,_xc_h.cols)
    while True:
        retval = _lib_alglib.alglib_minlmiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlmstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = func(_py_x, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_minlmstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_needfgh(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = hess(_py_x, _py_g, _py_h, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_h, _py_h, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlmstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)")
    return


def minlmoptimize_fj(state, func, jac, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlmstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_fi  = x_vector()
    _lib_alglib.x_minlmstate_get_fi(state.ptr, ctypes.byref(_xc_fi))
    _py_fi = create_real_vector(_xc_fi.cnt)
    _xc_j  = x_matrix()
    _lib_alglib.x_minlmstate_get_j(state.ptr, ctypes.byref(_xc_j))
    _py_j = create_real_matrix(_xc_j.rows,_xc_j.cols)
    while True:
        retval = _lib_alglib.alglib_minlmiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlmstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = func(_py_x, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_minlmstate_get_needfij(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            jac(_py_x, _py_fi, _py_j, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_j, _py_j, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlmstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)")
    return


def minlmoptimize_fgj(state, func, grad, jac, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minlmstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_minlmstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    _xc_fi  = x_vector()
    _lib_alglib.x_minlmstate_get_fi(state.ptr, ctypes.byref(_xc_fi))
    _py_fi = create_real_vector(_xc_fi.cnt)
    _xc_j  = x_matrix()
    _lib_alglib.x_minlmstate_get_j(state.ptr, ctypes.byref(_xc_j))
    _py_j = create_real_matrix(_xc_j.rows,_xc_j.cols)
    while True:
        retval = _lib_alglib.alglib_minlmiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minlmstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = func(_py_x, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_minlmstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_minlmstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_needfij(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            jac(_py_x, _py_fi, _py_j, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_j, _py_j, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minlmstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minlmstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_minlmresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmresults.restype = ctypes.c_int32
def minlmresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minlmreport()
    x_minlmreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmresults'")
        __r__x = list_from_x(__x)
        __r__rep = minlmreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minlmreport_clear(__rep)


_lib_alglib.alglib_minlmresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmresultsbuf.restype = ctypes.c_int32
def minlmresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minlmreport()
    x_minlmreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_minlmreport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = minlmreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minlmreport_clear(__rep)


_lib_alglib.alglib_minlmrestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minlmrestartfrom.restype = ctypes.c_int32
def minlmrestartfrom(state, x):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minlmrestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minlmrestartfrom'")
        return
    finally:
        x_vector_clear(__x)




class x_polynomialfitreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("taskrcond", ctypes.c_double),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("maxerror", ctypes.c_double)
        ]




class polynomialfitreport(object):
    def __init__(self):
        self.taskrcond = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.maxerror = 0


def x_polynomialfitreport_zero_fields(x):
    x.taskrcond = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.maxerror = 0
    return




def x_polynomialfitreport_clear(x):
    x_polynomialfitreport_zero_fields(x)
    return




def x_from_polynomialfitreport(x,v):
    x.taskrcond = float(v.taskrcond)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.maxerror = float(v.maxerror)
    return




def polynomialfitreport_from_x(x):
    r = polynomialfitreport()
    r.taskrcond = x.taskrcond
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.maxerror = x.maxerror
    return r




class x_barycentricfitreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("taskrcond", ctypes.c_double),
        ("dbest", c_ptrint_t),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("maxerror", ctypes.c_double)
        ]




class barycentricfitreport(object):
    def __init__(self):
        self.taskrcond = 0
        self.dbest = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.maxerror = 0


def x_barycentricfitreport_zero_fields(x):
    x.taskrcond = 0
    x.dbest = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.maxerror = 0
    return




def x_barycentricfitreport_clear(x):
    x_barycentricfitreport_zero_fields(x)
    return




def x_from_barycentricfitreport(x,v):
    x.taskrcond = float(v.taskrcond)
    x.dbest = int(v.dbest)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.maxerror = float(v.maxerror)
    return




def barycentricfitreport_from_x(x):
    r = barycentricfitreport()
    r.taskrcond = x.taskrcond
    r.dbest = x.dbest
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.maxerror = x.maxerror
    return r




class x_spline1dfitreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("taskrcond", ctypes.c_double),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("maxerror", ctypes.c_double)
        ]




class spline1dfitreport(object):
    def __init__(self):
        self.taskrcond = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.maxerror = 0


def x_spline1dfitreport_zero_fields(x):
    x.taskrcond = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.maxerror = 0
    return




def x_spline1dfitreport_clear(x):
    x_spline1dfitreport_zero_fields(x)
    return




def x_from_spline1dfitreport(x,v):
    x.taskrcond = float(v.taskrcond)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.maxerror = float(v.maxerror)
    return




def spline1dfitreport_from_x(x):
    r = spline1dfitreport()
    r.taskrcond = x.taskrcond
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.maxerror = x.maxerror
    return r




class x_lsfitreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("taskrcond", ctypes.c_double),
        ("rmserror", ctypes.c_double),
        ("avgerror", ctypes.c_double),
        ("avgrelerror", ctypes.c_double),
        ("maxerror", ctypes.c_double)
        ]




class lsfitreport(object):
    def __init__(self):
        self.taskrcond = 0
        self.rmserror = 0
        self.avgerror = 0
        self.avgrelerror = 0
        self.maxerror = 0


def x_lsfitreport_zero_fields(x):
    x.taskrcond = 0
    x.rmserror = 0
    x.avgerror = 0
    x.avgrelerror = 0
    x.maxerror = 0
    return




def x_lsfitreport_clear(x):
    x_lsfitreport_zero_fields(x)
    return




def x_from_lsfitreport(x,v):
    x.taskrcond = float(v.taskrcond)
    x.rmserror = float(v.rmserror)
    x.avgerror = float(v.avgerror)
    x.avgrelerror = float(v.avgrelerror)
    x.maxerror = float(v.maxerror)
    return




def lsfitreport_from_x(x):
    r = lsfitreport()
    r.taskrcond = x.taskrcond
    r.rmserror = x.rmserror
    r.avgerror = x.avgerror
    r.avgrelerror = x.avgrelerror
    r.maxerror = x.maxerror
    return r


_lib_alglib.x_obj_free_lsfitstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_lsfitstate.restype = None
_lib_alglib.x_lsfitstate_get_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_needf.restype = None
_lib_alglib.x_lsfitstate_set_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_set_needf.restype = None
_lib_alglib.x_lsfitstate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_needfg.restype = None
_lib_alglib.x_lsfitstate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_set_needfg.restype = None
_lib_alglib.x_lsfitstate_get_needfgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_needfgh.restype = None
_lib_alglib.x_lsfitstate_set_needfgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_set_needfgh.restype = None
_lib_alglib.x_lsfitstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_xupdated.restype = None
_lib_alglib.x_lsfitstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_set_xupdated.restype = None
_lib_alglib.x_lsfitstate_get_c.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_c.restype = None
_lib_alglib.x_lsfitstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_f.restype = None
_lib_alglib.x_lsfitstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_set_f.restype = None
_lib_alglib.x_lsfitstate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_g.restype = None
_lib_alglib.x_lsfitstate_get_h.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_h.restype = None
_lib_alglib.x_lsfitstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_lsfitstate_get_x.restype = None


class lsfitstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_lsfitstate(self.ptr)
_lib_alglib.alglib_polynomialfit.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialfit.restype = ctypes.c_int32
def polynomialfit(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        x,y,n,m = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,m = functionargs
        if safe_len("'polynomialfit': incorrect parameters",x)!=safe_len("'polynomialfit': incorrect parameters",y):
            raise RuntimeError("Error while calling 'polynomialfit': looks like one of arguments has wrong size")
        n = safe_len("'polynomialfit': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'polynomialfit': function must have 3 or 4 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __p = ctypes.c_void_p(0)
    __rep = x_polynomialfitreport()
    x_polynomialfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialfit(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__p), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialfit'")
        __r__info = __info.value
        __r__p = barycentricinterpolant(__p)
        __r__rep = polynomialfitreport_from_x(__rep)
        return (__r__info, __r__p, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_polynomialfitreport_clear(__rep)


_lib_alglib.alglib_polynomialfitwc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_polynomialfitwc.restype = ctypes.c_int32
def polynomialfitwc(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,w,n,xc,yc,dc,k,m = functionargs
    elif len(functionargs)==7:
        __friendly_form = True
        x,y,w,xc,yc,dc,m = functionargs
        if safe_len("'polynomialfitwc': incorrect parameters",x)!=safe_len("'polynomialfitwc': incorrect parameters",y) or safe_len("'polynomialfitwc': incorrect parameters",x)!=safe_len("'polynomialfitwc': incorrect parameters",w):
            raise RuntimeError("Error while calling 'polynomialfitwc': looks like one of arguments has wrong size")
        n = safe_len("'polynomialfitwc': incorrect parameters",x)
        if safe_len("'polynomialfitwc': incorrect parameters",xc)!=safe_len("'polynomialfitwc': incorrect parameters",yc) or safe_len("'polynomialfitwc': incorrect parameters",xc)!=safe_len("'polynomialfitwc': incorrect parameters",dc):
            raise RuntimeError("Error while calling 'polynomialfitwc': looks like one of arguments has wrong size")
        k = safe_len("'polynomialfitwc': incorrect parameters",xc)
    else:
        raise RuntimeError("Error while calling 'polynomialfitwc': function must have 7 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(xc):
        raise ValueError("'xc' parameter can't be cast to real_vector")
    __xc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(yc):
        raise ValueError("'yc' parameter can't be cast to real_vector")
    __yc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(dc):
        raise ValueError("'dc' parameter can't be cast to int_vector")
    __dc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __p = ctypes.c_void_p(0)
    __rep = x_polynomialfitreport()
    x_polynomialfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__xc, xc, DT_REAL, X_CREATE)
        x_from_list(__yc, yc, DT_REAL, X_CREATE)
        x_from_list(__dc, dc, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_polynomialfitwc(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__xc), ctypes.byref(__yc), ctypes.byref(__dc), ctypes.byref(__k), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__p), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'polynomialfitwc'")
        __r__info = __info.value
        __r__p = barycentricinterpolant(__p)
        __r__rep = polynomialfitreport_from_x(__rep)
        return (__r__info, __r__p, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__xc)
        x_vector_clear(__yc)
        x_vector_clear(__dc)
        x_polynomialfitreport_clear(__rep)


_lib_alglib.alglib_barycentricfitfloaterhormannwc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricfitfloaterhormannwc.restype = ctypes.c_int32
def barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, k, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(xc):
        raise ValueError("'xc' parameter can't be cast to real_vector")
    __xc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(yc):
        raise ValueError("'yc' parameter can't be cast to real_vector")
    __yc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(dc):
        raise ValueError("'dc' parameter can't be cast to int_vector")
    __dc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __b = ctypes.c_void_p(0)
    __rep = x_barycentricfitreport()
    x_barycentricfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__xc, xc, DT_REAL, X_CREATE)
        x_from_list(__yc, yc, DT_REAL, X_CREATE)
        x_from_list(__dc, dc, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricfitfloaterhormannwc(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__xc), ctypes.byref(__yc), ctypes.byref(__dc), ctypes.byref(__k), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__b), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricfitfloaterhormannwc'")
        __r__info = __info.value
        __r__b = barycentricinterpolant(__b)
        __r__rep = barycentricfitreport_from_x(__rep)
        return (__r__info, __r__b, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__xc)
        x_vector_clear(__yc)
        x_vector_clear(__dc)
        x_barycentricfitreport_clear(__rep)


_lib_alglib.alglib_barycentricfitfloaterhormann.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_barycentricfitfloaterhormann.restype = ctypes.c_int32
def barycentricfitfloaterhormann(x, y, n, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __b = ctypes.c_void_p(0)
    __rep = x_barycentricfitreport()
    x_barycentricfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_barycentricfitfloaterhormann(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__b), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'barycentricfitfloaterhormann'")
        __r__info = __info.value
        __r__b = barycentricinterpolant(__b)
        __r__rep = barycentricfitreport_from_x(__rep)
        return (__r__info, __r__b, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_barycentricfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfitpenalized.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfitpenalized.restype = ctypes.c_int32
def spline1dfitpenalized(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        x,y,n,m,rho = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        x,y,m,rho = functionargs
        if safe_len("'spline1dfitpenalized': incorrect parameters",x)!=safe_len("'spline1dfitpenalized': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dfitpenalized': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfitpenalized': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dfitpenalized': function must have 4 or 5 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __rho = ctypes.c_double(rho)
    if __rho.value!=rho:
        raise ValueError("Error while converting 'rho' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfitpenalized(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__rho), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfitpenalized'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfitpenalizedw.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfitpenalizedw.restype = ctypes.c_int32
def spline1dfitpenalizedw(*functionargs):
    if len(functionargs)==6:
        __friendly_form = False
        x,y,w,n,m,rho = functionargs
    elif len(functionargs)==5:
        __friendly_form = True
        x,y,w,m,rho = functionargs
        if safe_len("'spline1dfitpenalizedw': incorrect parameters",x)!=safe_len("'spline1dfitpenalizedw': incorrect parameters",y) or safe_len("'spline1dfitpenalizedw': incorrect parameters",x)!=safe_len("'spline1dfitpenalizedw': incorrect parameters",w):
            raise RuntimeError("Error while calling 'spline1dfitpenalizedw': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfitpenalizedw': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dfitpenalizedw': function must have 5 or 6 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __rho = ctypes.c_double(rho)
    if __rho.value!=rho:
        raise ValueError("Error while converting 'rho' parameter to 'ctypes.c_double'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfitpenalizedw(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__rho), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfitpenalizedw'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfitcubicwc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfitcubicwc.restype = ctypes.c_int32
def spline1dfitcubicwc(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,w,n,xc,yc,dc,k,m = functionargs
    elif len(functionargs)==7:
        __friendly_form = True
        x,y,w,xc,yc,dc,m = functionargs
        if safe_len("'spline1dfitcubicwc': incorrect parameters",x)!=safe_len("'spline1dfitcubicwc': incorrect parameters",y) or safe_len("'spline1dfitcubicwc': incorrect parameters",x)!=safe_len("'spline1dfitcubicwc': incorrect parameters",w):
            raise RuntimeError("Error while calling 'spline1dfitcubicwc': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfitcubicwc': incorrect parameters",x)
        if safe_len("'spline1dfitcubicwc': incorrect parameters",xc)!=safe_len("'spline1dfitcubicwc': incorrect parameters",yc) or safe_len("'spline1dfitcubicwc': incorrect parameters",xc)!=safe_len("'spline1dfitcubicwc': incorrect parameters",dc):
            raise RuntimeError("Error while calling 'spline1dfitcubicwc': looks like one of arguments has wrong size")
        k = safe_len("'spline1dfitcubicwc': incorrect parameters",xc)
    else:
        raise RuntimeError("Error while calling 'spline1dfitcubicwc': function must have 7 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(xc):
        raise ValueError("'xc' parameter can't be cast to real_vector")
    __xc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(yc):
        raise ValueError("'yc' parameter can't be cast to real_vector")
    __yc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(dc):
        raise ValueError("'dc' parameter can't be cast to int_vector")
    __dc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__xc, xc, DT_REAL, X_CREATE)
        x_from_list(__yc, yc, DT_REAL, X_CREATE)
        x_from_list(__dc, dc, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfitcubicwc(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__xc), ctypes.byref(__yc), ctypes.byref(__dc), ctypes.byref(__k), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfitcubicwc'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__xc)
        x_vector_clear(__yc)
        x_vector_clear(__dc)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfithermitewc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfithermitewc.restype = ctypes.c_int32
def spline1dfithermitewc(*functionargs):
    if len(functionargs)==9:
        __friendly_form = False
        x,y,w,n,xc,yc,dc,k,m = functionargs
    elif len(functionargs)==7:
        __friendly_form = True
        x,y,w,xc,yc,dc,m = functionargs
        if safe_len("'spline1dfithermitewc': incorrect parameters",x)!=safe_len("'spline1dfithermitewc': incorrect parameters",y) or safe_len("'spline1dfithermitewc': incorrect parameters",x)!=safe_len("'spline1dfithermitewc': incorrect parameters",w):
            raise RuntimeError("Error while calling 'spline1dfithermitewc': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfithermitewc': incorrect parameters",x)
        if safe_len("'spline1dfithermitewc': incorrect parameters",xc)!=safe_len("'spline1dfithermitewc': incorrect parameters",yc) or safe_len("'spline1dfithermitewc': incorrect parameters",xc)!=safe_len("'spline1dfithermitewc': incorrect parameters",dc):
            raise RuntimeError("Error while calling 'spline1dfithermitewc': looks like one of arguments has wrong size")
        k = safe_len("'spline1dfithermitewc': incorrect parameters",xc)
    else:
        raise RuntimeError("Error while calling 'spline1dfithermitewc': function must have 7 or 9 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(xc):
        raise ValueError("'xc' parameter can't be cast to real_vector")
    __xc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(yc):
        raise ValueError("'yc' parameter can't be cast to real_vector")
    __yc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(dc):
        raise ValueError("'dc' parameter can't be cast to int_vector")
    __dc = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__xc, xc, DT_REAL, X_CREATE)
        x_from_list(__yc, yc, DT_REAL, X_CREATE)
        x_from_list(__dc, dc, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfithermitewc(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__n), ctypes.byref(__xc), ctypes.byref(__yc), ctypes.byref(__dc), ctypes.byref(__k), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfithermitewc'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__xc)
        x_vector_clear(__yc)
        x_vector_clear(__dc)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfitcubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfitcubic.restype = ctypes.c_int32
def spline1dfitcubic(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        x,y,n,m = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,m = functionargs
        if safe_len("'spline1dfitcubic': incorrect parameters",x)!=safe_len("'spline1dfitcubic': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dfitcubic': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfitcubic': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dfitcubic': function must have 3 or 4 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfitcubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfitcubic'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_spline1dfithermite.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline1dfithermite.restype = ctypes.c_int32
def spline1dfithermite(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        x,y,n,m = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,m = functionargs
        if safe_len("'spline1dfithermite': incorrect parameters",x)!=safe_len("'spline1dfithermite': incorrect parameters",y):
            raise RuntimeError("Error while calling 'spline1dfithermite': looks like one of arguments has wrong size")
        n = safe_len("'spline1dfithermite': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'spline1dfithermite': function must have 3 or 4 parameters")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __s = ctypes.c_void_p(0)
    __rep = x_spline1dfitreport()
    x_spline1dfitreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline1dfithermite(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__s), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline1dfithermite'")
        __r__info = __info.value
        __r__s = spline1dinterpolant(__s)
        __r__rep = spline1dfitreport_from_x(__rep)
        return (__r__info, __r__s, __r__rep)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_spline1dfitreport_clear(__rep)


_lib_alglib.alglib_lsfitlinearw.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitlinearw.restype = ctypes.c_int32
def lsfitlinearw(*functionargs):
    if len(functionargs)==5:
        __friendly_form = False
        y,w,fmatrix,n,m = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        y,w,fmatrix = functionargs
        if safe_len("'lsfitlinearw': incorrect parameters",y)!=safe_len("'lsfitlinearw': incorrect parameters",w) or safe_len("'lsfitlinearw': incorrect parameters",y)!=safe_rows("'lsfitlinearw': incorrect parameters",fmatrix):
            raise RuntimeError("Error while calling 'lsfitlinearw': looks like one of arguments has wrong size")
        n = safe_len("'lsfitlinearw': incorrect parameters",y)
        m = safe_cols("'lsfitlinearw': incorrect parameters",fmatrix)
    else:
        raise RuntimeError("Error while calling 'lsfitlinearw': function must have 3 or 5 parameters")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(fmatrix):
        raise ValueError("'fmatrix' parameter can't be cast to real_matrix")
    __fmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_lsfitreport()
    x_lsfitreport_zero_fields(__rep)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_listlist(__fmatrix, fmatrix, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitlinearw(ctypes.byref(_error_msg), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__fmatrix), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitlinearw'")
        __r__info = __info.value
        __r__c = list_from_x(__c)
        __r__rep = lsfitreport_from_x(__rep)
        return (__r__info, __r__c, __r__rep)
    finally:
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_matrix_clear(__fmatrix)
        x_vector_clear(__c)
        x_lsfitreport_clear(__rep)


_lib_alglib.alglib_lsfitlinearwc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitlinearwc.restype = ctypes.c_int32
def lsfitlinearwc(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        y,w,fmatrix,cmatrix,n,m,k = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        y,w,fmatrix,cmatrix = functionargs
        if safe_len("'lsfitlinearwc': incorrect parameters",y)!=safe_len("'lsfitlinearwc': incorrect parameters",w) or safe_len("'lsfitlinearwc': incorrect parameters",y)!=safe_rows("'lsfitlinearwc': incorrect parameters",fmatrix):
            raise RuntimeError("Error while calling 'lsfitlinearwc': looks like one of arguments has wrong size")
        n = safe_len("'lsfitlinearwc': incorrect parameters",y)
        if safe_cols("'lsfitlinearwc': incorrect parameters",fmatrix)!=safe_cols("'lsfitlinearwc': incorrect parameters",cmatrix)-1:
            raise RuntimeError("Error while calling 'lsfitlinearwc': looks like one of arguments has wrong size")
        m = safe_cols("'lsfitlinearwc': incorrect parameters",fmatrix)
        k = safe_rows("'lsfitlinearwc': incorrect parameters",cmatrix)
    else:
        raise RuntimeError("Error while calling 'lsfitlinearwc': function must have 4 or 7 parameters")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(fmatrix):
        raise ValueError("'fmatrix' parameter can't be cast to real_matrix")
    __fmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(cmatrix):
        raise ValueError("'cmatrix' parameter can't be cast to real_matrix")
    __cmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_lsfitreport()
    x_lsfitreport_zero_fields(__rep)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_listlist(__fmatrix, fmatrix, DT_REAL, X_CREATE)
        x_from_listlist(__cmatrix, cmatrix, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitlinearwc(ctypes.byref(_error_msg), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__fmatrix), ctypes.byref(__cmatrix), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitlinearwc'")
        __r__info = __info.value
        __r__c = list_from_x(__c)
        __r__rep = lsfitreport_from_x(__rep)
        return (__r__info, __r__c, __r__rep)
    finally:
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_matrix_clear(__fmatrix)
        x_matrix_clear(__cmatrix)
        x_vector_clear(__c)
        x_lsfitreport_clear(__rep)


_lib_alglib.alglib_lsfitlinear.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitlinear.restype = ctypes.c_int32
def lsfitlinear(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        y,fmatrix,n,m = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        y,fmatrix = functionargs
        if safe_len("'lsfitlinear': incorrect parameters",y)!=safe_rows("'lsfitlinear': incorrect parameters",fmatrix):
            raise RuntimeError("Error while calling 'lsfitlinear': looks like one of arguments has wrong size")
        n = safe_len("'lsfitlinear': incorrect parameters",y)
        m = safe_cols("'lsfitlinear': incorrect parameters",fmatrix)
    else:
        raise RuntimeError("Error while calling 'lsfitlinear': function must have 2 or 4 parameters")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(fmatrix):
        raise ValueError("'fmatrix' parameter can't be cast to real_matrix")
    __fmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_lsfitreport()
    x_lsfitreport_zero_fields(__rep)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_listlist(__fmatrix, fmatrix, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitlinear(ctypes.byref(_error_msg), ctypes.byref(__y), ctypes.byref(__fmatrix), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitlinear'")
        __r__info = __info.value
        __r__c = list_from_x(__c)
        __r__rep = lsfitreport_from_x(__rep)
        return (__r__info, __r__c, __r__rep)
    finally:
        x_vector_clear(__y)
        x_matrix_clear(__fmatrix)
        x_vector_clear(__c)
        x_lsfitreport_clear(__rep)


_lib_alglib.alglib_lsfitlinearc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitlinearc.restype = ctypes.c_int32
def lsfitlinearc(*functionargs):
    if len(functionargs)==6:
        __friendly_form = False
        y,fmatrix,cmatrix,n,m,k = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        y,fmatrix,cmatrix = functionargs
        if safe_len("'lsfitlinearc': incorrect parameters",y)!=safe_rows("'lsfitlinearc': incorrect parameters",fmatrix):
            raise RuntimeError("Error while calling 'lsfitlinearc': looks like one of arguments has wrong size")
        n = safe_len("'lsfitlinearc': incorrect parameters",y)
        if safe_cols("'lsfitlinearc': incorrect parameters",fmatrix)!=safe_cols("'lsfitlinearc': incorrect parameters",cmatrix)-1:
            raise RuntimeError("Error while calling 'lsfitlinearc': looks like one of arguments has wrong size")
        m = safe_cols("'lsfitlinearc': incorrect parameters",fmatrix)
        k = safe_rows("'lsfitlinearc': incorrect parameters",cmatrix)
    else:
        raise RuntimeError("Error while calling 'lsfitlinearc': function must have 3 or 6 parameters")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(fmatrix):
        raise ValueError("'fmatrix' parameter can't be cast to real_matrix")
    __fmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(cmatrix):
        raise ValueError("'cmatrix' parameter can't be cast to real_matrix")
    __cmatrix = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __info = c_ptrint_t(0)
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_lsfitreport()
    x_lsfitreport_zero_fields(__rep)
    try:
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_listlist(__fmatrix, fmatrix, DT_REAL, X_CREATE)
        x_from_listlist(__cmatrix, cmatrix, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitlinearc(ctypes.byref(_error_msg), ctypes.byref(__y), ctypes.byref(__fmatrix), ctypes.byref(__cmatrix), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitlinearc'")
        __r__info = __info.value
        __r__c = list_from_x(__c)
        __r__rep = lsfitreport_from_x(__rep)
        return (__r__info, __r__c, __r__rep)
    finally:
        x_vector_clear(__y)
        x_matrix_clear(__fmatrix)
        x_matrix_clear(__cmatrix)
        x_vector_clear(__c)
        x_lsfitreport_clear(__rep)


_lib_alglib.alglib_lsfitcreatewf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatewf.restype = ctypes.c_int32
def lsfitcreatewf(*functionargs):
    if len(functionargs)==8:
        __friendly_form = False
        x,y,w,c,n,m,k,diffstep = functionargs
    elif len(functionargs)==5:
        __friendly_form = True
        x,y,w,c,diffstep = functionargs
        if safe_rows("'lsfitcreatewf': incorrect parameters",x)!=safe_len("'lsfitcreatewf': incorrect parameters",y) or safe_rows("'lsfitcreatewf': incorrect parameters",x)!=safe_len("'lsfitcreatewf': incorrect parameters",w):
            raise RuntimeError("Error while calling 'lsfitcreatewf': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatewf': incorrect parameters",x)
        m = safe_cols("'lsfitcreatewf': incorrect parameters",x)
        k = safe_len("'lsfitcreatewf': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatewf': function must have 5 or 8 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __diffstep = ctypes.c_double(diffstep)
    if __diffstep.value!=diffstep:
        raise ValueError("Error while converting 'diffstep' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatewf(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__diffstep), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatewf'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitcreatef.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatef.restype = ctypes.c_int32
def lsfitcreatef(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,c,n,m,k,diffstep = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        x,y,c,diffstep = functionargs
        if safe_rows("'lsfitcreatef': incorrect parameters",x)!=safe_len("'lsfitcreatef': incorrect parameters",y):
            raise RuntimeError("Error while calling 'lsfitcreatef': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatef': incorrect parameters",x)
        m = safe_cols("'lsfitcreatef': incorrect parameters",x)
        k = safe_len("'lsfitcreatef': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatef': function must have 4 or 7 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __diffstep = ctypes.c_double(diffstep)
    if __diffstep.value!=diffstep:
        raise ValueError("Error while converting 'diffstep' parameter to 'ctypes.c_double'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatef(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__diffstep), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatef'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitcreatewfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatewfg.restype = ctypes.c_int32
def lsfitcreatewfg(*functionargs):
    if len(functionargs)==8:
        __friendly_form = False
        x,y,w,c,n,m,k,cheapfg = functionargs
    elif len(functionargs)==5:
        __friendly_form = True
        x,y,w,c,cheapfg = functionargs
        if safe_rows("'lsfitcreatewfg': incorrect parameters",x)!=safe_len("'lsfitcreatewfg': incorrect parameters",y) or safe_rows("'lsfitcreatewfg': incorrect parameters",x)!=safe_len("'lsfitcreatewfg': incorrect parameters",w):
            raise RuntimeError("Error while calling 'lsfitcreatewfg': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatewfg': incorrect parameters",x)
        m = safe_cols("'lsfitcreatewfg': incorrect parameters",x)
        k = safe_len("'lsfitcreatewfg': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatewfg': function must have 5 or 8 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __cheapfg = ctypes.c_uint8(cheapfg)
    if __cheapfg.value!=0:
        __cheapfg = ctypes.c_uint8(1)
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatewfg(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__cheapfg), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatewfg'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitcreatefg.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatefg.restype = ctypes.c_int32
def lsfitcreatefg(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,c,n,m,k,cheapfg = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        x,y,c,cheapfg = functionargs
        if safe_rows("'lsfitcreatefg': incorrect parameters",x)!=safe_len("'lsfitcreatefg': incorrect parameters",y):
            raise RuntimeError("Error while calling 'lsfitcreatefg': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatefg': incorrect parameters",x)
        m = safe_cols("'lsfitcreatefg': incorrect parameters",x)
        k = safe_len("'lsfitcreatefg': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatefg': function must have 4 or 7 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __cheapfg = ctypes.c_uint8(cheapfg)
    if __cheapfg.value!=0:
        __cheapfg = ctypes.c_uint8(1)
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatefg(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__cheapfg), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatefg'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitcreatewfgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatewfgh.restype = ctypes.c_int32
def lsfitcreatewfgh(*functionargs):
    if len(functionargs)==7:
        __friendly_form = False
        x,y,w,c,n,m,k = functionargs
    elif len(functionargs)==4:
        __friendly_form = True
        x,y,w,c = functionargs
        if safe_rows("'lsfitcreatewfgh': incorrect parameters",x)!=safe_len("'lsfitcreatewfgh': incorrect parameters",y) or safe_rows("'lsfitcreatewfgh': incorrect parameters",x)!=safe_len("'lsfitcreatewfgh': incorrect parameters",w):
            raise RuntimeError("Error while calling 'lsfitcreatewfgh': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatewfgh': incorrect parameters",x)
        m = safe_cols("'lsfitcreatewfgh': incorrect parameters",x)
        k = safe_len("'lsfitcreatewfgh': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatewfgh': function must have 4 or 7 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(w):
        raise ValueError("'w' parameter can't be cast to real_vector")
    __w = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__w, w, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatewfgh(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__w), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatewfgh'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__w)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitcreatefgh.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitcreatefgh.restype = ctypes.c_int32
def lsfitcreatefgh(*functionargs):
    if len(functionargs)==6:
        __friendly_form = False
        x,y,c,n,m,k = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,y,c = functionargs
        if safe_rows("'lsfitcreatefgh': incorrect parameters",x)!=safe_len("'lsfitcreatefgh': incorrect parameters",y):
            raise RuntimeError("Error while calling 'lsfitcreatefgh': looks like one of arguments has wrong size")
        n = safe_rows("'lsfitcreatefgh': incorrect parameters",x)
        m = safe_cols("'lsfitcreatefgh': incorrect parameters",x)
        k = safe_len("'lsfitcreatefgh': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'lsfitcreatefgh': function must have 3 or 6 parameters")
    if not is_real_matrix(x):
        raise ValueError("'x' parameter can't be cast to real_matrix")
    __x = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __state = ctypes.c_void_p(0)
    try:
        x_from_listlist(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitcreatefgh(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__k), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitcreatefgh'")
        __r__state = lsfitstate(__state)
        return __r__state
    finally:
        x_matrix_clear(__x)
        x_vector_clear(__y)
        x_vector_clear(__c)


_lib_alglib.alglib_lsfitsetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitsetcond.restype = ctypes.c_int32
def lsfitsetcond(state, epsf, epsx, maxits):
    pass
    __state = state.ptr
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitsetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsf), ctypes.byref(__epsx), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitsetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_lsfitsetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitsetstpmax.restype = ctypes.c_int32
def lsfitsetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitsetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitsetstpmax'")
        return
    finally:
        pass


_lib_alglib.alglib_lsfitsetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitsetxrep.restype = ctypes.c_int32
def lsfitsetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitsetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitsetxrep'")
        return
    finally:
        pass




def lsfitfit_f(state, func, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_c  = x_vector()
    _lib_alglib.x_lsfitstate_get_c(state.ptr, ctypes.byref(_xc_c))
    _py_c = create_real_vector(_xc_c.cnt)
    _xc_x  = x_vector()
    _lib_alglib.x_lsfitstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    while True:
        retval = _lib_alglib.alglib_lsfititeration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfititeration'")
        if not _xc_result:
            break
        _lib_alglib.x_lsfitstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = func(_py_c, _py_x, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_lsfitstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_c, _py_c)
                _lib_alglib.x_lsfitstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_c, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'lsfitfit' (some derivatives were not provided?)")
    return


def lsfitfit_fg(state, func, grad, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_c  = x_vector()
    _lib_alglib.x_lsfitstate_get_c(state.ptr, ctypes.byref(_xc_c))
    _py_c = create_real_vector(_xc_c.cnt)
    _xc_x  = x_vector()
    _lib_alglib.x_lsfitstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_lsfitstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    while True:
        retval = _lib_alglib.alglib_lsfititeration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfititeration'")
        if not _xc_result:
            break
        _lib_alglib.x_lsfitstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = func(_py_c, _py_x, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_lsfitstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = grad(_py_c, _py_x, _py_g, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_lsfitstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_c, _py_c)
                _lib_alglib.x_lsfitstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_c, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'lsfitfit' (some derivatives were not provided?)")
    return


def lsfitfit_fgh(state, func, grad, hess, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_c  = x_vector()
    _lib_alglib.x_lsfitstate_get_c(state.ptr, ctypes.byref(_xc_c))
    _py_c = create_real_vector(_xc_c.cnt)
    _xc_x  = x_vector()
    _lib_alglib.x_lsfitstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_lsfitstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    _xc_h  = x_matrix()
    _lib_alglib.x_lsfitstate_get_h(state.ptr, ctypes.byref(_xc_h))
    _py_h = create_real_matrix(_xc_h.rows,_xc_h.cols)
    while True:
        retval = _lib_alglib.alglib_lsfititeration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfititeration'")
        if not _xc_result:
            break
        _lib_alglib.x_lsfitstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = func(_py_c, _py_x, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_lsfitstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = grad(_py_c, _py_x, _py_g, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_lsfitstate_get_needfgh(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_c, _py_c)
            copy_x_to_list(_xc_x, _py_x)
            _xc_f.value = hess(_py_c, _py_x, _py_g, _py_h, param)
            _lib_alglib.x_lsfitstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_h, _py_h, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_lsfitstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_c, _py_c)
                _lib_alglib.x_lsfitstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_c, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'lsfitfit' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_lsfitresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_lsfitresults.restype = ctypes.c_int32
def lsfitresults(state):
    pass
    __state = state.ptr
    __info = c_ptrint_t(0)
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_lsfitreport()
    x_lsfitreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_lsfitresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__info), ctypes.byref(__c), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'lsfitresults'")
        __r__info = __info.value
        __r__c = list_from_x(__c)
        __r__rep = lsfitreport_from_x(__rep)
        return (__r__info, __r__c, __r__rep)
    finally:
        x_vector_clear(__c)
        x_lsfitreport_clear(__rep)


_lib_alglib.x_obj_free_pspline2interpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_pspline2interpolant.restype = None


class pspline2interpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_pspline2interpolant(self.ptr)
_lib_alglib.x_obj_free_pspline3interpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_pspline3interpolant.restype = None


class pspline3interpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_pspline3interpolant(self.ptr)
_lib_alglib.alglib_pspline2build.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2build.restype = ctypes.c_int32
def pspline2build(xy, n, st, pt):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __st = c_ptrint_t(st)
    if __st.value!=st:
        raise ValueError("Error while converting 'st' parameter to 'c_ptrint_t'")
    __pt = c_ptrint_t(pt)
    if __pt.value!=pt:
        raise ValueError("Error while converting 'pt' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2build(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__st), ctypes.byref(__pt), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2build'")
        __r__p = pspline2interpolant(__p)
        return __r__p
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_pspline3build.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3build.restype = ctypes.c_int32
def pspline3build(xy, n, st, pt):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __st = c_ptrint_t(st)
    if __st.value!=st:
        raise ValueError("Error while converting 'st' parameter to 'c_ptrint_t'")
    __pt = c_ptrint_t(pt)
    if __pt.value!=pt:
        raise ValueError("Error while converting 'pt' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3build(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__st), ctypes.byref(__pt), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3build'")
        __r__p = pspline3interpolant(__p)
        return __r__p
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_pspline2buildperiodic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2buildperiodic.restype = ctypes.c_int32
def pspline2buildperiodic(xy, n, st, pt):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __st = c_ptrint_t(st)
    if __st.value!=st:
        raise ValueError("Error while converting 'st' parameter to 'c_ptrint_t'")
    __pt = c_ptrint_t(pt)
    if __pt.value!=pt:
        raise ValueError("Error while converting 'pt' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2buildperiodic(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__st), ctypes.byref(__pt), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2buildperiodic'")
        __r__p = pspline2interpolant(__p)
        return __r__p
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_pspline3buildperiodic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3buildperiodic.restype = ctypes.c_int32
def pspline3buildperiodic(xy, n, st, pt):
    pass
    if not is_real_matrix(xy):
        raise ValueError("'xy' parameter can't be cast to real_matrix")
    __xy = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __st = c_ptrint_t(st)
    if __st.value!=st:
        raise ValueError("Error while converting 'st' parameter to 'c_ptrint_t'")
    __pt = c_ptrint_t(pt)
    if __pt.value!=pt:
        raise ValueError("Error while converting 'pt' parameter to 'c_ptrint_t'")
    __p = ctypes.c_void_p(0)
    try:
        x_from_listlist(__xy, xy, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3buildperiodic(ctypes.byref(_error_msg), ctypes.byref(__xy), ctypes.byref(__n), ctypes.byref(__st), ctypes.byref(__pt), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3buildperiodic'")
        __r__p = pspline3interpolant(__p)
        return __r__p
    finally:
        x_matrix_clear(__xy)


_lib_alglib.alglib_pspline2parametervalues.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2parametervalues.restype = ctypes.c_int32
def pspline2parametervalues(p):
    pass
    __p = p.ptr
    __n = c_ptrint_t(0)
    __t = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2parametervalues(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2parametervalues'")
        __r__n = __n.value
        __r__t = list_from_x(__t)
        return (__r__n, __r__t)
    finally:
        x_vector_clear(__t)


_lib_alglib.alglib_pspline3parametervalues.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3parametervalues.restype = ctypes.c_int32
def pspline3parametervalues(p):
    pass
    __p = p.ptr
    __n = c_ptrint_t(0)
    __t = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3parametervalues(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__n), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3parametervalues'")
        __r__n = __n.value
        __r__t = list_from_x(__t)
        return (__r__n, __r__t)
    finally:
        x_vector_clear(__t)


_lib_alglib.alglib_pspline2calc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2calc.restype = ctypes.c_int32
def pspline2calc(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2calc(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2calc'")
        __r__x = __x.value
        __r__y = __y.value
        return (__r__x, __r__y)
    finally:
        pass


_lib_alglib.alglib_pspline3calc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3calc.restype = ctypes.c_int32
def pspline3calc(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __z = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3calc(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3calc'")
        __r__x = __x.value
        __r__y = __y.value
        __r__z = __z.value
        return (__r__x, __r__y, __r__z)
    finally:
        pass


_lib_alglib.alglib_pspline2tangent.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2tangent.restype = ctypes.c_int32
def pspline2tangent(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2tangent(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2tangent'")
        __r__x = __x.value
        __r__y = __y.value
        return (__r__x, __r__y)
    finally:
        pass


_lib_alglib.alglib_pspline3tangent.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3tangent.restype = ctypes.c_int32
def pspline3tangent(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __z = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3tangent(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3tangent'")
        __r__x = __x.value
        __r__y = __y.value
        __r__z = __z.value
        return (__r__x, __r__y, __r__z)
    finally:
        pass


_lib_alglib.alglib_pspline2diff.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2diff.restype = ctypes.c_int32
def pspline2diff(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __dx = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __dy = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2diff(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__dx), ctypes.byref(__y), ctypes.byref(__dy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2diff'")
        __r__x = __x.value
        __r__dx = __dx.value
        __r__y = __y.value
        __r__dy = __dy.value
        return (__r__x, __r__dx, __r__y, __r__dy)
    finally:
        pass


_lib_alglib.alglib_pspline3diff.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3diff.restype = ctypes.c_int32
def pspline3diff(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __dx = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __dy = ctypes.c_double(0)
    __z = ctypes.c_double(0)
    __dz = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3diff(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__dx), ctypes.byref(__y), ctypes.byref(__dy), ctypes.byref(__z), ctypes.byref(__dz))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3diff'")
        __r__x = __x.value
        __r__dx = __dx.value
        __r__y = __y.value
        __r__dy = __dy.value
        __r__z = __z.value
        __r__dz = __dz.value
        return (__r__x, __r__dx, __r__y, __r__dy, __r__z, __r__dz)
    finally:
        pass


_lib_alglib.alglib_pspline2diff2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2diff2.restype = ctypes.c_int32
def pspline2diff2(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __dx = ctypes.c_double(0)
    __d2x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __dy = ctypes.c_double(0)
    __d2y = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2diff2(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__dx), ctypes.byref(__d2x), ctypes.byref(__y), ctypes.byref(__dy), ctypes.byref(__d2y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2diff2'")
        __r__x = __x.value
        __r__dx = __dx.value
        __r__d2x = __d2x.value
        __r__y = __y.value
        __r__dy = __dy.value
        __r__d2y = __d2y.value
        return (__r__x, __r__dx, __r__d2x, __r__y, __r__dy, __r__d2y)
    finally:
        pass


_lib_alglib.alglib_pspline3diff2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3diff2.restype = ctypes.c_int32
def pspline3diff2(p, t):
    pass
    __p = p.ptr
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(0)
    __dx = ctypes.c_double(0)
    __d2x = ctypes.c_double(0)
    __y = ctypes.c_double(0)
    __dy = ctypes.c_double(0)
    __d2y = ctypes.c_double(0)
    __z = ctypes.c_double(0)
    __dz = ctypes.c_double(0)
    __d2z = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3diff2(ctypes.byref(_error_msg), ctypes.byref(__p), ctypes.byref(__t), ctypes.byref(__x), ctypes.byref(__dx), ctypes.byref(__d2x), ctypes.byref(__y), ctypes.byref(__dy), ctypes.byref(__d2y), ctypes.byref(__z), ctypes.byref(__dz), ctypes.byref(__d2z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3diff2'")
        __r__x = __x.value
        __r__dx = __dx.value
        __r__d2x = __d2x.value
        __r__y = __y.value
        __r__dy = __dy.value
        __r__d2y = __d2y.value
        __r__z = __z.value
        __r__dz = __dz.value
        __r__d2z = __d2z.value
        return (__r__x, __r__dx, __r__d2x, __r__y, __r__dy, __r__d2y, __r__z, __r__dz, __r__d2z)
    finally:
        pass


_lib_alglib.alglib_pspline2arclength.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline2arclength.restype = ctypes.c_int32
def pspline2arclength(p, a, b):
    pass
    __result = ctypes.c_double(0)
    __p = p.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline2arclength(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__p), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline2arclength'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_pspline3arclength.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pspline3arclength.restype = ctypes.c_int32
def pspline3arclength(p, a, b):
    pass
    __result = ctypes.c_double(0)
    __p = p.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pspline3arclength(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__p), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pspline3arclength'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.x_obj_free_spline2dinterpolant.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_spline2dinterpolant.restype = None


class spline2dinterpolant(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_spline2dinterpolant(self.ptr)
_lib_alglib.alglib_spline2dbuildbilinear.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dbuildbilinear.restype = ctypes.c_int32
def spline2dbuildbilinear(x, y, f, m, n):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(f):
        raise ValueError("'f' parameter can't be cast to real_matrix")
    __f = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_listlist(__f, f, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dbuildbilinear(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__f), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dbuildbilinear'")
        __r__c = spline2dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_matrix_clear(__f)


_lib_alglib.alglib_spline2dbuildbicubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dbuildbicubic.restype = ctypes.c_int32
def spline2dbuildbicubic(x, y, f, m, n):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_matrix(f):
        raise ValueError("'f' parameter can't be cast to real_matrix")
    __f = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        x_from_listlist(__f, f, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dbuildbicubic(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__f), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dbuildbicubic'")
        __r__c = spline2dinterpolant(__c)
        return __r__c
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)
        x_matrix_clear(__f)


_lib_alglib.alglib_spline2dcalc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dcalc.restype = ctypes.c_int32
def spline2dcalc(c, x, y):
    pass
    __result = ctypes.c_double(0)
    __c = c.ptr
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dcalc(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__x), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dcalc'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_spline2ddiff.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2ddiff.restype = ctypes.c_int32
def spline2ddiff(c, x, y):
    pass
    __c = c.ptr
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    __f = ctypes.c_double(0)
    __fx = ctypes.c_double(0)
    __fy = ctypes.c_double(0)
    __fxy = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2ddiff(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__x), ctypes.byref(__y), ctypes.byref(__f), ctypes.byref(__fx), ctypes.byref(__fy), ctypes.byref(__fxy))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2ddiff'")
        __r__f = __f.value
        __r__fx = __fx.value
        __r__fy = __fy.value
        __r__fxy = __fxy.value
        return (__r__f, __r__fx, __r__fy, __r__fxy)
    finally:
        pass


_lib_alglib.alglib_spline2dunpack.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dunpack.restype = ctypes.c_int32
def spline2dunpack(c):
    pass
    __c = c.ptr
    __m = c_ptrint_t(0)
    __n = c_ptrint_t(0)
    __tbl = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dunpack(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__m), ctypes.byref(__n), ctypes.byref(__tbl))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dunpack'")
        __r__m = __m.value
        __r__n = __n.value
        __r__tbl = listlist_from_x(__tbl)
        return (__r__m, __r__n, __r__tbl)
    finally:
        x_matrix_clear(__tbl)


_lib_alglib.alglib_spline2dlintransxy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dlintransxy.restype = ctypes.c_int32
def spline2dlintransxy(c, ax, bx, ay, by):
    pass
    __c = c.ptr
    __ax = ctypes.c_double(ax)
    if __ax.value!=ax:
        raise ValueError("Error while converting 'ax' parameter to 'ctypes.c_double'")
    __bx = ctypes.c_double(bx)
    if __bx.value!=bx:
        raise ValueError("Error while converting 'bx' parameter to 'ctypes.c_double'")
    __ay = ctypes.c_double(ay)
    if __ay.value!=ay:
        raise ValueError("Error while converting 'ay' parameter to 'ctypes.c_double'")
    __by = ctypes.c_double(by)
    if __by.value!=by:
        raise ValueError("Error while converting 'by' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dlintransxy(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__ax), ctypes.byref(__bx), ctypes.byref(__ay), ctypes.byref(__by))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dlintransxy'")
        return
    finally:
        pass


_lib_alglib.alglib_spline2dlintransf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dlintransf.restype = ctypes.c_int32
def spline2dlintransf(c, a, b):
    pass
    __c = c.ptr
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dlintransf(ctypes.byref(_error_msg), ctypes.byref(__c), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dlintransf'")
        return
    finally:
        pass


_lib_alglib.alglib_spline2dresamplebicubic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dresamplebicubic.restype = ctypes.c_int32
def spline2dresamplebicubic(a, oldheight, oldwidth, newheight, newwidth):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __oldheight = c_ptrint_t(oldheight)
    if __oldheight.value!=oldheight:
        raise ValueError("Error while converting 'oldheight' parameter to 'c_ptrint_t'")
    __oldwidth = c_ptrint_t(oldwidth)
    if __oldwidth.value!=oldwidth:
        raise ValueError("Error while converting 'oldwidth' parameter to 'c_ptrint_t'")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __newheight = c_ptrint_t(newheight)
    if __newheight.value!=newheight:
        raise ValueError("Error while converting 'newheight' parameter to 'c_ptrint_t'")
    __newwidth = c_ptrint_t(newwidth)
    if __newwidth.value!=newwidth:
        raise ValueError("Error while converting 'newwidth' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dresamplebicubic(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__oldheight), ctypes.byref(__oldwidth), ctypes.byref(__b), ctypes.byref(__newheight), ctypes.byref(__newwidth))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dresamplebicubic'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_spline2dresamplebilinear.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spline2dresamplebilinear.restype = ctypes.c_int32
def spline2dresamplebilinear(a, oldheight, oldwidth, newheight, newwidth):
    pass
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __oldheight = c_ptrint_t(oldheight)
    if __oldheight.value!=oldheight:
        raise ValueError("Error while converting 'oldheight' parameter to 'c_ptrint_t'")
    __oldwidth = c_ptrint_t(oldwidth)
    if __oldwidth.value!=oldwidth:
        raise ValueError("Error while converting 'oldwidth' parameter to 'c_ptrint_t'")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __newheight = c_ptrint_t(newheight)
    if __newheight.value!=newheight:
        raise ValueError("Error while converting 'newheight' parameter to 'c_ptrint_t'")
    __newwidth = c_ptrint_t(newwidth)
    if __newwidth.value!=newwidth:
        raise ValueError("Error while converting 'newwidth' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spline2dresamplebilinear(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__oldheight), ctypes.byref(__oldwidth), ctypes.byref(__b), ctypes.byref(__newheight), ctypes.byref(__newwidth))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spline2dresamplebilinear'")
        __r__b = listlist_from_x(__b)
        return __r__b
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)


_lib_alglib.alglib_rmatrixludet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixludet.restype = ctypes.c_int32
def rmatrixludet(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,pivots,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,pivots = functionargs
        if safe_rows("'rmatrixludet': incorrect parameters",a)!=safe_cols("'rmatrixludet': incorrect parameters",a) or safe_rows("'rmatrixludet': incorrect parameters",a)!=safe_len("'rmatrixludet': incorrect parameters",pivots):
            raise RuntimeError("Error while calling 'rmatrixludet': looks like one of arguments has wrong size")
        n = safe_rows("'rmatrixludet': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'rmatrixludet': function must have 2 or 3 parameters")
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(pivots):
        raise ValueError("'pivots' parameter can't be cast to int_vector")
    __pivots = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_list(__pivots, pivots, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixludet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__pivots), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixludet'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)


_lib_alglib.alglib_rmatrixdet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixdet.restype = ctypes.c_int32
def rmatrixdet(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_rows("'rmatrixdet': incorrect parameters",a)!=safe_cols("'rmatrixdet': incorrect parameters",a):
            raise RuntimeError("Error while calling 'rmatrixdet': looks like one of arguments has wrong size")
        n = safe_rows("'rmatrixdet': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'rmatrixdet': function must have 1 or 2 parameters")
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixdet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixdet'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_cmatrixludet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixludet.restype = ctypes.c_int32
def cmatrixludet(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,pivots,n = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        a,pivots = functionargs
        if safe_rows("'cmatrixludet': incorrect parameters",a)!=safe_cols("'cmatrixludet': incorrect parameters",a) or safe_rows("'cmatrixludet': incorrect parameters",a)!=safe_len("'cmatrixludet': incorrect parameters",pivots):
            raise RuntimeError("Error while calling 'cmatrixludet': looks like one of arguments has wrong size")
        n = safe_rows("'cmatrixludet': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'cmatrixludet': function must have 2 or 3 parameters")
    __result = x_complex(x=0,y=0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(pivots):
        raise ValueError("'pivots' parameter can't be cast to int_vector")
    __pivots = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        x_from_list(__pivots, pivots, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixludet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__pivots), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixludet'")
        __r__result = complex(__result.x,__result.y)
        return __r__result
    finally:
        x_matrix_clear(__a)
        x_vector_clear(__pivots)


_lib_alglib.alglib_cmatrixdet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_cmatrixdet.restype = ctypes.c_int32
def cmatrixdet(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_rows("'cmatrixdet': incorrect parameters",a)!=safe_cols("'cmatrixdet': incorrect parameters",a):
            raise RuntimeError("Error while calling 'cmatrixdet': looks like one of arguments has wrong size")
        n = safe_rows("'cmatrixdet': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'cmatrixdet': function must have 1 or 2 parameters")
    __result = x_complex(x=0,y=0)
    if not is_complex_matrix(a):
        raise ValueError("'a' parameter can't be cast to complex_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_COMPLEX, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_cmatrixdet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'cmatrixdet'")
        __r__result = complex(__result.x,__result.y)
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_spdmatrixcholeskydet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixcholeskydet.restype = ctypes.c_int32
def spdmatrixcholeskydet(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        a,n = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_rows("'spdmatrixcholeskydet': incorrect parameters",a)!=safe_cols("'spdmatrixcholeskydet': incorrect parameters",a):
            raise RuntimeError("Error while calling 'spdmatrixcholeskydet': looks like one of arguments has wrong size")
        n = safe_rows("'spdmatrixcholeskydet': incorrect parameters",a)
    else:
        raise RuntimeError("Error while calling 'spdmatrixcholeskydet': function must have 1 or 2 parameters")
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixcholeskydet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixcholeskydet'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_spdmatrixdet.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spdmatrixdet.restype = ctypes.c_int32
def spdmatrixdet(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        a,n,isupper = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        a, = functionargs
        if safe_rows("'spdmatrixdet': incorrect parameters",a)!=safe_cols("'spdmatrixdet': incorrect parameters",a):
            raise RuntimeError("Error while calling 'spdmatrixdet': looks like one of arguments has wrong size")
        n = safe_rows("'spdmatrixdet': incorrect parameters",a)
        isupper = False
    else:
        raise RuntimeError("Error while calling 'spdmatrixdet': function must have 1 or 3 parameters")
    __result = ctypes.c_double(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to symmetric")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isupper = ctypes.c_uint8(isupper)
    if __isupper.value!=0:
        __isupper = ctypes.c_uint8(1)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        if __friendly_form:
            if not x_is_symmetric(__a):
                raise ValueError("'a' parameter is not symmetric matrix")
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spdmatrixdet(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isupper))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spdmatrixdet'")
        __r__result = __result.value
        return __r__result
    finally:
        x_matrix_clear(__a)


_lib_alglib.alglib_smatrixgevd.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixgevd.restype = ctypes.c_int32
def smatrixgevd(a, n, isuppera, b, isupperb, zneeded, problemtype):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isuppera = ctypes.c_uint8(isuppera)
    if __isuppera.value!=0:
        __isuppera = ctypes.c_uint8(1)
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __isupperb = ctypes.c_uint8(isupperb)
    if __isupperb.value!=0:
        __isupperb = ctypes.c_uint8(1)
    __zneeded = c_ptrint_t(zneeded)
    if __zneeded.value!=zneeded:
        raise ValueError("Error while converting 'zneeded' parameter to 'c_ptrint_t'")
    __problemtype = c_ptrint_t(problemtype)
    if __problemtype.value!=problemtype:
        raise ValueError("Error while converting 'problemtype' parameter to 'c_ptrint_t'")
    __d = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __z = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixgevd(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isuppera), ctypes.byref(__b), ctypes.byref(__isupperb), ctypes.byref(__zneeded), ctypes.byref(__problemtype), ctypes.byref(__d), ctypes.byref(__z))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixgevd'")
        __r__result = __result.value!=0
        __r__d = list_from_x(__d)
        __r__z = listlist_from_x(__z)
        return (__r__result, __r__d, __r__z)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_vector_clear(__d)
        x_matrix_clear(__z)


_lib_alglib.alglib_smatrixgevdreduce.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_smatrixgevdreduce.restype = ctypes.c_int32
def smatrixgevdreduce(a, n, isuppera, b, isupperb, problemtype):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __isuppera = ctypes.c_uint8(isuppera)
    if __isuppera.value!=0:
        __isuppera = ctypes.c_uint8(1)
    if not is_real_matrix(b):
        raise ValueError("'b' parameter can't be cast to real_matrix")
    __b = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __isupperb = ctypes.c_uint8(isupperb)
    if __isupperb.value!=0:
        __isupperb = ctypes.c_uint8(1)
    __problemtype = c_ptrint_t(problemtype)
    if __problemtype.value!=problemtype:
        raise ValueError("Error while converting 'problemtype' parameter to 'c_ptrint_t'")
    __r = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __isupperr = ctypes.c_uint8(0)
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        x_from_listlist(__b, b, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_smatrixgevdreduce(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__isuppera), ctypes.byref(__b), ctypes.byref(__isupperb), ctypes.byref(__problemtype), ctypes.byref(__r), ctypes.byref(__isupperr))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'smatrixgevdreduce'")
        __r__result = __result.value!=0
        __r__a = listlist_from_x(__a)
        __r__r = listlist_from_x(__r)
        __r__isupperr = __isupperr.value!=0
        return (__r__result, __r__a, __r__r, __r__isupperr)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__b)
        x_matrix_clear(__r)


_lib_alglib.alglib_rmatrixinvupdatesimple.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixinvupdatesimple.restype = ctypes.c_int32
def rmatrixinvupdatesimple(inva, n, updrow, updcolumn, updval):
    pass
    if not is_real_matrix(inva):
        raise ValueError("'inva' parameter can't be cast to real_matrix")
    __inva = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __updrow = c_ptrint_t(updrow)
    if __updrow.value!=updrow:
        raise ValueError("Error while converting 'updrow' parameter to 'c_ptrint_t'")
    __updcolumn = c_ptrint_t(updcolumn)
    if __updcolumn.value!=updcolumn:
        raise ValueError("Error while converting 'updcolumn' parameter to 'c_ptrint_t'")
    __updval = ctypes.c_double(updval)
    if __updval.value!=updval:
        raise ValueError("Error while converting 'updval' parameter to 'ctypes.c_double'")
    try:
        x_from_listlist(__inva, inva, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixinvupdatesimple(ctypes.byref(_error_msg), ctypes.byref(__inva), ctypes.byref(__n), ctypes.byref(__updrow), ctypes.byref(__updcolumn), ctypes.byref(__updval))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixinvupdatesimple'")
        __r__inva = listlist_from_x(__inva)
        return __r__inva
    finally:
        x_matrix_clear(__inva)


_lib_alglib.alglib_rmatrixinvupdaterow.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixinvupdaterow.restype = ctypes.c_int32
def rmatrixinvupdaterow(inva, n, updrow, v):
    pass
    if not is_real_matrix(inva):
        raise ValueError("'inva' parameter can't be cast to real_matrix")
    __inva = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __updrow = c_ptrint_t(updrow)
    if __updrow.value!=updrow:
        raise ValueError("Error while converting 'updrow' parameter to 'c_ptrint_t'")
    if not is_real_vector(v):
        raise ValueError("'v' parameter can't be cast to real_vector")
    __v = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__inva, inva, DT_REAL, X_CREATE)
        x_from_list(__v, v, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixinvupdaterow(ctypes.byref(_error_msg), ctypes.byref(__inva), ctypes.byref(__n), ctypes.byref(__updrow), ctypes.byref(__v))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixinvupdaterow'")
        __r__inva = listlist_from_x(__inva)
        return __r__inva
    finally:
        x_matrix_clear(__inva)
        x_vector_clear(__v)


_lib_alglib.alglib_rmatrixinvupdatecolumn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixinvupdatecolumn.restype = ctypes.c_int32
def rmatrixinvupdatecolumn(inva, n, updcolumn, u):
    pass
    if not is_real_matrix(inva):
        raise ValueError("'inva' parameter can't be cast to real_matrix")
    __inva = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __updcolumn = c_ptrint_t(updcolumn)
    if __updcolumn.value!=updcolumn:
        raise ValueError("Error while converting 'updcolumn' parameter to 'c_ptrint_t'")
    if not is_real_vector(u):
        raise ValueError("'u' parameter can't be cast to real_vector")
    __u = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__inva, inva, DT_REAL, X_CREATE)
        x_from_list(__u, u, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixinvupdatecolumn(ctypes.byref(_error_msg), ctypes.byref(__inva), ctypes.byref(__n), ctypes.byref(__updcolumn), ctypes.byref(__u))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixinvupdatecolumn'")
        __r__inva = listlist_from_x(__inva)
        return __r__inva
    finally:
        x_matrix_clear(__inva)
        x_vector_clear(__u)


_lib_alglib.alglib_rmatrixinvupdateuv.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixinvupdateuv.restype = ctypes.c_int32
def rmatrixinvupdateuv(inva, n, u, v):
    pass
    if not is_real_matrix(inva):
        raise ValueError("'inva' parameter can't be cast to real_matrix")
    __inva = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(u):
        raise ValueError("'u' parameter can't be cast to real_vector")
    __u = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(v):
        raise ValueError("'v' parameter can't be cast to real_vector")
    __v = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__inva, inva, DT_REAL, X_CREATE)
        x_from_list(__u, u, DT_REAL, X_CREATE)
        x_from_list(__v, v, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixinvupdateuv(ctypes.byref(_error_msg), ctypes.byref(__inva), ctypes.byref(__n), ctypes.byref(__u), ctypes.byref(__v))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixinvupdateuv'")
        __r__inva = listlist_from_x(__inva)
        return __r__inva
    finally:
        x_matrix_clear(__inva)
        x_vector_clear(__u)
        x_vector_clear(__v)


_lib_alglib.alglib_rmatrixschur.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_rmatrixschur.restype = ctypes.c_int32
def rmatrixschur(a, n):
    pass
    __result = ctypes.c_uint8(0)
    if not is_real_matrix(a):
        raise ValueError("'a' parameter can't be cast to real_matrix")
    __a = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __s = x_matrix(rows=0,cols=0,stride=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_listlist(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_rmatrixschur(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__s))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'rmatrixschur'")
        __r__result = __result.value!=0
        __r__a = listlist_from_x(__a)
        __r__s = listlist_from_x(__s)
        return (__r__result, __r__a, __r__s)
    finally:
        x_matrix_clear(__a)
        x_matrix_clear(__s)


_lib_alglib.x_obj_free_minasastate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_minasastate.restype = None
_lib_alglib.x_minasastate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_get_needfg.restype = None
_lib_alglib.x_minasastate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_set_needfg.restype = None
_lib_alglib.x_minasastate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_get_xupdated.restype = None
_lib_alglib.x_minasastate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_set_xupdated.restype = None
_lib_alglib.x_minasastate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_get_f.restype = None
_lib_alglib.x_minasastate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_set_f.restype = None
_lib_alglib.x_minasastate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_get_g.restype = None
_lib_alglib.x_minasastate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minasastate_get_x.restype = None


class minasastate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_minasastate(self.ptr)


class x_minasareport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("iterationscount", c_ptrint_t),
        ("nfev", c_ptrint_t),
        ("terminationtype", c_ptrint_t),
        ("activeconstraints", c_ptrint_t)
        ]




class minasareport(object):
    def __init__(self):
        self.iterationscount = 0
        self.nfev = 0
        self.terminationtype = 0
        self.activeconstraints = 0


def x_minasareport_zero_fields(x):
    x.iterationscount = 0
    x.nfev = 0
    x.terminationtype = 0
    x.activeconstraints = 0
    return




def x_minasareport_clear(x):
    x_minasareport_zero_fields(x)
    return




def x_from_minasareport(x,v):
    x.iterationscount = int(v.iterationscount)
    x.nfev = int(v.nfev)
    x.terminationtype = int(v.terminationtype)
    x.activeconstraints = int(v.activeconstraints)
    return




def minasareport_from_x(x):
    r = minasareport()
    r.iterationscount = x.iterationscount
    r.nfev = x.nfev
    r.terminationtype = x.terminationtype
    r.activeconstraints = x.activeconstraints
    return r


_lib_alglib.alglib_minasacreate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasacreate.restype = ctypes.c_int32
def minasacreate(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        n,x,bndl,bndu = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        x,bndl,bndu = functionargs
        if safe_len("'minasacreate': incorrect parameters",x)!=safe_len("'minasacreate': incorrect parameters",bndl) or safe_len("'minasacreate': incorrect parameters",x)!=safe_len("'minasacreate': incorrect parameters",bndu):
            raise RuntimeError("Error while calling 'minasacreate': looks like one of arguments has wrong size")
        n = safe_len("'minasacreate': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minasacreate': function must have 3 or 4 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bndl):
        raise ValueError("'bndl' parameter can't be cast to real_vector")
    __bndl = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bndu):
        raise ValueError("'bndu' parameter can't be cast to real_vector")
    __bndu = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__bndl, bndl, DT_REAL, X_CREATE)
        x_from_list(__bndu, bndu, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasacreate(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__bndl), ctypes.byref(__bndu), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasacreate'")
        __r__state = minasastate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)
        x_vector_clear(__bndl)
        x_vector_clear(__bndu)


_lib_alglib.alglib_minasasetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasasetcond.restype = ctypes.c_int32
def minasasetcond(state, epsg, epsf, epsx, maxits):
    pass
    __state = state.ptr
    __epsg = ctypes.c_double(epsg)
    if __epsg.value!=epsg:
        raise ValueError("Error while converting 'epsg' parameter to 'ctypes.c_double'")
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasasetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsg), ctypes.byref(__epsf), ctypes.byref(__epsx), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasasetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_minasasetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasasetxrep.restype = ctypes.c_int32
def minasasetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasasetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasasetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_minasasetalgorithm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasasetalgorithm.restype = ctypes.c_int32
def minasasetalgorithm(state, algotype):
    pass
    __state = state.ptr
    __algotype = c_ptrint_t(algotype)
    if __algotype.value!=algotype:
        raise ValueError("Error while converting 'algotype' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasasetalgorithm(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__algotype))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasasetalgorithm'")
        return
    finally:
        pass


_lib_alglib.alglib_minasasetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasasetstpmax.restype = ctypes.c_int32
def minasasetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasasetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasasetstpmax'")
        return
    finally:
        pass




def minasaoptimize_g(state, grad, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minasastate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_minasastate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    while True:
        retval = _lib_alglib.alglib_minasaiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasaiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minasastate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_minasastate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minasastate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minasastate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minasaoptimize' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_minasaresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasaresults.restype = ctypes.c_int32
def minasaresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minasareport()
    x_minasareport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasaresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasaresults'")
        __r__x = list_from_x(__x)
        __r__rep = minasareport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minasareport_clear(__rep)


_lib_alglib.alglib_minasaresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasaresultsbuf.restype = ctypes.c_int32
def minasaresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minasareport()
    x_minasareport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_minasareport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasaresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasaresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = minasareport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minasareport_clear(__rep)


_lib_alglib.alglib_minasarestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minasarestartfrom.restype = ctypes.c_int32
def minasarestartfrom(state, x, bndl, bndu):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bndl):
        raise ValueError("'bndl' parameter can't be cast to real_vector")
    __bndl = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bndu):
        raise ValueError("'bndu' parameter can't be cast to real_vector")
    __bndu = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__bndl, bndl, DT_REAL, X_CREATE)
        x_from_list(__bndu, bndu, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minasarestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__bndl), ctypes.byref(__bndu))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minasarestartfrom'")
        return
    finally:
        x_vector_clear(__x)
        x_vector_clear(__bndl)
        x_vector_clear(__bndu)


_lib_alglib.x_obj_free_mincgstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_mincgstate.restype = None
_lib_alglib.x_mincgstate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_get_needfg.restype = None
_lib_alglib.x_mincgstate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_set_needfg.restype = None
_lib_alglib.x_mincgstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_get_xupdated.restype = None
_lib_alglib.x_mincgstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_set_xupdated.restype = None
_lib_alglib.x_mincgstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_get_f.restype = None
_lib_alglib.x_mincgstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_set_f.restype = None
_lib_alglib.x_mincgstate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_get_g.restype = None
_lib_alglib.x_mincgstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_mincgstate_get_x.restype = None


class mincgstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_mincgstate(self.ptr)


class x_mincgreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("iterationscount", c_ptrint_t),
        ("nfev", c_ptrint_t),
        ("terminationtype", c_ptrint_t)
        ]




class mincgreport(object):
    def __init__(self):
        self.iterationscount = 0
        self.nfev = 0
        self.terminationtype = 0


def x_mincgreport_zero_fields(x):
    x.iterationscount = 0
    x.nfev = 0
    x.terminationtype = 0
    return




def x_mincgreport_clear(x):
    x_mincgreport_zero_fields(x)
    return




def x_from_mincgreport(x,v):
    x.iterationscount = int(v.iterationscount)
    x.nfev = int(v.nfev)
    x.terminationtype = int(v.terminationtype)
    return




def mincgreport_from_x(x):
    r = mincgreport()
    r.iterationscount = x.iterationscount
    r.nfev = x.nfev
    r.terminationtype = x.terminationtype
    return r


_lib_alglib.alglib_mincgcreate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgcreate.restype = ctypes.c_int32
def mincgcreate(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        n,x = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'mincgcreate': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'mincgcreate': function must have 1 or 2 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgcreate(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgcreate'")
        __r__state = mincgstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_mincgsetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgsetcond.restype = ctypes.c_int32
def mincgsetcond(state, epsg, epsf, epsx, maxits):
    pass
    __state = state.ptr
    __epsg = ctypes.c_double(epsg)
    if __epsg.value!=epsg:
        raise ValueError("Error while converting 'epsg' parameter to 'ctypes.c_double'")
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgsetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsg), ctypes.byref(__epsf), ctypes.byref(__epsx), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgsetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_mincgsetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgsetxrep.restype = ctypes.c_int32
def mincgsetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgsetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgsetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_mincgsetcgtype.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgsetcgtype.restype = ctypes.c_int32
def mincgsetcgtype(state, cgtype):
    pass
    __state = state.ptr
    __cgtype = c_ptrint_t(cgtype)
    if __cgtype.value!=cgtype:
        raise ValueError("Error while converting 'cgtype' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgsetcgtype(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__cgtype))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgsetcgtype'")
        return
    finally:
        pass


_lib_alglib.alglib_mincgsetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgsetstpmax.restype = ctypes.c_int32
def mincgsetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgsetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgsetstpmax'")
        return
    finally:
        pass


_lib_alglib.alglib_mincgsuggeststep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgsuggeststep.restype = ctypes.c_int32
def mincgsuggeststep(state, stp):
    pass
    __state = state.ptr
    __stp = ctypes.c_double(stp)
    if __stp.value!=stp:
        raise ValueError("Error while converting 'stp' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgsuggeststep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stp))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgsuggeststep'")
        return
    finally:
        pass




def mincgoptimize_g(state, grad, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_mincgstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_mincgstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    while True:
        retval = _lib_alglib.alglib_mincgiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_mincgstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_mincgstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_mincgstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_mincgstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'mincgoptimize' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_mincgresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgresults.restype = ctypes.c_int32
def mincgresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_mincgreport()
    x_mincgreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgresults'")
        __r__x = list_from_x(__x)
        __r__rep = mincgreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_mincgreport_clear(__rep)


_lib_alglib.alglib_mincgresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgresultsbuf.restype = ctypes.c_int32
def mincgresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_mincgreport()
    x_mincgreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_mincgreport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = mincgreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_mincgreport_clear(__rep)


_lib_alglib.alglib_mincgrestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mincgrestartfrom.restype = ctypes.c_int32
def mincgrestartfrom(state, x):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mincgrestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mincgrestartfrom'")
        return
    finally:
        x_vector_clear(__x)


_lib_alglib.x_obj_free_minbleicstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_minbleicstate.restype = None
_lib_alglib.x_minbleicstate_get_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_get_needfg.restype = None
_lib_alglib.x_minbleicstate_set_needfg.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_set_needfg.restype = None
_lib_alglib.x_minbleicstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_get_xupdated.restype = None
_lib_alglib.x_minbleicstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_set_xupdated.restype = None
_lib_alglib.x_minbleicstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_get_f.restype = None
_lib_alglib.x_minbleicstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_set_f.restype = None
_lib_alglib.x_minbleicstate_get_g.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_get_g.restype = None
_lib_alglib.x_minbleicstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_minbleicstate_get_x.restype = None


class minbleicstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_minbleicstate(self.ptr)


class x_minbleicreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("inneriterationscount", c_ptrint_t),
        ("outeriterationscount", c_ptrint_t),
        ("nfev", c_ptrint_t),
        ("terminationtype", c_ptrint_t),
        ("debugeqerr", ctypes.c_double),
        ("debugfs", ctypes.c_double),
        ("debugff", ctypes.c_double),
        ("debugdx", ctypes.c_double)
        ]




class minbleicreport(object):
    def __init__(self):
        self.inneriterationscount = 0
        self.outeriterationscount = 0
        self.nfev = 0
        self.terminationtype = 0
        self.debugeqerr = 0
        self.debugfs = 0
        self.debugff = 0
        self.debugdx = 0


def x_minbleicreport_zero_fields(x):
    x.inneriterationscount = 0
    x.outeriterationscount = 0
    x.nfev = 0
    x.terminationtype = 0
    x.debugeqerr = 0
    x.debugfs = 0
    x.debugff = 0
    x.debugdx = 0
    return




def x_minbleicreport_clear(x):
    x_minbleicreport_zero_fields(x)
    return




def x_from_minbleicreport(x,v):
    x.inneriterationscount = int(v.inneriterationscount)
    x.outeriterationscount = int(v.outeriterationscount)
    x.nfev = int(v.nfev)
    x.terminationtype = int(v.terminationtype)
    x.debugeqerr = float(v.debugeqerr)
    x.debugfs = float(v.debugfs)
    x.debugff = float(v.debugff)
    x.debugdx = float(v.debugdx)
    return




def minbleicreport_from_x(x):
    r = minbleicreport()
    r.inneriterationscount = x.inneriterationscount
    r.outeriterationscount = x.outeriterationscount
    r.nfev = x.nfev
    r.terminationtype = x.terminationtype
    r.debugeqerr = x.debugeqerr
    r.debugfs = x.debugfs
    r.debugff = x.debugff
    r.debugdx = x.debugdx
    return r


_lib_alglib.alglib_minbleiccreate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleiccreate.restype = ctypes.c_int32
def minbleiccreate(*functionargs):
    if len(functionargs)==2:
        __friendly_form = False
        n,x = functionargs
    elif len(functionargs)==1:
        __friendly_form = True
        x, = functionargs
        n = safe_len("'minbleiccreate': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'minbleiccreate': function must have 1 or 2 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleiccreate(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleiccreate'")
        __r__state = minbleicstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_minbleicsetbc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetbc.restype = ctypes.c_int32
def minbleicsetbc(state, bndl, bndu):
    pass
    __state = state.ptr
    if not is_real_vector(bndl):
        raise ValueError("'bndl' parameter can't be cast to real_vector")
    __bndl = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_real_vector(bndu):
        raise ValueError("'bndu' parameter can't be cast to real_vector")
    __bndu = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__bndl, bndl, DT_REAL, X_CREATE)
        x_from_list(__bndu, bndu, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetbc(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__bndl), ctypes.byref(__bndu))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetbc'")
        return
    finally:
        x_vector_clear(__bndl)
        x_vector_clear(__bndu)


_lib_alglib.alglib_minbleicsetlc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetlc.restype = ctypes.c_int32
def minbleicsetlc(*functionargs):
    if len(functionargs)==4:
        __friendly_form = False
        state,c,ct,k = functionargs
    elif len(functionargs)==3:
        __friendly_form = True
        state,c,ct = functionargs
        if safe_rows("'minbleicsetlc': incorrect parameters",c)!=safe_len("'minbleicsetlc': incorrect parameters",ct):
            raise RuntimeError("Error while calling 'minbleicsetlc': looks like one of arguments has wrong size")
        k = safe_rows("'minbleicsetlc': incorrect parameters",c)
    else:
        raise RuntimeError("Error while calling 'minbleicsetlc': function must have 3 or 4 parameters")
    __state = state.ptr
    if not is_real_matrix(c):
        raise ValueError("'c' parameter can't be cast to real_matrix")
    __c = x_matrix(rows=0,cols=0,stride=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    if not is_int_vector(ct):
        raise ValueError("'ct' parameter can't be cast to int_vector")
    __ct = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    try:
        x_from_listlist(__c, c, DT_REAL, X_CREATE)
        x_from_list(__ct, ct, DT_INT, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetlc(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__c), ctypes.byref(__ct), ctypes.byref(__k))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetlc'")
        return
    finally:
        x_matrix_clear(__c)
        x_vector_clear(__ct)


_lib_alglib.alglib_minbleicsetinnercond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetinnercond.restype = ctypes.c_int32
def minbleicsetinnercond(state, epsg, epsf, epsx):
    pass
    __state = state.ptr
    __epsg = ctypes.c_double(epsg)
    if __epsg.value!=epsg:
        raise ValueError("Error while converting 'epsg' parameter to 'ctypes.c_double'")
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetinnercond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsg), ctypes.byref(__epsf), ctypes.byref(__epsx))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetinnercond'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetoutercond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetoutercond.restype = ctypes.c_int32
def minbleicsetoutercond(state, epsx, epsi):
    pass
    __state = state.ptr
    __epsx = ctypes.c_double(epsx)
    if __epsx.value!=epsx:
        raise ValueError("Error while converting 'epsx' parameter to 'ctypes.c_double'")
    __epsi = ctypes.c_double(epsi)
    if __epsi.value!=epsi:
        raise ValueError("Error while converting 'epsi' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetoutercond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsx), ctypes.byref(__epsi))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetoutercond'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetbarrierwidth.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetbarrierwidth.restype = ctypes.c_int32
def minbleicsetbarrierwidth(state, mu):
    pass
    __state = state.ptr
    __mu = ctypes.c_double(mu)
    if __mu.value!=mu:
        raise ValueError("Error while converting 'mu' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetbarrierwidth(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__mu))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetbarrierwidth'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetbarrierdecay.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetbarrierdecay.restype = ctypes.c_int32
def minbleicsetbarrierdecay(state, mudecay):
    pass
    __state = state.ptr
    __mudecay = ctypes.c_double(mudecay)
    if __mudecay.value!=mudecay:
        raise ValueError("Error while converting 'mudecay' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetbarrierdecay(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__mudecay))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetbarrierdecay'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetmaxits.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetmaxits.restype = ctypes.c_int32
def minbleicsetmaxits(state, maxits):
    pass
    __state = state.ptr
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetmaxits(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetmaxits'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetxrep.restype = ctypes.c_int32
def minbleicsetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_minbleicsetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicsetstpmax.restype = ctypes.c_int32
def minbleicsetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicsetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicsetstpmax'")
        return
    finally:
        pass




def minbleicoptimize_g(state, grad, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_minbleicstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_g  = x_vector()
    _lib_alglib.x_minbleicstate_get_g(state.ptr, ctypes.byref(_xc_g))
    _py_g = create_real_vector(_xc_g.cnt)
    while True:
        retval = _lib_alglib.alglib_minbleiciteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleiciteration'")
        if not _xc_result:
            break
        _lib_alglib.x_minbleicstate_get_needfg(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = grad(_py_x, _py_g, param)
            _lib_alglib.x_minbleicstate_set_f(state.ptr, ctypes.byref(_xc_f))
            x_from_list(_xc_g, _py_g, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_minbleicstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_minbleicstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'minbleicoptimize' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_minbleicresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicresults.restype = ctypes.c_int32
def minbleicresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minbleicreport()
    x_minbleicreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicresults'")
        __r__x = list_from_x(__x)
        __r__rep = minbleicreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minbleicreport_clear(__rep)


_lib_alglib.alglib_minbleicresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicresultsbuf.restype = ctypes.c_int32
def minbleicresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_minbleicreport()
    x_minbleicreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_minbleicreport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = minbleicreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_minbleicreport_clear(__rep)


_lib_alglib.alglib_minbleicrestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_minbleicrestartfrom.restype = ctypes.c_int32
def minbleicrestartfrom(state, x):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_minbleicrestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'minbleicrestartfrom'")
        return
    finally:
        x_vector_clear(__x)


_lib_alglib.x_obj_free_nleqstate.argtypes = [ctypes.c_void_p]
_lib_alglib.x_obj_free_nleqstate.restype = None
_lib_alglib.x_nleqstate_get_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_needf.restype = None
_lib_alglib.x_nleqstate_set_needf.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_set_needf.restype = None
_lib_alglib.x_nleqstate_get_needfij.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_needfij.restype = None
_lib_alglib.x_nleqstate_set_needfij.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_set_needfij.restype = None
_lib_alglib.x_nleqstate_get_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_xupdated.restype = None
_lib_alglib.x_nleqstate_set_xupdated.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_set_xupdated.restype = None
_lib_alglib.x_nleqstate_get_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_f.restype = None
_lib_alglib.x_nleqstate_set_f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_set_f.restype = None
_lib_alglib.x_nleqstate_get_fi.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_fi.restype = None
_lib_alglib.x_nleqstate_get_j.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_j.restype = None
_lib_alglib.x_nleqstate_get_x.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.x_nleqstate_get_x.restype = None


class nleqstate(object):
    def __init__(self,ptr):
        self.ptr = ptr
    def __del__(self):
        _lib_alglib.x_obj_free_nleqstate(self.ptr)


class x_nleqreport(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
        ("iterationscount", c_ptrint_t),
        ("nfunc", c_ptrint_t),
        ("njac", c_ptrint_t),
        ("terminationtype", c_ptrint_t)
        ]




class nleqreport(object):
    def __init__(self):
        self.iterationscount = 0
        self.nfunc = 0
        self.njac = 0
        self.terminationtype = 0


def x_nleqreport_zero_fields(x):
    x.iterationscount = 0
    x.nfunc = 0
    x.njac = 0
    x.terminationtype = 0
    return




def x_nleqreport_clear(x):
    x_nleqreport_zero_fields(x)
    return




def x_from_nleqreport(x,v):
    x.iterationscount = int(v.iterationscount)
    x.nfunc = int(v.nfunc)
    x.njac = int(v.njac)
    x.terminationtype = int(v.terminationtype)
    return




def nleqreport_from_x(x):
    r = nleqreport()
    r.iterationscount = x.iterationscount
    r.nfunc = x.nfunc
    r.njac = x.njac
    r.terminationtype = x.terminationtype
    return r


_lib_alglib.alglib_nleqcreatelm.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqcreatelm.restype = ctypes.c_int32
def nleqcreatelm(*functionargs):
    if len(functionargs)==3:
        __friendly_form = False
        n,m,x = functionargs
    elif len(functionargs)==2:
        __friendly_form = True
        m,x = functionargs
        n = safe_len("'nleqcreatelm': incorrect parameters",x)
    else:
        raise RuntimeError("Error while calling 'nleqcreatelm': function must have 2 or 3 parameters")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __state = ctypes.c_void_p(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqcreatelm(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__m), ctypes.byref(__x), ctypes.byref(__state))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqcreatelm'")
        __r__state = nleqstate(__state)
        return __r__state
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_nleqsetcond.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqsetcond.restype = ctypes.c_int32
def nleqsetcond(state, epsf, maxits):
    pass
    __state = state.ptr
    __epsf = ctypes.c_double(epsf)
    if __epsf.value!=epsf:
        raise ValueError("Error while converting 'epsf' parameter to 'ctypes.c_double'")
    __maxits = c_ptrint_t(maxits)
    if __maxits.value!=maxits:
        raise ValueError("Error while converting 'maxits' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqsetcond(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__epsf), ctypes.byref(__maxits))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqsetcond'")
        return
    finally:
        pass


_lib_alglib.alglib_nleqsetxrep.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqsetxrep.restype = ctypes.c_int32
def nleqsetxrep(state, needxrep):
    pass
    __state = state.ptr
    __needxrep = ctypes.c_uint8(needxrep)
    if __needxrep.value!=0:
        __needxrep = ctypes.c_uint8(1)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqsetxrep(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__needxrep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqsetxrep'")
        return
    finally:
        pass


_lib_alglib.alglib_nleqsetstpmax.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqsetstpmax.restype = ctypes.c_int32
def nleqsetstpmax(state, stpmax):
    pass
    __state = state.ptr
    __stpmax = ctypes.c_double(stpmax)
    if __stpmax.value!=stpmax:
        raise ValueError("Error while converting 'stpmax' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqsetstpmax(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__stpmax))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqsetstpmax'")
        return
    finally:
        pass




def nleqsolve_fj(state, func, jac, rep = None, param = None):
    _xc_result = ctypes.c_uint8(0)
    _xc_msg = ctypes.c_char_p()
    _xc_x  = x_vector()
    _lib_alglib.x_nleqstate_get_x(state.ptr, ctypes.byref(_xc_x))
    _py_x = create_real_vector(_xc_x.cnt)
    _xc_flag = ctypes.c_uint8()
    _xc_f = ctypes.c_double()
    _xc_fi  = x_vector()
    _lib_alglib.x_nleqstate_get_fi(state.ptr, ctypes.byref(_xc_fi))
    _py_fi = create_real_vector(_xc_fi.cnt)
    _xc_j  = x_matrix()
    _lib_alglib.x_nleqstate_get_j(state.ptr, ctypes.byref(_xc_j))
    _py_j = create_real_matrix(_xc_j.rows,_xc_j.cols)
    while True:
        retval = _lib_alglib.alglib_nleqiteration(ctypes.byref(_xc_msg), ctypes.byref(_xc_result), ctypes.byref(state.ptr))
        if retval!=0:
            if retval==X_ASSERTION_FAILED:
                raise RuntimeError(_xc_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqiteration'")
        if not _xc_result:
            break
        _lib_alglib.x_nleqstate_get_needf(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            _xc_f.value = func(_py_x, param)
            _lib_alglib.x_nleqstate_set_f(state.ptr, ctypes.byref(_xc_f))
            continue
        _lib_alglib.x_nleqstate_get_needfij(state.ptr, ctypes.byref(_xc_flag))
        if  _xc_flag.value!=0:
            copy_x_to_list(_xc_x, _py_x)

            jac(_py_x, _py_fi, _py_j, param)
            x_from_list(_xc_fi, _py_fi, DT_REAL, X_REWRITE)
            x_from_listlist(_xc_j, _py_j, DT_REAL, X_REWRITE)
            continue
        _lib_alglib.x_nleqstate_get_xupdated(state.ptr, ctypes.byref(_xc_flag))
        if _xc_flag.value!=0 :
            if not (rep is None):
                copy_x_to_list(_xc_x, _py_x)
                _lib_alglib.x_nleqstate_get_f(state.ptr, ctypes.byref(_xc_f))
                rep(_py_x, _xc_f.value, param)
            continue
        raise RuntimeError("ALGLIB: error in 'nleqsolve' (some derivatives were not provided?)")
    return


_lib_alglib.alglib_nleqresults.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqresults.restype = ctypes.c_int32
def nleqresults(state):
    pass
    __state = state.ptr
    __x = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_nleqreport()
    x_nleqreport_zero_fields(__rep)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqresults(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqresults'")
        __r__x = list_from_x(__x)
        __r__rep = nleqreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_nleqreport_clear(__rep)


_lib_alglib.alglib_nleqresultsbuf.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqresultsbuf.restype = ctypes.c_int32
def nleqresultsbuf(state, x, rep):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __rep = x_nleqreport()
    x_nleqreport_zero_fields(__rep)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_nleqreport(__rep, rep)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqresultsbuf(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x), ctypes.byref(__rep))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqresultsbuf'")
        __r__x = list_from_x(__x)
        __r__rep = nleqreport_from_x(__rep)
        return (__r__x, __r__rep)
    finally:
        x_vector_clear(__x)
        x_nleqreport_clear(__rep)


_lib_alglib.alglib_nleqrestartfrom.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_nleqrestartfrom.restype = ctypes.c_int32
def nleqrestartfrom(state, x):
    pass
    __state = state.ptr
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_nleqrestartfrom(ctypes.byref(_error_msg), ctypes.byref(__state), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'nleqrestartfrom'")
        return
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_airy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_airy.restype = ctypes.c_int32
def airy(x):
    pass
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __ai = ctypes.c_double(0)
    __aip = ctypes.c_double(0)
    __bi = ctypes.c_double(0)
    __bip = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_airy(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__ai), ctypes.byref(__aip), ctypes.byref(__bi), ctypes.byref(__bip))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'airy'")
        __r__ai = __ai.value
        __r__aip = __aip.value
        __r__bi = __bi.value
        __r__bip = __bip.value
        return (__r__ai, __r__aip, __r__bi, __r__bip)
    finally:
        pass


_lib_alglib.alglib_besselj0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselj0.restype = ctypes.c_int32
def besselj0(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselj0(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselj0'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besselj1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselj1.restype = ctypes.c_int32
def besselj1(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselj1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselj1'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besseljn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besseljn.restype = ctypes.c_int32
def besseljn(n, x):
    pass
    __result = ctypes.c_double(0)
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besseljn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besseljn'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_bessely0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_bessely0.restype = ctypes.c_int32
def bessely0(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_bessely0(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'bessely0'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_bessely1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_bessely1.restype = ctypes.c_int32
def bessely1(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_bessely1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'bessely1'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besselyn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselyn.restype = ctypes.c_int32
def besselyn(n, x):
    pass
    __result = ctypes.c_double(0)
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselyn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselyn'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besseli0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besseli0.restype = ctypes.c_int32
def besseli0(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besseli0(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besseli0'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besseli1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besseli1.restype = ctypes.c_int32
def besseli1(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besseli1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besseli1'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besselk0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselk0.restype = ctypes.c_int32
def besselk0(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselk0(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselk0'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besselk1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselk1.restype = ctypes.c_int32
def besselk1(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselk1(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselk1'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_besselkn.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_besselkn.restype = ctypes.c_int32
def besselkn(nn, x):
    pass
    __result = ctypes.c_double(0)
    __nn = c_ptrint_t(nn)
    if __nn.value!=nn:
        raise ValueError("Error while converting 'nn' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_besselkn(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__nn), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'besselkn'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_beta.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_beta.restype = ctypes.c_int32
def beta(a, b):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_beta(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'beta'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_incompletebeta.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_incompletebeta.restype = ctypes.c_int32
def incompletebeta(a, b, x):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_incompletebeta(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'incompletebeta'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invincompletebeta.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invincompletebeta.restype = ctypes.c_int32
def invincompletebeta(a, b, y):
    pass
    __result = ctypes.c_double(0)
    __a = ctypes.c_double(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'ctypes.c_double'")
    __b = ctypes.c_double(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'ctypes.c_double'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invincompletebeta(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invincompletebeta'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_binomialdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_binomialdistribution.restype = ctypes.c_int32
def binomialdistribution(k, n, p):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_double(p)
    if __p.value!=p:
        raise ValueError("Error while converting 'p' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_binomialdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'binomialdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_binomialcdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_binomialcdistribution.restype = ctypes.c_int32
def binomialcdistribution(k, n, p):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_double(p)
    if __p.value!=p:
        raise ValueError("Error while converting 'p' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_binomialcdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'binomialcdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invbinomialdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invbinomialdistribution.restype = ctypes.c_int32
def invbinomialdistribution(k, n, y):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invbinomialdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__n), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invbinomialdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_chebyshevcalculate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_chebyshevcalculate.restype = ctypes.c_int32
def chebyshevcalculate(r, n, x):
    pass
    __result = ctypes.c_double(0)
    __r = c_ptrint_t(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_chebyshevcalculate(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'chebyshevcalculate'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_chebyshevsum.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_chebyshevsum.restype = ctypes.c_int32
def chebyshevsum(c, r, n, x):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __r = c_ptrint_t(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'c_ptrint_t'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_chebyshevsum(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'chebyshevsum'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_chebyshevcoefficients.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_chebyshevcoefficients.restype = ctypes.c_int32
def chebyshevcoefficients(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_chebyshevcoefficients(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'chebyshevcoefficients'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_fromchebyshev.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fromchebyshev.restype = ctypes.c_int32
def fromchebyshev(a, n):
    pass
    if not is_real_vector(a):
        raise ValueError("'a' parameter can't be cast to real_vector")
    __a = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __b = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        x_from_list(__a, a, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fromchebyshev(ctypes.byref(_error_msg), ctypes.byref(__a), ctypes.byref(__n), ctypes.byref(__b))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fromchebyshev'")
        __r__b = list_from_x(__b)
        return __r__b
    finally:
        x_vector_clear(__a)
        x_vector_clear(__b)


_lib_alglib.alglib_chisquaredistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_chisquaredistribution.restype = ctypes.c_int32
def chisquaredistribution(v, x):
    pass
    __result = ctypes.c_double(0)
    __v = ctypes.c_double(v)
    if __v.value!=v:
        raise ValueError("Error while converting 'v' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_chisquaredistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__v), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'chisquaredistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_chisquarecdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_chisquarecdistribution.restype = ctypes.c_int32
def chisquarecdistribution(v, x):
    pass
    __result = ctypes.c_double(0)
    __v = ctypes.c_double(v)
    if __v.value!=v:
        raise ValueError("Error while converting 'v' parameter to 'ctypes.c_double'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_chisquarecdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__v), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'chisquarecdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invchisquaredistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invchisquaredistribution.restype = ctypes.c_int32
def invchisquaredistribution(v, y):
    pass
    __result = ctypes.c_double(0)
    __v = ctypes.c_double(v)
    if __v.value!=v:
        raise ValueError("Error while converting 'v' parameter to 'ctypes.c_double'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invchisquaredistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__v), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invchisquaredistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_dawsonintegral.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_dawsonintegral.restype = ctypes.c_int32
def dawsonintegral(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_dawsonintegral(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'dawsonintegral'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_ellipticintegralk.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_ellipticintegralk.restype = ctypes.c_int32
def ellipticintegralk(m):
    pass
    __result = ctypes.c_double(0)
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_ellipticintegralk(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'ellipticintegralk'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_ellipticintegralkhighprecision.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_ellipticintegralkhighprecision.restype = ctypes.c_int32
def ellipticintegralkhighprecision(m1):
    pass
    __result = ctypes.c_double(0)
    __m1 = ctypes.c_double(m1)
    if __m1.value!=m1:
        raise ValueError("Error while converting 'm1' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_ellipticintegralkhighprecision(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__m1))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'ellipticintegralkhighprecision'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_incompleteellipticintegralk.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_incompleteellipticintegralk.restype = ctypes.c_int32
def incompleteellipticintegralk(phi, m):
    pass
    __result = ctypes.c_double(0)
    __phi = ctypes.c_double(phi)
    if __phi.value!=phi:
        raise ValueError("Error while converting 'phi' parameter to 'ctypes.c_double'")
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_incompleteellipticintegralk(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__phi), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'incompleteellipticintegralk'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_ellipticintegrale.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_ellipticintegrale.restype = ctypes.c_int32
def ellipticintegrale(m):
    pass
    __result = ctypes.c_double(0)
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_ellipticintegrale(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'ellipticintegrale'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_incompleteellipticintegrale.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_incompleteellipticintegrale.restype = ctypes.c_int32
def incompleteellipticintegrale(phi, m):
    pass
    __result = ctypes.c_double(0)
    __phi = ctypes.c_double(phi)
    if __phi.value!=phi:
        raise ValueError("Error while converting 'phi' parameter to 'ctypes.c_double'")
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_incompleteellipticintegrale(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__phi), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'incompleteellipticintegrale'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_exponentialintegralei.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_exponentialintegralei.restype = ctypes.c_int32
def exponentialintegralei(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_exponentialintegralei(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'exponentialintegralei'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_exponentialintegralen.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_exponentialintegralen.restype = ctypes.c_int32
def exponentialintegralen(x, n):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_exponentialintegralen(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x), ctypes.byref(__n))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'exponentialintegralen'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_fdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fdistribution.restype = ctypes.c_int32
def fdistribution(a, b, x):
    pass
    __result = ctypes.c_double(0)
    __a = c_ptrint_t(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'c_ptrint_t'")
    __b = c_ptrint_t(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_fcdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fcdistribution.restype = ctypes.c_int32
def fcdistribution(a, b, x):
    pass
    __result = ctypes.c_double(0)
    __a = c_ptrint_t(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'c_ptrint_t'")
    __b = c_ptrint_t(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fcdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fcdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invfdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invfdistribution.restype = ctypes.c_int32
def invfdistribution(a, b, y):
    pass
    __result = ctypes.c_double(0)
    __a = c_ptrint_t(a)
    if __a.value!=a:
        raise ValueError("Error while converting 'a' parameter to 'c_ptrint_t'")
    __b = c_ptrint_t(b)
    if __b.value!=b:
        raise ValueError("Error while converting 'b' parameter to 'c_ptrint_t'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invfdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__a), ctypes.byref(__b), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invfdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_fresnelintegral.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_fresnelintegral.restype = ctypes.c_int32
def fresnelintegral(x, c, s):
    pass
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __c = ctypes.c_double(c)
    if __c.value!=c:
        raise ValueError("Error while converting 'c' parameter to 'ctypes.c_double'")
    __s = ctypes.c_double(s)
    if __s.value!=s:
        raise ValueError("Error while converting 's' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_fresnelintegral(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__c), ctypes.byref(__s))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'fresnelintegral'")
        __r__c = __c.value
        __r__s = __s.value
        return (__r__c, __r__s)
    finally:
        pass


_lib_alglib.alglib_hermitecalculate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hermitecalculate.restype = ctypes.c_int32
def hermitecalculate(n, x):
    pass
    __result = ctypes.c_double(0)
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hermitecalculate(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hermitecalculate'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_hermitesum.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hermitesum.restype = ctypes.c_int32
def hermitesum(c, n, x):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hermitesum(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hermitesum'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_hermitecoefficients.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hermitecoefficients.restype = ctypes.c_int32
def hermitecoefficients(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hermitecoefficients(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hermitecoefficients'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_jacobianellipticfunctions.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_jacobianellipticfunctions.restype = ctypes.c_int32
def jacobianellipticfunctions(u, m):
    pass
    __u = ctypes.c_double(u)
    if __u.value!=u:
        raise ValueError("Error while converting 'u' parameter to 'ctypes.c_double'")
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    __sn = ctypes.c_double(0)
    __cn = ctypes.c_double(0)
    __dn = ctypes.c_double(0)
    __ph = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_jacobianellipticfunctions(ctypes.byref(_error_msg), ctypes.byref(__u), ctypes.byref(__m), ctypes.byref(__sn), ctypes.byref(__cn), ctypes.byref(__dn), ctypes.byref(__ph))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'jacobianellipticfunctions'")
        __r__sn = __sn.value
        __r__cn = __cn.value
        __r__dn = __dn.value
        __r__ph = __ph.value
        return (__r__sn, __r__cn, __r__dn, __r__ph)
    finally:
        pass


_lib_alglib.alglib_laguerrecalculate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_laguerrecalculate.restype = ctypes.c_int32
def laguerrecalculate(n, x):
    pass
    __result = ctypes.c_double(0)
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_laguerrecalculate(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'laguerrecalculate'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_laguerresum.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_laguerresum.restype = ctypes.c_int32
def laguerresum(c, n, x):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_laguerresum(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'laguerresum'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_laguerrecoefficients.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_laguerrecoefficients.restype = ctypes.c_int32
def laguerrecoefficients(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_laguerrecoefficients(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'laguerrecoefficients'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_legendrecalculate.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_legendrecalculate.restype = ctypes.c_int32
def legendrecalculate(n, x):
    pass
    __result = ctypes.c_double(0)
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_legendrecalculate(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'legendrecalculate'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_legendresum.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_legendresum.restype = ctypes.c_int32
def legendresum(c, n, x):
    pass
    __result = ctypes.c_double(0)
    if not is_real_vector(c):
        raise ValueError("'c' parameter can't be cast to real_vector")
    __c = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        x_from_list(__c, c, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_legendresum(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__c), ctypes.byref(__n), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'legendresum'")
        __r__result = __result.value
        return __r__result
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_legendrecoefficients.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_legendrecoefficients.restype = ctypes.c_int32
def legendrecoefficients(n):
    pass
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __c = x_vector(cnt=0,datatype=DT_REAL,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_legendrecoefficients(ctypes.byref(_error_msg), ctypes.byref(__n), ctypes.byref(__c))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'legendrecoefficients'")
        __r__c = list_from_x(__c)
        return __r__c
    finally:
        x_vector_clear(__c)


_lib_alglib.alglib_poissondistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_poissondistribution.restype = ctypes.c_int32
def poissondistribution(k, m):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_poissondistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'poissondistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_poissoncdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_poissoncdistribution.restype = ctypes.c_int32
def poissoncdistribution(k, m):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __m = ctypes.c_double(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_poissoncdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__m))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'poissoncdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invpoissondistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invpoissondistribution.restype = ctypes.c_int32
def invpoissondistribution(k, y):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __y = ctypes.c_double(y)
    if __y.value!=y:
        raise ValueError("Error while converting 'y' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invpoissondistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__y))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invpoissondistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_psi.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_psi.restype = ctypes.c_int32
def psi(x):
    pass
    __result = ctypes.c_double(0)
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_psi(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__x))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'psi'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_studenttdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_studenttdistribution.restype = ctypes.c_int32
def studenttdistribution(k, t):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __t = ctypes.c_double(t)
    if __t.value!=t:
        raise ValueError("Error while converting 't' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_studenttdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__t))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'studenttdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_invstudenttdistribution.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_invstudenttdistribution.restype = ctypes.c_int32
def invstudenttdistribution(k, p):
    pass
    __result = ctypes.c_double(0)
    __k = c_ptrint_t(k)
    if __k.value!=k:
        raise ValueError("Error while converting 'k' parameter to 'c_ptrint_t'")
    __p = ctypes.c_double(p)
    if __p.value!=p:
        raise ValueError("Error while converting 'p' parameter to 'ctypes.c_double'")
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_invstudenttdistribution(ctypes.byref(_error_msg), ctypes.byref(__result), ctypes.byref(__k), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'invstudenttdistribution'")
        __r__result = __result.value
        return __r__result
    finally:
        pass


_lib_alglib.alglib_sinecosineintegrals.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_sinecosineintegrals.restype = ctypes.c_int32
def sinecosineintegrals(x):
    pass
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __si = ctypes.c_double(0)
    __ci = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_sinecosineintegrals(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__si), ctypes.byref(__ci))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'sinecosineintegrals'")
        __r__si = __si.value
        __r__ci = __ci.value
        return (__r__si, __r__ci)
    finally:
        pass


_lib_alglib.alglib_hyperbolicsinecosineintegrals.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_hyperbolicsinecosineintegrals.restype = ctypes.c_int32
def hyperbolicsinecosineintegrals(x):
    pass
    __x = ctypes.c_double(x)
    if __x.value!=x:
        raise ValueError("Error while converting 'x' parameter to 'ctypes.c_double'")
    __shi = ctypes.c_double(0)
    __chi = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_hyperbolicsinecosineintegrals(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__shi), ctypes.byref(__chi))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'hyperbolicsinecosineintegrals'")
        __r__shi = __shi.value
        __r__chi = __chi.value
        return (__r__shi, __r__chi)
    finally:
        pass


_lib_alglib.alglib_pearsoncorrelationsignificance.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_pearsoncorrelationsignificance.restype = ctypes.c_int32
def pearsoncorrelationsignificance(r, n):
    pass
    __r = ctypes.c_double(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_pearsoncorrelationsignificance(ctypes.byref(_error_msg), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'pearsoncorrelationsignificance'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        pass


_lib_alglib.alglib_spearmanrankcorrelationsignificance.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_spearmanrankcorrelationsignificance.restype = ctypes.c_int32
def spearmanrankcorrelationsignificance(r, n):
    pass
    __r = ctypes.c_double(r)
    if __r.value!=r:
        raise ValueError("Error while converting 'r' parameter to 'ctypes.c_double'")
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        pass
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_spearmanrankcorrelationsignificance(ctypes.byref(_error_msg), ctypes.byref(__r), ctypes.byref(__n), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'spearmanrankcorrelationsignificance'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        pass


_lib_alglib.alglib_jarqueberatest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_jarqueberatest.restype = ctypes.c_int32
def jarqueberatest(x, n):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __p = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_jarqueberatest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__p))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'jarqueberatest'")
        __r__p = __p.value
        return __r__p
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_mannwhitneyutest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_mannwhitneyutest.restype = ctypes.c_int32
def mannwhitneyutest(x, n, y, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_mannwhitneyutest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__y), ctypes.byref(__m), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'mannwhitneyutest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_onesamplesigntest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_onesamplesigntest.restype = ctypes.c_int32
def onesamplesigntest(x, n, median):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __median = ctypes.c_double(median)
    if __median.value!=median:
        raise ValueError("Error while converting 'median' parameter to 'ctypes.c_double'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_onesamplesigntest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__median), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'onesamplesigntest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_studentttest1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_studentttest1.restype = ctypes.c_int32
def studentttest1(x, n, mean):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __mean = ctypes.c_double(mean)
    if __mean.value!=mean:
        raise ValueError("Error while converting 'mean' parameter to 'ctypes.c_double'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_studentttest1(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__mean), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'studentttest1'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_studentttest2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_studentttest2.restype = ctypes.c_int32
def studentttest2(x, n, y, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_studentttest2(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__y), ctypes.byref(__m), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'studentttest2'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_unequalvariancettest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_unequalvariancettest.restype = ctypes.c_int32
def unequalvariancettest(x, n, y, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_unequalvariancettest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__y), ctypes.byref(__m), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'unequalvariancettest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_ftest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_ftest.restype = ctypes.c_int32
def ftest(x, n, y, m):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    if not is_real_vector(y):
        raise ValueError("'y' parameter can't be cast to real_vector")
    __y = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __m = c_ptrint_t(m)
    if __m.value!=m:
        raise ValueError("Error while converting 'm' parameter to 'c_ptrint_t'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        x_from_list(__y, y, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_ftest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__y), ctypes.byref(__m), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'ftest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)
        x_vector_clear(__y)


_lib_alglib.alglib_onesamplevariancetest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_onesamplevariancetest.restype = ctypes.c_int32
def onesamplevariancetest(x, n, variance):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __variance = ctypes.c_double(variance)
    if __variance.value!=variance:
        raise ValueError("Error while converting 'variance' parameter to 'ctypes.c_double'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_onesamplevariancetest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__variance), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'onesamplevariancetest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)


_lib_alglib.alglib_wilcoxonsignedranktest.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib_alglib.alglib_wilcoxonsignedranktest.restype = ctypes.c_int32
def wilcoxonsignedranktest(x, n, e):
    pass
    if not is_real_vector(x):
        raise ValueError("'x' parameter can't be cast to real_vector")
    __x = x_vector(cnt=0,datatype=0,owner=OWN_CALLER,last_action=0,ptr=x_multiptr(p_ptr=0))
    __n = c_ptrint_t(n)
    if __n.value!=n:
        raise ValueError("Error while converting 'n' parameter to 'c_ptrint_t'")
    __e = ctypes.c_double(e)
    if __e.value!=e:
        raise ValueError("Error while converting 'e' parameter to 'ctypes.c_double'")
    __bothtails = ctypes.c_double(0)
    __lefttail = ctypes.c_double(0)
    __righttail = ctypes.c_double(0)
    try:
        x_from_list(__x, x, DT_REAL, X_CREATE)
        _error_msg = ctypes.c_char_p(0)
        __x__retval =  _lib_alglib.alglib_wilcoxonsignedranktest(ctypes.byref(_error_msg), ctypes.byref(__x), ctypes.byref(__n), ctypes.byref(__e), ctypes.byref(__bothtails), ctypes.byref(__lefttail), ctypes.byref(__righttail))
        if __x__retval!=0:
            if __x__retval==X_ASSERTION_FAILED:
                raise RuntimeError(_error_msg.value)
            else:
                raise RuntimeError("Error while calling 'wilcoxonsignedranktest'")
        __r__bothtails = __bothtails.value
        __r__lefttail = __lefttail.value
        __r__righttail = __righttail.value
        return (__r__bothtails, __r__lefttail, __r__righttail)
    finally:
        x_vector_clear(__x)


