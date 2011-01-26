import xalglib
import math
import cmath
import sys
import random

def doc_print_test(v, t, typecode, threshold = 0):
    if typecode=="bool":
        if (not xalglib.is_bool(v)) or (not xalglib.is_bool(t)):
            return False
        return v==t
    if typecode=="int":
        if (not xalglib.is_int(v)) or (not xalglib.is_int(t)):
            return False
        return v==t
    if typecode=="real":
        if (not xalglib.is_real(v)) or (not xalglib.is_real(t)):
            return False
        return abs(v-t)<=threshold    
    if typecode=="complex":
        if (not xalglib.is_complex(v)) or (not xalglib.is_complex(t)):
            return False
        return abs(v-t)<=threshold
    
    if typecode=="bool_vector":
        if (not xalglib.is_bool_vector(t)) or (not xalglib.is_bool_vector(v)):
            return False
        if xalglib.safe_len("",v)!=xalglib.safe_len("",t):
            return False
        cnt = xalglib.safe_len("",t)
        for idx in range(cnt):
            if v[idx]!=t[idx]:
                return False
        return True
    
    if typecode=="bool_matrix":
        if (not xalglib.is_bool_matrix(t)) or (not xalglib.is_bool_matrix(v)):
            return False
        if xalglib.safe_rows("",v)!=xalglib.safe_rows("",t):
            return False
        if xalglib.safe_cols("",v)!=xalglib.safe_cols("",t):
            return False
        cols = xalglib.safe_cols("",t)
        rows = xalglib.safe_rows("",t)
        for idx0 in range(rows):
            for idx1 in range(cols):
                if v[idx0][idx1]!=t[idx0][idx1]:
                    return False
        return True
    
    if typecode=="int_vector":
        if (not xalglib.is_int_vector(t)) or (not xalglib.is_int_vector(v)):
            return False
        if xalglib.safe_len("",v)!=xalglib.safe_len("",t):
            return False
        cnt = xalglib.safe_len("",t)
        for idx in range(cnt):
            if v[idx]!=t[idx]:
                return False
        return True
    
    if typecode=="int_matrix":
        if (not xalglib.is_int_matrix(t)) or (not xalglib.is_int_matrix(v)):
            return False
        if xalglib.safe_rows("",v)!=xalglib.safe_rows("",t):
            return False
        if xalglib.safe_cols("",v)!=xalglib.safe_cols("",t):
            return False
        cols = xalglib.safe_cols("",t)
        rows = xalglib.safe_rows("",t)
        for idx0 in range(rows):
            for idx1 in range(cols):
                if v[idx0][idx1]!=t[idx0][idx1]:
                    return False
        return True
            
    if typecode=="real_vector":
        if (not xalglib.is_real_vector(t)) or (not xalglib.is_real_vector(v)):
            return False
        if xalglib.safe_len("",v)!=xalglib.safe_len("",t):
            return False
        cnt = xalglib.safe_len("",t)
        for idx in range(cnt):
            if abs(v[idx]-t[idx])>threshold:
                return False
        return True
    
    if typecode=="real_matrix":
        if (not xalglib.is_real_matrix(t)) or (not xalglib.is_real_matrix(v)):
            return False
        if xalglib.safe_rows("",v)!=xalglib.safe_rows("",t):
            return False
        if xalglib.safe_cols("",v)!=xalglib.safe_cols("",t):
            return False
        cols = xalglib.safe_cols("",t)
        rows = xalglib.safe_rows("",t)
        for idx0 in range(rows):
            for idx1 in range(cols):
                if abs(v[idx0][idx1]-t[idx0][idx1])>threshold:
                    return False
        return True
            
    if typecode=="complex_vector":
        if (not xalglib.is_complex_vector(t)) or (not xalglib.is_complex_vector(v)):
            return False
        if xalglib.safe_len("",v)!=xalglib.safe_len("",t):
            return False
        cnt = xalglib.safe_len("",t)
        for idx in range(cnt):
            if abs(v[idx]-t[idx])>threshold:
                return False
        return True
    
    if typecode=="complex_matrix":
        if (not xalglib.is_complex_matrix(t)) or (not xalglib.is_complex_matrix(v)):
            return False
        if xalglib.safe_rows("",v)!=xalglib.safe_rows("",t):
            return False
        if xalglib.safe_cols("",v)!=xalglib.safe_cols("",t):
            return False
        cols = xalglib.safe_cols("",t)
        rows = xalglib.safe_rows("",t)
        for idx0 in range(rows):
            for idx1 in range(cols):
                if abs(v[idx0][idx1]-t[idx0][idx1])>threshold:
                    return False
        return True
    raise ValueError("doc_print_test: unknown test value type")

def spoil_vec_by_adding_element(x):
    x.append(0)
    return

def spoil_vec_by_deleting_element(x):
    if len(x)!=0:
        x.pop()
    return

def spoil_mat_by_adding_row(x):
    cols = xalglib.safe_cols("internal error",x)
    if cols!=0:
        x.append([0]*cols)
    return

def spoil_mat_by_deleting_row(x):
    if len(x)!=0:
        x.pop()
    return

def spoil_mat_by_adding_col(x):
    for row in x:
        row.append(0)
    return

def spoil_mat_by_deleting_col(x):
    for row in x:
        row.pop()
    return

def spoil_vec_by_nan(x):
    lenx = xalglib.safe_len("internal error",x)
    if lenx!=0:
        i = random.randint(0,lenx-1)
        x[i] = float("nan")
    return

def spoil_vec_by_posinf(x):
    lenx = xalglib.safe_len("internal error",x)
    if lenx!=0:
        i = random.randint(0,lenx-1)
        x[i] = float("+inf")
    return

def spoil_vec_by_neginf(x):
    lenx = xalglib.safe_len("internal error",x)
    if lenx!=0:
        i = random.randint(0,lenx-1)
        x[i] = float("-inf")
    return

def spoil_mat_by_nan(x):
    rows = xalglib.safe_rows("internal error",x)
    cols = xalglib.safe_cols("internal error",x)
    if cols*rows!=0:
        x[random.randint(0,rows-1)][random.randint(0,cols-1)] = float("nan")
    return

def spoil_mat_by_posinf(x):
    rows = xalglib.safe_rows("internal error",x)
    cols = xalglib.safe_cols("internal error",x)
    if cols*rows!=0:
        x[random.randint(0,rows-1)][random.randint(0,cols-1)] = float("+inf")
    return

def spoil_mat_by_neginf(x):
    rows = xalglib.safe_rows("internal error",x)
    cols = xalglib.safe_cols("internal error",x)
    if cols*rows!=0:
        x[random.randint(0,rows-1)][random.randint(0,cols-1)] = float("-inf")
    return

def function1_func(x, param):
    #
    # this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    #
    return 100*(x[0]+3)**4 + (x[1]-3)**4

def function1_grad(x, grad, param):
    #
    # this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    # and its derivatives df/d0 and df/dx1
    #
    func = 100*(x[0]+3)**4 + (x[1]-3)**4
    grad[0] = 400*(x[0]+3)**3
    grad[1] = 4*(x[1]-3)**3
    return func

def function1_hess(x, grad, hess, param):
    #
    # this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    # its derivatives df/d0 and df/dx1
    # and its Hessian.
    #
    func = 100*(x[0]+3)**4 + (x[1]-3)**4
    grad[0] = 400*(x[0]+3)**3
    grad[1] = 4*(x[1]-3)**3
    hess[0][0] = 1200*(x[0]+3)**2
    hess[0][1] = 0
    hess[1][0] = 0
    hess[1][1] = 12*(x[1]-3)**2
    return func

def function1_fvec(x, fi, param):
    #
    # this callback calculates
    # f0(x0,x1) = 100*(x0+3)^4,
    # f1(x0,x1) = (x1-3)^4
    #
    fi[0] = 10*(x[0]+3)**2
    fi[1] = (x[1]-3)**2
    return

def function1_jac(x, fi, jac, param):
    #
    # this callback calculates
    # f0(x0,x1) = 100*(x0+3)^4,
    # f1(x0,x1) = (x1-3)^4
    # and Jacobian matrix J = [dfi/dxj]
    #
    fi[0] = 10*(x[0]+3)**2
    fi[1] = (x[1]-3)**2
    jac[0][0] = 20*(x[0]+3)
    jac[0][1] = 0
    jac[1][0] = 0
    jac[1][1] = 2*(x[1]-3)
    return

def function2_func(x, param):
    #
    # this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
    #
    return (x[0]*x[0]+1)**2 + (x[1]-1)**2

def function2_grad(x, grad, param):
    #
    # this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
    # and its derivatives df/d0 and df/dx1
    #
    func = (x[0]*x[0]+1)**2 + (x[1]-1)**2
    grad[0] = 4*(x[0]*x[0]+1)*x[0]
    grad[1] = 2*(x[1]-1)
    return func

def function2_hess(x, grad, hess, param):
    #
    # this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
    # its gradient and Hessian
    #
    func = (x[0]*x[0]+1)**2 + (x[1]-1)**2
    grad[0] = 4*(x[0]*x[0]+1)*x[0]
    grad[1] = 2*(x[1]-1)
    hess[0][0] = 12*x[0]*x[0]+4
    hess[0][1] = 0
    hess[1][0] = 0
    hess[1][1] = 2
    return func

def function2_fvec(x, fi, param):
    #
    # this callback calculates
    # f0(x0,x1) = x0^2+1
    # f1(x0,x1) = x1-1
    #
    fi[0] = x[0]*x[0]+1
    fi[1] = x[1]-1
    return

def function2_jac(x, fi, jac, param):
    #
    # this callback calculates
    # f0(x0,x1) = x0^2+1
    # f1(x0,x1) = x1-1
    # and Jacobian matrix J = [dfi/dxj]
    #
    fi[0] = x[0]*x[0]+1
    fi[1] = x[1]-1
    jac[0][0] = 2*x[0]
    jac[0][1] = 0
    jac[1][0] = 0
    jac[1][1] = 1
    return

def function_cx_1_func(c, x, param):
    #
    # this callback calculates f(c,x)=exp(-c0*sqr(x0))
    # where x is a position on X-axis and c is adjustable parameter
    #
    return math.exp(-c[0]*(x[0])**2)

def function_cx_1_grad(c, x, grad, param):
    #
    # this callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient G={df/dc[i]}
    # where x is a position on X-axis and c is adjustable parameter.
    # IMPORTANT: gradient is calculated with respect to C, not to X
    #
    func = math.exp(-c[0]*(x[0])**2)
    grad[0] = -((x[0])**2)*func
    return func

def function_cx_1_hess(c, x, grad, hess, param):
    #
    # this callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])}
    # where x is a position on X-axis and c is adjustable parameter.
    # IMPORTANT: gradient/Hessian are calculated with respect to C, not to X
    #
    func = math.exp(-c[0]*(x[0])**2)
    grad[0] = -((x[0])**2)*func
    hess[0][0] = ((x[0])**4)*func
    return func

def ode_function_1_diff(y, x, dy, param):
    #
    # this callback calculates f(y[],x)=-y[0]
    #
    dy[0] = -y[0]
    return

def int_function_1_func(x, xminusa, bminusx, param):
    #
    # this callback calculates f(x)=exp(x)
    #
    return math.exp(x)

_TotalResult = True
sys.stdout.write("Python-ALGLIB communication tests. Please wait...\n")
try:
    #
    # TEST nneighbor_d_1
    #      Nearest neighbor search, KNN queries
    #
    sys.stdout.write("0/78\n")
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            a = [[0,0],[0,1],[1,0],[1,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            nx = 2
            ny = 0
            normtype = 2
            r = [[]]
            kdt = xalglib.kdtreebuild(a, nx, ny, normtype)
            x = [-1,0]
            k = xalglib.kdtreequeryknn(kdt, x, 1)
            _TestResult = _TestResult and doc_print_test(k, 1, "int")
            r = xalglib.kdtreequeryresultsx(kdt, r)
            _TestResult = _TestResult and doc_print_test(r, [[0,0]], "real_matrix", 0.05)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("nneighbor_d_1                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST nneighbor_t_2
    #      Subsequent queries; buffered functions must use previously allocated storage (if large enough), so buffer may contain some info from previous call
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            a = [[0,0],[0,1],[1,0],[1,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            nx = 2
            ny = 0
            normtype = 2
            rx = [[]]
            kdt = xalglib.kdtreebuild(a, nx, ny, normtype)
            x = [+2,0]
            k = xalglib.kdtreequeryknn(kdt, x, 2, True)
            _TestResult = _TestResult and doc_print_test(k, 2, "int")
            rx = xalglib.kdtreequeryresultsx(kdt, rx)
            _TestResult = _TestResult and doc_print_test(rx, [[1,0],[1,1]], "real_matrix", 0.05)
            x = [-2,0]
            k = xalglib.kdtreequeryknn(kdt, x, 1, True)
            _TestResult = _TestResult and doc_print_test(k, 1, "int")
            rx = xalglib.kdtreequeryresultsx(kdt, rx)
            _TestResult = _TestResult and doc_print_test(rx, [[0,0],[1,1]], "real_matrix", 0.05)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("nneighbor_t_2                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_d_base
    #      Basic functionality (moments, adev, median, percentile)
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            x = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)

            #
            # Here we demonstrate calculation of sample moments
            # (mean, variance, skewness, kurtosis)
            #
            mean, variance, skewness, kurtosis = xalglib.samplemoments(x)
            _TestResult = _TestResult and doc_print_test(mean, 28.5, "real", 0.01)
            _TestResult = _TestResult and doc_print_test(variance, 801.1667, "real", 0.01)
            _TestResult = _TestResult and doc_print_test(skewness, 0.5751, "real", 0.01)
            _TestResult = _TestResult and doc_print_test(kurtosis, -1.2666, "real", 0.01)

            #
            # Average deviation
            #
            adev = xalglib.sampleadev(x)
            _TestResult = _TestResult and doc_print_test(adev, 23.2, "real", 0.01)

            #
            # Median and percentile
            #
            v = xalglib.samplemedian(x)
            _TestResult = _TestResult and doc_print_test(v, 20.5, "real", 0.01)
            p = 0.5
            if _spoil_scenario==3:
                p = float("nan")
            if _spoil_scenario==4:
                p = float("+inf")
            if _spoil_scenario==5:
                p = float("-inf")
            v = xalglib.samplepercentile(x, p)
            _TestResult = _TestResult and doc_print_test(v, 20.5, "real", 0.01)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_d_base                  FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_d_c2
    #      Correlation (covariance) between two random variables
    #
    _TestResult = True
    for _spoil_scenario in range(-1,10):
        try:
            #
            # We have two samples - x and y, and want to measure dependency between them
            #
            x = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)

            #
            # Three dependency measures are calculated:
            # * covariation
            # * Pearson correlation
            # * Spearman rank correlation
            #
            v = xalglib.cov2(x, y)
            _TestResult = _TestResult and doc_print_test(v, 82.5, "real", 0.001)
            v = xalglib.pearsoncorr2(x, y)
            _TestResult = _TestResult and doc_print_test(v, 0.9627, "real", 0.001)
            v = xalglib.spearmancorr2(x, y)
            _TestResult = _TestResult and doc_print_test(v, 1.000, "real", 0.001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_d_c2                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_d_cm
    #      Correlation (covariance) between components of random vector
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            #
            # X is a sample matrix:
            # * I-th row corresponds to I-th observation
            # * J-th column corresponds to J-th variable
            #
            x = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(x)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(x)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(x)

            #
            # Three dependency measures are calculated:
            # * covariation
            # * Pearson correlation
            # * Spearman rank correlation
            #
            # Result is stored into C, with C[i,j] equal to correlation
            # (covariance) between I-th and J-th variables of X.
            #
            c = xalglib.covm(x)
            _TestResult = _TestResult and doc_print_test(c, [[1.80,0.60,-1.40],[0.60,0.70,-0.80],[-1.40,-0.80,14.70]], "real_matrix", 0.01)
            c = xalglib.pearsoncorrm(x)
            _TestResult = _TestResult and doc_print_test(c, [[1.000,0.535,-0.272],[0.535,1.000,-0.249],[-0.272,-0.249,1.000]], "real_matrix", 0.01)
            c = xalglib.spearmancorrm(x)
            _TestResult = _TestResult and doc_print_test(c, [[1.000,0.556,-0.306],[0.556,1.000,-0.750],[-0.306,-0.750,1.000]], "real_matrix", 0.01)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_d_cm                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_d_cm2
    #      Correlation (covariance) between two random vectors
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            #
            # X and Y are sample matrices:
            # * I-th row corresponds to I-th observation
            # * J-th column corresponds to J-th variable
            #
            x = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(x)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(x)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(x)
            y = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==3:
                spoil_mat_by_nan(y)
            if _spoil_scenario==4:
                spoil_mat_by_posinf(y)
            if _spoil_scenario==5:
                spoil_mat_by_neginf(y)

            #
            # Three dependency measures are calculated:
            # * covariation
            # * Pearson correlation
            # * Spearman rank correlation
            #
            # Result is stored into C, with C[i,j] equal to correlation
            # (covariance) between I-th variable of X and J-th variable of Y.
            #
            c = xalglib.covm2(x, y)
            _TestResult = _TestResult and doc_print_test(c, [[4.100,-3.250],[2.450,-1.500],[13.450,-5.750]], "real_matrix", 0.01)
            c = xalglib.pearsoncorrm2(x, y)
            _TestResult = _TestResult and doc_print_test(c, [[0.519,-0.699],[0.497,-0.518],[0.596,-0.433]], "real_matrix", 0.01)
            c = xalglib.spearmancorrm2(x, y)
            _TestResult = _TestResult and doc_print_test(c, [[0.541,-0.649],[0.216,-0.433],[0.433,-0.135]], "real_matrix", 0.01)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_d_cm2                   FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_t_base
    #      Tests ability to detect errors in inputs
    #
    _TestResult = True
    for _spoil_scenario in range(-1,34):
        try:

            #
            # first, we test short form of functions
            #
            x1 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x1)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x1)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x1)
            mean, variance, skewness, kurtosis = xalglib.samplemoments(x1)
            x2 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==3:
                spoil_vec_by_nan(x2)
            if _spoil_scenario==4:
                spoil_vec_by_posinf(x2)
            if _spoil_scenario==5:
                spoil_vec_by_neginf(x2)
            adev = xalglib.sampleadev(x2)
            x3 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==6:
                spoil_vec_by_nan(x3)
            if _spoil_scenario==7:
                spoil_vec_by_posinf(x3)
            if _spoil_scenario==8:
                spoil_vec_by_neginf(x3)
            v = xalglib.samplemedian(x3)
            x4 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==9:
                spoil_vec_by_nan(x4)
            if _spoil_scenario==10:
                spoil_vec_by_posinf(x4)
            if _spoil_scenario==11:
                spoil_vec_by_neginf(x4)
            p = 0.5
            if _spoil_scenario==12:
                p = float("nan")
            if _spoil_scenario==13:
                p = float("+inf")
            if _spoil_scenario==14:
                p = float("-inf")
            v = xalglib.samplepercentile(x4, p)

            #
            # and then we test full form
            #
            x5 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==15:
                spoil_vec_by_nan(x5)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(x5)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(x5)
            if _spoil_scenario==18:
                spoil_vec_by_deleting_element(x5)
            mean, variance, skewness, kurtosis = xalglib.samplemoments(x5, 10)
            x6 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==19:
                spoil_vec_by_nan(x6)
            if _spoil_scenario==20:
                spoil_vec_by_posinf(x6)
            if _spoil_scenario==21:
                spoil_vec_by_neginf(x6)
            if _spoil_scenario==22:
                spoil_vec_by_deleting_element(x6)
            adev = xalglib.sampleadev(x6, 10)
            x7 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==23:
                spoil_vec_by_nan(x7)
            if _spoil_scenario==24:
                spoil_vec_by_posinf(x7)
            if _spoil_scenario==25:
                spoil_vec_by_neginf(x7)
            if _spoil_scenario==26:
                spoil_vec_by_deleting_element(x7)
            v = xalglib.samplemedian(x7, 10)
            x8 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==27:
                spoil_vec_by_nan(x8)
            if _spoil_scenario==28:
                spoil_vec_by_posinf(x8)
            if _spoil_scenario==29:
                spoil_vec_by_neginf(x8)
            if _spoil_scenario==30:
                spoil_vec_by_deleting_element(x8)
            p = 0.5
            if _spoil_scenario==31:
                p = float("nan")
            if _spoil_scenario==32:
                p = float("+inf")
            if _spoil_scenario==33:
                p = float("-inf")
            v = xalglib.samplepercentile(x8, 10, p)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_t_base                  FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST basestat_t_covcorr
    #      Tests ability to detect errors in inputs
    #
    _TestResult = True
    for _spoil_scenario in range(-1,126):
        try:

            #
            # 2-sample short-form cov/corr are tested
            #
            x1 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x1)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x1)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x1)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x1)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x1)
            y1 = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y1)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y1)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y1)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y1)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y1)
            v = xalglib.cov2(x1, y1)
            x2 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==10:
                spoil_vec_by_nan(x2)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(x2)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(x2)
            if _spoil_scenario==13:
                spoil_vec_by_adding_element(x2)
            if _spoil_scenario==14:
                spoil_vec_by_deleting_element(x2)
            y2 = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==15:
                spoil_vec_by_nan(y2)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(y2)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(y2)
            if _spoil_scenario==18:
                spoil_vec_by_adding_element(y2)
            if _spoil_scenario==19:
                spoil_vec_by_deleting_element(y2)
            v = xalglib.pearsoncorr2(x2, y2)
            x3 = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==20:
                spoil_vec_by_nan(x3)
            if _spoil_scenario==21:
                spoil_vec_by_posinf(x3)
            if _spoil_scenario==22:
                spoil_vec_by_neginf(x3)
            if _spoil_scenario==23:
                spoil_vec_by_adding_element(x3)
            if _spoil_scenario==24:
                spoil_vec_by_deleting_element(x3)
            y3 = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==25:
                spoil_vec_by_nan(y3)
            if _spoil_scenario==26:
                spoil_vec_by_posinf(y3)
            if _spoil_scenario==27:
                spoil_vec_by_neginf(y3)
            if _spoil_scenario==28:
                spoil_vec_by_adding_element(y3)
            if _spoil_scenario==29:
                spoil_vec_by_deleting_element(y3)
            v = xalglib.spearmancorr2(x3, y3)

            #
            # 2-sample full-form cov/corr are tested
            #
            x1a = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==30:
                spoil_vec_by_nan(x1a)
            if _spoil_scenario==31:
                spoil_vec_by_posinf(x1a)
            if _spoil_scenario==32:
                spoil_vec_by_neginf(x1a)
            if _spoil_scenario==33:
                spoil_vec_by_deleting_element(x1a)
            y1a = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==34:
                spoil_vec_by_nan(y1a)
            if _spoil_scenario==35:
                spoil_vec_by_posinf(y1a)
            if _spoil_scenario==36:
                spoil_vec_by_neginf(y1a)
            if _spoil_scenario==37:
                spoil_vec_by_deleting_element(y1a)
            v = xalglib.cov2(x1a, y1a, 10)
            x2a = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==38:
                spoil_vec_by_nan(x2a)
            if _spoil_scenario==39:
                spoil_vec_by_posinf(x2a)
            if _spoil_scenario==40:
                spoil_vec_by_neginf(x2a)
            if _spoil_scenario==41:
                spoil_vec_by_deleting_element(x2a)
            y2a = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==42:
                spoil_vec_by_nan(y2a)
            if _spoil_scenario==43:
                spoil_vec_by_posinf(y2a)
            if _spoil_scenario==44:
                spoil_vec_by_neginf(y2a)
            if _spoil_scenario==45:
                spoil_vec_by_deleting_element(y2a)
            v = xalglib.pearsoncorr2(x2a, y2a, 10)
            x3a = [0,1,4,9,16,25,36,49,64,81]
            if _spoil_scenario==46:
                spoil_vec_by_nan(x3a)
            if _spoil_scenario==47:
                spoil_vec_by_posinf(x3a)
            if _spoil_scenario==48:
                spoil_vec_by_neginf(x3a)
            if _spoil_scenario==49:
                spoil_vec_by_deleting_element(x3a)
            y3a = [0,1,2,3,4,5,6,7,8,9]
            if _spoil_scenario==50:
                spoil_vec_by_nan(y3a)
            if _spoil_scenario==51:
                spoil_vec_by_posinf(y3a)
            if _spoil_scenario==52:
                spoil_vec_by_neginf(y3a)
            if _spoil_scenario==53:
                spoil_vec_by_deleting_element(y3a)
            v = xalglib.spearmancorr2(x3a, y3a, 10)

            #
            # vector short-form cov/corr are tested.
            #
            x4 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==54:
                spoil_mat_by_nan(x4)
            if _spoil_scenario==55:
                spoil_mat_by_posinf(x4)
            if _spoil_scenario==56:
                spoil_mat_by_neginf(x4)
            c = xalglib.covm(x4)
            x5 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==57:
                spoil_mat_by_nan(x5)
            if _spoil_scenario==58:
                spoil_mat_by_posinf(x5)
            if _spoil_scenario==59:
                spoil_mat_by_neginf(x5)
            c = xalglib.pearsoncorrm(x5)
            x6 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==60:
                spoil_mat_by_nan(x6)
            if _spoil_scenario==61:
                spoil_mat_by_posinf(x6)
            if _spoil_scenario==62:
                spoil_mat_by_neginf(x6)
            c = xalglib.spearmancorrm(x6)

            #
            # vector full-form cov/corr are tested.
            #
            x7 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==63:
                spoil_mat_by_nan(x7)
            if _spoil_scenario==64:
                spoil_mat_by_posinf(x7)
            if _spoil_scenario==65:
                spoil_mat_by_neginf(x7)
            if _spoil_scenario==66:
                spoil_mat_by_deleting_row(x7)
            if _spoil_scenario==67:
                spoil_mat_by_deleting_col(x7)
            c = xalglib.covm(x7, 5, 3)
            x8 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==68:
                spoil_mat_by_nan(x8)
            if _spoil_scenario==69:
                spoil_mat_by_posinf(x8)
            if _spoil_scenario==70:
                spoil_mat_by_neginf(x8)
            if _spoil_scenario==71:
                spoil_mat_by_deleting_row(x8)
            if _spoil_scenario==72:
                spoil_mat_by_deleting_col(x8)
            c = xalglib.pearsoncorrm(x8, 5, 3)
            x9 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==73:
                spoil_mat_by_nan(x9)
            if _spoil_scenario==74:
                spoil_mat_by_posinf(x9)
            if _spoil_scenario==75:
                spoil_mat_by_neginf(x9)
            if _spoil_scenario==76:
                spoil_mat_by_deleting_row(x9)
            if _spoil_scenario==77:
                spoil_mat_by_deleting_col(x9)
            c = xalglib.spearmancorrm(x9, 5, 3)

            #
            # cross-vector short-form cov/corr are tested.
            #
            x10 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==78:
                spoil_mat_by_nan(x10)
            if _spoil_scenario==79:
                spoil_mat_by_posinf(x10)
            if _spoil_scenario==80:
                spoil_mat_by_neginf(x10)
            y10 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==81:
                spoil_mat_by_nan(y10)
            if _spoil_scenario==82:
                spoil_mat_by_posinf(y10)
            if _spoil_scenario==83:
                spoil_mat_by_neginf(y10)
            c = xalglib.covm2(x10, y10)
            x11 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==84:
                spoil_mat_by_nan(x11)
            if _spoil_scenario==85:
                spoil_mat_by_posinf(x11)
            if _spoil_scenario==86:
                spoil_mat_by_neginf(x11)
            y11 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==87:
                spoil_mat_by_nan(y11)
            if _spoil_scenario==88:
                spoil_mat_by_posinf(y11)
            if _spoil_scenario==89:
                spoil_mat_by_neginf(y11)
            c = xalglib.pearsoncorrm2(x11, y11)
            x12 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==90:
                spoil_mat_by_nan(x12)
            if _spoil_scenario==91:
                spoil_mat_by_posinf(x12)
            if _spoil_scenario==92:
                spoil_mat_by_neginf(x12)
            y12 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==93:
                spoil_mat_by_nan(y12)
            if _spoil_scenario==94:
                spoil_mat_by_posinf(y12)
            if _spoil_scenario==95:
                spoil_mat_by_neginf(y12)
            c = xalglib.spearmancorrm2(x12, y12)

            #
            # cross-vector full-form cov/corr are tested.
            #
            x13 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==96:
                spoil_mat_by_nan(x13)
            if _spoil_scenario==97:
                spoil_mat_by_posinf(x13)
            if _spoil_scenario==98:
                spoil_mat_by_neginf(x13)
            if _spoil_scenario==99:
                spoil_mat_by_deleting_row(x13)
            if _spoil_scenario==100:
                spoil_mat_by_deleting_col(x13)
            y13 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==101:
                spoil_mat_by_nan(y13)
            if _spoil_scenario==102:
                spoil_mat_by_posinf(y13)
            if _spoil_scenario==103:
                spoil_mat_by_neginf(y13)
            if _spoil_scenario==104:
                spoil_mat_by_deleting_row(y13)
            if _spoil_scenario==105:
                spoil_mat_by_deleting_col(y13)
            c = xalglib.covm2(x13, y13, 5, 3, 2)
            x14 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==106:
                spoil_mat_by_nan(x14)
            if _spoil_scenario==107:
                spoil_mat_by_posinf(x14)
            if _spoil_scenario==108:
                spoil_mat_by_neginf(x14)
            if _spoil_scenario==109:
                spoil_mat_by_deleting_row(x14)
            if _spoil_scenario==110:
                spoil_mat_by_deleting_col(x14)
            y14 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==111:
                spoil_mat_by_nan(y14)
            if _spoil_scenario==112:
                spoil_mat_by_posinf(y14)
            if _spoil_scenario==113:
                spoil_mat_by_neginf(y14)
            if _spoil_scenario==114:
                spoil_mat_by_deleting_row(y14)
            if _spoil_scenario==115:
                spoil_mat_by_deleting_col(y14)
            c = xalglib.pearsoncorrm2(x14, y14, 5, 3, 2)
            x15 = [[1,0,1],[1,1,0],[-1,1,0],[-2,-1,1],[-1,0,9]]
            if _spoil_scenario==116:
                spoil_mat_by_nan(x15)
            if _spoil_scenario==117:
                spoil_mat_by_posinf(x15)
            if _spoil_scenario==118:
                spoil_mat_by_neginf(x15)
            if _spoil_scenario==119:
                spoil_mat_by_deleting_row(x15)
            if _spoil_scenario==120:
                spoil_mat_by_deleting_col(x15)
            y15 = [[2,3],[2,1],[-1,6],[-9,9],[7,1]]
            if _spoil_scenario==121:
                spoil_mat_by_nan(y15)
            if _spoil_scenario==122:
                spoil_mat_by_posinf(y15)
            if _spoil_scenario==123:
                spoil_mat_by_neginf(y15)
            if _spoil_scenario==124:
                spoil_mat_by_deleting_row(y15)
            if _spoil_scenario==125:
                spoil_mat_by_deleting_col(y15)
            c = xalglib.spearmancorrm2(x15, y15, 5, 3, 2)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("basestat_t_covcorr               FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_d_r1
    #      Real matrix inverse
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            a = [[1,-1],[1,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(a)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(a)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(a)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(a)
            a, info, rep = xalglib.rmatrixinverse(a)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(a, [[0.5,0.5],[-0.5,0.5]], "real_matrix", 0.00005)
            _TestResult = _TestResult and doc_print_test(rep.r1, 0.5, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(rep.rinf, 0.5, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matinv_d_r1                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_d_c1
    #      Complex matrix inverse
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            a = [[1j,-1],[1j,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(a)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(a)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(a)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(a)
            a, info, rep = xalglib.cmatrixinverse(a)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(a, [[-0.5j,-0.5j],[-0.5,0.5]], "complex_matrix", 0.00005)
            _TestResult = _TestResult and doc_print_test(rep.r1, 0.5, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(rep.rinf, 0.5, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matinv_d_c1                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_d_spd1
    #      SPD matrix inverse
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            a = [[2,1],[1,2]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(a)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(a)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(a)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(a)
            a, info, rep = xalglib.spdmatrixinverse(a)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(a, [[0.666666,-0.333333],[-0.333333,0.666666]], "real_matrix", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matinv_d_spd1                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_d_hpd1
    #      HPD matrix inverse
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            a = [[2,1],[1,2]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(a)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(a)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(a)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(a)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(a)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(a)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(a)
            a, info, rep = xalglib.hpdmatrixinverse(a)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(a, [[0.666666,-0.333333],[-0.333333,0.666666]], "complex_matrix", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matinv_d_hpd1                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_t_r1
    #      Real matrix inverse: singular matrix
    #
    _TestResult = True
    try:
        a = [[1,-1],[-2,2]]
        a, info, rep = xalglib.rmatrixinverse(a)
        _TestResult = _TestResult and doc_print_test(info, -3, "int")
        _TestResult = _TestResult and doc_print_test(rep.r1, 0.0, "real", 0.00005)
        _TestResult = _TestResult and doc_print_test(rep.rinf, 0.0, "real", 0.00005)
    except (RuntimeError, ValueError):
        _TestResult = False
    except:
        raise
    if not _TestResult:
        sys.stdout.write("matinv_t_r1                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_t_c1
    #      Complex matrix inverse: singular matrix
    #
    _TestResult = True
    try:
        a = [[1j,-1j],[-2,2]]
        a, info, rep = xalglib.cmatrixinverse(a)
        _TestResult = _TestResult and doc_print_test(info, -3, "int")
        _TestResult = _TestResult and doc_print_test(rep.r1, 0.0, "real", 0.00005)
        _TestResult = _TestResult and doc_print_test(rep.rinf, 0.0, "real", 0.00005)
    except (RuntimeError, ValueError):
        _TestResult = False
    except:
        raise
    if not _TestResult:
        sys.stdout.write("matinv_t_c1                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_e_spd1
    #      Attempt to use SPD function on nonsymmetrix matrix
    #
    _TestResult = True
    try:
        a = [[1,0],[1,1]]
        a, info, rep = xalglib.spdmatrixinverse(a)
        _TestResult = False
    except (RuntimeError, ValueError):
        pass
    except:
        raise
    if not _TestResult:
        sys.stdout.write("matinv_e_spd1                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matinv_e_hpd1
    #      Attempt to use SPD function on nonsymmetrix matrix
    #
    _TestResult = True
    try:
        a = [[1,0],[1,1]]
        a, info, rep = xalglib.hpdmatrixinverse(a)
        _TestResult = False
    except (RuntimeError, ValueError):
        pass
    except:
        raise
    if not _TestResult:
        sys.stdout.write("matinv_e_hpd1                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlbfgs_d_1
    #      Nonlinear optimization by L-BFGS
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # using LBFGS method.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.minlbfgscreate(1, x)
            xalglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlbfgsoptimize_g(state, function1_grad)
            x, rep = xalglib.minlbfgsresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlbfgs_d_1                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlbfgs_d_2
    #      Nonlinear optimization with additional settings and restarts
    #
    _TestResult = True
    for _spoil_scenario in range(-1,18):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # using LBFGS method.
            #
            # Several advanced techniques are demonstrated:
            # * upper limit on step size
            # * restart from new point
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            stpmax = 0.1
            if _spoil_scenario==12:
                stpmax = float("nan")
            if _spoil_scenario==13:
                stpmax = float("+inf")
            if _spoil_scenario==14:
                stpmax = float("-inf")
            maxits = 0

            # first run
            state = xalglib.minlbfgscreate(1, x)
            xalglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlbfgssetstpmax(state, stpmax)
            xalglib.minlbfgsoptimize_g(state, function1_grad)
            x, rep = xalglib.minlbfgsresults(state)

            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)

            # second run - algorithm is restarted
            x = [10,10]
            if _spoil_scenario==15:
                spoil_vec_by_nan(x)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(x)
            xalglib.minlbfgsrestartfrom(state, x)
            xalglib.minlbfgsoptimize_g(state, function1_grad)
            x, rep = xalglib.minlbfgsresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlbfgs_d_2                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST odesolver_d1
    #      Solving y'=-y with ODE solver
    #
    _TestResult = True
    for _spoil_scenario in range(-1,13):
        try:
            y = [1]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            x = [0, 1, 2, 3]
            if _spoil_scenario==4:
                spoil_vec_by_nan(x)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(x)
            eps = 0.00001
            if _spoil_scenario==7:
                eps = float("nan")
            if _spoil_scenario==8:
                eps = float("+inf")
            if _spoil_scenario==9:
                eps = float("-inf")
            h = 0
            if _spoil_scenario==10:
                h = float("nan")
            if _spoil_scenario==11:
                h = float("+inf")
            if _spoil_scenario==12:
                h = float("-inf")
            s = xalglib.odesolverrkck(y, x, eps, h)
            xalglib.odesolversolve(s, ode_function_1_diff);
            m, xtbl, ytbl, rep = xalglib.odesolverresults(s)
            _TestResult = _TestResult and doc_print_test(m, 4, "int")
            _TestResult = _TestResult and doc_print_test(xtbl, [0, 1, 2, 3], "real_vector", 0.005)
            _TestResult = _TestResult and doc_print_test(ytbl, [[1], [0.367], [0.135], [0.050]], "real_matrix", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("odesolver_d1                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST fft_complex_d1
    #      Complex FFT: simple example
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            #
            # first we demonstrate forward FFT:
            # [1i,1i,1i,1i] is converted to [4i, 0, 0, 0]
            #
            z = [1j,1j,1j,1j]
            if _spoil_scenario==0:
                spoil_vec_by_nan(z)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(z)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(z)
            z = xalglib.fftc1d(z)
            _TestResult = _TestResult and doc_print_test(z, [4j,0,0,0], "complex_vector", 0.0001)

            #
            # now we convert [4i, 0, 0, 0] back to [1i,1i,1i,1i]
            # with backward FFT
            #
            z = xalglib.fftc1dinv(z)
            _TestResult = _TestResult and doc_print_test(z, [1j,1j,1j,1j], "complex_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("fft_complex_d1                   FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST fft_complex_d2
    #      Complex FFT: advanced example
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            #
            # first we demonstrate forward FFT:
            # [0,1,0,1i] is converted to [1+1i, -1-1i, -1-1i, 1+1i]
            #
            z = [0,1,0,1j]
            if _spoil_scenario==0:
                spoil_vec_by_nan(z)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(z)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(z)
            z = xalglib.fftc1d(z)
            _TestResult = _TestResult and doc_print_test(z, [1+1j, -1-1j, -1-1j, 1+1j], "complex_vector", 0.0001)

            #
            # now we convert result back with backward FFT
            #
            z = xalglib.fftc1dinv(z)
            _TestResult = _TestResult and doc_print_test(z, [0,1,0,1j], "complex_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("fft_complex_d2                   FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST fft_real_d1
    #      Real FFT: simple example
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            #
            # first we demonstrate forward FFT:
            # [1,1,1,1] is converted to [4, 0, 0, 0]
            #
            x = [1,1,1,1]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            f = xalglib.fftr1d(x)
            _TestResult = _TestResult and doc_print_test(f, [4,0,0,0], "complex_vector", 0.0001)

            #
            # now we convert [4, 0, 0, 0] back to [1,1,1,1]
            # with backward FFT
            #
            x2 = xalglib.fftr1dinv(f)
            _TestResult = _TestResult and doc_print_test(x2, [1,1,1,1], "real_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("fft_real_d1                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST fft_real_d2
    #      Real FFT: advanced example
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            #
            # first we demonstrate forward FFT:
            # [1,2,3,4] is converted to [10, -2+2i, -2, -2-2i]
            #
            # note that output array is self-adjoint:
            # * f[0] = conj(f[0])
            # * f[1] = conj(f[3])
            # * f[2] = conj(f[2])
            #
            x = [1,2,3,4]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            f = xalglib.fftr1d(x)
            _TestResult = _TestResult and doc_print_test(f, [10, -2+2j, -2, -2-2j], "complex_vector", 0.0001)

            #
            # now we convert [10, -2+2i, -2, -2-2i] back to [1,2,3,4]
            #
            x2 = xalglib.fftr1dinv(f)
            _TestResult = _TestResult and doc_print_test(x2, [1,2,3,4], "real_vector", 0.0001)

            #
            # remember that F is self-adjoint? It means that we can pass just half
            # (slightly larger than half) of F to inverse real FFT and still get our result.
            #
            # I.e. instead [10, -2+2i, -2, -2-2i] we pass just [10, -2+2i, -2] and everything works!
            #
            # NOTE: in this case we should explicitly pass array length (which is 4) to ALGLIB;
            # if not, it will automatically use array length to determine FFT size and
            # will erroneously make half-length FFT.
            #
            f = [10, -2+2j, -2]
            x2 = xalglib.fftr1dinv(f, 4)
            _TestResult = _TestResult and doc_print_test(x2, [1,2,3,4], "real_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("fft_real_d2                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST fft_complex_e1
    #      error detection in backward FFT
    #
    _TestResult = True
    for _spoil_scenario in range(-1,3):
        try:
            z = [0,2,0,-2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(z)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(z)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(z)
            z = xalglib.fftc1dinv(z)
            _TestResult = _TestResult and doc_print_test(z, [0,1j,0,-1j], "complex_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("fft_complex_e1                   FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST autogk_d1
    #      Integrating f=exp(x) by adaptive integrator
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            #
            # This example demonstrates integration of f=exp(x) on [0,1]:
            # * first, autogkstate is initialized
            # * then we call integration function
            # * and finally we obtain results with autogkresults() call
            #
            a = 0
            if _spoil_scenario==0:
                a = float("nan")
            if _spoil_scenario==1:
                a = float("+inf")
            if _spoil_scenario==2:
                a = float("-inf")
            b = 1
            if _spoil_scenario==3:
                b = float("nan")
            if _spoil_scenario==4:
                b = float("+inf")
            if _spoil_scenario==5:
                b = float("-inf")

            s = xalglib.autogksmooth(a, b)
            xalglib.autogkintegrate(s, int_function_1_func);
            v, rep = xalglib.autogkresults(s)

            _TestResult = _TestResult and doc_print_test(v, 1.7182, "real", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("autogk_d1                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_d_calcdiff
    #      Interpolation and differentiation using barycentric representation
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # Here we demonstrate polynomial interpolation and differentiation
            # of y=x^2-x sampled at [0,1,2]. Barycentric representation of polynomial is used.
            #
            x = [0,1,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [0,0,2]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==10:
                t = float("+inf")
            if _spoil_scenario==11:
                t = float("-inf")

            # barycentric model is created
            p = xalglib.polynomialbuild(x, y)

            # barycentric interpolation is demonstrated
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)

            # barycentric differentation is demonstrated
            v, dv = xalglib.barycentricdiff1(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(dv, -3.0, "real", 0.00005)

            # second derivatives with barycentric representation
            v, dv = xalglib.barycentricdiff1(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(dv, -3.0, "real", 0.00005)
            v, dv, d2v = xalglib.barycentricdiff2(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(dv, -3.0, "real", 0.00005)
            _TestResult = _TestResult and doc_print_test(d2v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_d_calcdiff                FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_d_conv
    #      Conversion between power basis and barycentric representation
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            #
            # Here we demonstrate conversion of y=x^2-x
            # between power basis and barycentric representation.
            #
            a = [0,-1,+1]
            if _spoil_scenario==0:
                spoil_vec_by_nan(a)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(a)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(a)
            t = 2
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")

            #
            # a=[0,-1,+1] is decomposition of y=x^2-x in the power basis:
            #
            #     y = 0 - 1*x + 1*x^2
            #
            # We convert it to the barycentric form.
            #
            p = xalglib.polynomialpow2bar(a)

            # now we have barycentric interpolation; we can use it for interpolation
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.005)

            # we can also convert back from barycentric representation to power basis
            a2 = xalglib.polynomialbar2pow(p)
            _TestResult = _TestResult and doc_print_test(a2, [0,-1,+1], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_d_conv                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_d_spec
    #      Polynomial interpolation on special grids (equidistant, Chebyshev I/II)
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            #
            # Temporaries:
            # * values of y=x^2-x sampled at three special grids:
            #   * equdistant grid spanning [0,2],     x[i] = 2*i/(N-1), i=0..N-1
            #   * Chebyshev-I grid spanning [-1,+1],  x[i] = 1 + Cos(PI*(2*i+1)/(2*n)), i=0..N-1
            #   * Chebyshev-II grid spanning [-1,+1], x[i] = 1 + Cos(PI*i/(n-1)), i=0..N-1
            # * barycentric interpolants for these three grids
            # * vectors to store coefficients of quadratic representation
            #
            y_eqdist = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y_eqdist)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y_eqdist)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y_eqdist)
            y_cheb1 = [-0.116025,0.000000,1.616025]
            if _spoil_scenario==3:
                spoil_vec_by_nan(y_cheb1)
            if _spoil_scenario==4:
                spoil_vec_by_posinf(y_cheb1)
            if _spoil_scenario==5:
                spoil_vec_by_neginf(y_cheb1)
            y_cheb2 = [0,0,2]
            if _spoil_scenario==6:
                spoil_vec_by_nan(y_cheb2)
            if _spoil_scenario==7:
                spoil_vec_by_posinf(y_cheb2)
            if _spoil_scenario==8:
                spoil_vec_by_neginf(y_cheb2)

            #
            # First, we demonstrate construction of barycentric interpolants on
            # special grids. We unpack power representation to ensure that
            # interpolant was built correctly.
            #
            # In all three cases we should get same quadratic function.
            #
            p_eqdist = xalglib.polynomialbuildeqdist(0.0, 2.0, y_eqdist)
            a_eqdist = xalglib.polynomialbar2pow(p_eqdist)
            _TestResult = _TestResult and doc_print_test(a_eqdist, [0,-1,+1], "real_vector", 0.00005)

            p_cheb1 = xalglib.polynomialbuildcheb1(-1, +1, y_cheb1)
            a_cheb1 = xalglib.polynomialbar2pow(p_cheb1)
            _TestResult = _TestResult and doc_print_test(a_cheb1, [0,-1,+1], "real_vector", 0.00005)

            p_cheb2 = xalglib.polynomialbuildcheb2(-1, +1, y_cheb2)
            a_cheb2 = xalglib.polynomialbar2pow(p_cheb2)
            _TestResult = _TestResult and doc_print_test(a_cheb2, [0,-1,+1], "real_vector", 0.00005)

            #
            # Now we demonstrate polynomial interpolation without construction 
            # of the barycentricinterpolant structure.
            #
            # We calculate interpolant value at x=-2.
            # In all three cases we should get same f=6
            #
            t = -2
            if _spoil_scenario==9:
                t = float("+inf")
            if _spoil_scenario==10:
                t = float("-inf")
            v = xalglib.polynomialcalceqdist(0.0, 2.0, y_eqdist, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)

            v = xalglib.polynomialcalccheb1(-1, +1, y_cheb1, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)

            v = xalglib.polynomialcalccheb2(-1, +1, y_cheb2, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_d_spec                    FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_1
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,10):
        try:
            x = [0,1,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x)
            y = [0,0,2]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==8:
                t = float("+inf")
            if _spoil_scenario==9:
                t = float("-inf")
            p = xalglib.polynomialbuild(x, y, 3)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_1                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_2
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            p = xalglib.polynomialbuildeqdist(0.0, 2.0, y, 3)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_2                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_3
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            y = [-0.116025,0.000000,1.616025]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            p = xalglib.polynomialbuildcheb1(-1.0, +1.0, y, 3)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_3                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_4
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -2
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            a = -1
            if _spoil_scenario==6:
                a = float("nan")
            if _spoil_scenario==7:
                a = float("+inf")
            if _spoil_scenario==8:
                a = float("-inf")
            b = +1
            if _spoil_scenario==9:
                b = float("nan")
            if _spoil_scenario==10:
                b = float("+inf")
            if _spoil_scenario==11:
                b = float("-inf")
            p = xalglib.polynomialbuildcheb2(a, b, y, 3)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_4                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_5
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            v = xalglib.polynomialcalceqdist(0.0, 2.0, y, 3, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_5                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_6
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            y = [-0.116025,0.000000,1.616025]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -1
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            a = -1
            if _spoil_scenario==6:
                a = float("nan")
            if _spoil_scenario==7:
                a = float("+inf")
            if _spoil_scenario==8:
                a = float("-inf")
            b = +1
            if _spoil_scenario==9:
                b = float("nan")
            if _spoil_scenario==10:
                b = float("+inf")
            if _spoil_scenario==11:
                b = float("-inf")
            v = xalglib.polynomialcalccheb1(a, b, y, 3, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_6                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_7
    #      Polynomial interpolation, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(y)
            t = -2
            if _spoil_scenario==4:
                t = float("+inf")
            if _spoil_scenario==5:
                t = float("-inf")
            a = -1
            if _spoil_scenario==6:
                a = float("nan")
            if _spoil_scenario==7:
                a = float("+inf")
            if _spoil_scenario==8:
                a = float("-inf")
            b = +1
            if _spoil_scenario==9:
                b = float("nan")
            if _spoil_scenario==10:
                b = float("+inf")
            if _spoil_scenario==11:
                b = float("-inf")
            v = xalglib.polynomialcalccheb2(a, b, y, 3, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_7                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_8
    #      Polynomial interpolation: y=x^2-x, equidistant grid, barycentric form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -1
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            p = xalglib.polynomialbuildeqdist(0.0, 2.0, y)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_8                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_9
    #      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind), barycentric form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            y = [-0.116025,0.000000,1.616025]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -1
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            a = -1
            if _spoil_scenario==5:
                a = float("nan")
            if _spoil_scenario==6:
                a = float("+inf")
            if _spoil_scenario==7:
                a = float("-inf")
            b = +1
            if _spoil_scenario==8:
                b = float("nan")
            if _spoil_scenario==9:
                b = float("+inf")
            if _spoil_scenario==10:
                b = float("-inf")
            p = xalglib.polynomialbuildcheb1(a, b, y)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_9                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_10
    #      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind), barycentric form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -2
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            a = -1
            if _spoil_scenario==5:
                a = float("nan")
            if _spoil_scenario==6:
                a = float("+inf")
            if _spoil_scenario==7:
                a = float("-inf")
            b = +1
            if _spoil_scenario==8:
                b = float("nan")
            if _spoil_scenario==9:
                b = float("+inf")
            if _spoil_scenario==10:
                b = float("-inf")
            p = xalglib.polynomialbuildcheb2(a, b, y)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_10                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_11
    #      Polynomial interpolation: y=x^2-x, equidistant grid
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -1
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            v = xalglib.polynomialcalceqdist(0.0, 2.0, y, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_11                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_12
    #      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind)
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            y = [-0.116025,0.000000,1.616025]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -1
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            a = -1
            if _spoil_scenario==5:
                a = float("nan")
            if _spoil_scenario==6:
                a = float("+inf")
            if _spoil_scenario==7:
                a = float("-inf")
            b = +1
            if _spoil_scenario==8:
                b = float("nan")
            if _spoil_scenario==9:
                b = float("+inf")
            if _spoil_scenario==10:
                b = float("-inf")
            v = xalglib.polynomialcalccheb1(a, b, y, t)
            _TestResult = _TestResult and doc_print_test(v, 2.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_12                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST polint_t_13
    #      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind)
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            y = [0,0,2]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            t = -2
            if _spoil_scenario==3:
                t = float("+inf")
            if _spoil_scenario==4:
                t = float("-inf")
            a = -1
            if _spoil_scenario==5:
                a = float("nan")
            if _spoil_scenario==6:
                a = float("+inf")
            if _spoil_scenario==7:
                a = float("-inf")
            b = +1
            if _spoil_scenario==8:
                b = float("nan")
            if _spoil_scenario==9:
                b = float("+inf")
            if _spoil_scenario==10:
                b = float("-inf")
            v = xalglib.polynomialcalccheb2(a, b, y, t)
            _TestResult = _TestResult and doc_print_test(v, 6.0, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("polint_t_13                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST spline1d_d_linear
    #      Piecewise linear spline interpolation
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # We use piecewise linear spline to interpolate f(x)=x^2 sampled 
            # at 5 equidistant nodes on [-1,+1].
            #
            x = [-1.0,-0.5,0.0,+0.5,+1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [+1.0,0.25,0.0,0.25,+1.0]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            t = 0.25
            if _spoil_scenario==10:
                t = float("+inf")
            if _spoil_scenario==11:
                t = float("-inf")

            # build spline
            s = xalglib.spline1dbuildlinear(x, y)

            # calculate S(0.25) - it is quite different from 0.25^2=0.0625
            v = xalglib.spline1dcalc(s, t)
            _TestResult = _TestResult and doc_print_test(v, 0.125, "real", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("spline1d_d_linear                FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST spline1d_d_cubic
    #      Cubic spline interpolation
    #
    _TestResult = True
    for _spoil_scenario in range(-1,10):
        try:
            #
            # We use cubic spline to interpolate f(x)=x^2 sampled 
            # at 5 equidistant nodes on [-1,+1].
            #
            # First, we use default boundary conditions ("parabolically terminated
            # spline") because cubic spline built with such boundary conditions 
            # will exactly reproduce any quadratic f(x).
            #
            # Then we try to use natural boundary conditions
            #     d2S(-1)/dx^2 = 0.0
            #     d2S(+1)/dx^2 = 0.0
            # and see that such spline interpolated f(x) with small error.
            #
            x = [-1.0,-0.5,0.0,+0.5,+1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x)
            y = [+1.0,0.25,0.0,0.25,+1.0]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)
            t = 0.25
            if _spoil_scenario==8:
                t = float("+inf")
            if _spoil_scenario==9:
                t = float("-inf")
            natural_bound_type = 2
            #
            # Test exact boundary conditions: build S(x), calculare S(0.25)
            # (almost same as original function)
            #
            s = xalglib.spline1dbuildcubic(x, y)
            v = xalglib.spline1dcalc(s, t)
            _TestResult = _TestResult and doc_print_test(v, 0.0625, "real", 0.00001)

            #
            # Test natural boundary conditions: build S(x), calculare S(0.25)
            # (small interpolation error)
            #
            s = xalglib.spline1dbuildcubic(x, y, 5, natural_bound_type, 0.0, natural_bound_type, 0.0)
            v = xalglib.spline1dcalc(s, t)
            _TestResult = _TestResult and doc_print_test(v, 0.0580, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("spline1d_d_cubic                 FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST spline1d_d_griddiff
    #      Differentiation on the grid using cubic splines
    #
    _TestResult = True
    for _spoil_scenario in range(-1,10):
        try:
            #
            # We use cubic spline to do grid differentiation, i.e. having
            # values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
            # we calculate derivatives of cubic spline at nodes WITHOUT
            # CONSTRUCTION OF SPLINE OBJECT.
            #
            # There are efficient functions spline1dgriddiffcubic() and
            # spline1dgriddiff2cubic() for such calculations.
            #
            # We use default boundary conditions ("parabolically terminated
            # spline") because cubic spline built with such boundary conditions 
            # will exactly reproduce any quadratic f(x).
            #
            # Actually, we could use natural conditions, but we feel that 
            # spline which exactly reproduces f() will show us more 
            # understandable results.
            #
            x = [-1.0,-0.5,0.0,+0.5,+1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [+1.0,0.25,0.0,0.25,+1.0]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)

            #
            # We calculate first derivatives: they must be equal to 2*x
            #
            d1 = xalglib.spline1dgriddiffcubic(x, y)
            _TestResult = _TestResult and doc_print_test(d1, [-2.0, -1.0, 0.0, +1.0, +2.0], "real_vector", 0.0001)

            #
            # Now test griddiff2, which returns first AND second derivatives.
            # First derivative is 2*x, second is equal to 2.0
            #
            d1, d2 = xalglib.spline1dgriddiff2cubic(x, y)
            _TestResult = _TestResult and doc_print_test(d1, [-2.0, -1.0, 0.0, +1.0, +2.0], "real_vector", 0.0001)
            _TestResult = _TestResult and doc_print_test(d2, [ 2.0,  2.0, 2.0,  2.0,  2.0], "real_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("spline1d_d_griddiff              FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST spline1d_d_convdiff
    #      Resampling using cubic splines
    #
    _TestResult = True
    for _spoil_scenario in range(-1,11):
        try:
            #
            # We use cubic spline to do resampling, i.e. having
            # values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
            # we calculate values/derivatives of cubic spline on 
            # another grid (equidistant with 9 nodes on [-1,+1])
            # WITHOUT CONSTRUCTION OF SPLINE OBJECT.
            #
            # There are efficient functions spline1dconvcubic(),
            # spline1dconvdiffcubic() and spline1dconvdiff2cubic() 
            # for such calculations.
            #
            # We use default boundary conditions ("parabolically terminated
            # spline") because cubic spline built with such boundary conditions 
            # will exactly reproduce any quadratic f(x).
            #
            # Actually, we could use natural conditions, but we feel that 
            # spline which exactly reproduces f() will show us more 
            # understandable results.
            #
            x_old = [-1.0,-0.5,0.0,+0.5,+1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x_old)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x_old)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x_old)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x_old)
            y_old = [+1.0,0.25,0.0,0.25,+1.0]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y_old)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y_old)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y_old)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y_old)
            x_new = [-1.00,-0.75,-0.50,-0.25,0.00,+0.25,+0.50,+0.75,+1.00]
            if _spoil_scenario==8:
                spoil_vec_by_nan(x_new)
            if _spoil_scenario==9:
                spoil_vec_by_posinf(x_new)
            if _spoil_scenario==10:
                spoil_vec_by_neginf(x_new)

            #
            # First, conversion without differentiation.
            #
            #
            y_new = xalglib.spline1dconvcubic(x_old, y_old, x_new)
            _TestResult = _TestResult and doc_print_test(y_new, [1.0000, 0.5625, 0.2500, 0.0625, 0.0000, 0.0625, 0.2500, 0.5625, 1.0000], "real_vector", 0.0001)

            #
            # Then, conversion with differentiation (first derivatives only)
            #
            #
            y_new, d1_new = xalglib.spline1dconvdiffcubic(x_old, y_old, x_new)
            _TestResult = _TestResult and doc_print_test(y_new, [1.0000, 0.5625, 0.2500, 0.0625, 0.0000, 0.0625, 0.2500, 0.5625, 1.0000], "real_vector", 0.0001)
            _TestResult = _TestResult and doc_print_test(d1_new, [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0], "real_vector", 0.0001)

            #
            # Finally, conversion with first and second derivatives
            #
            #
            y_new, d1_new, d2_new = xalglib.spline1dconvdiff2cubic(x_old, y_old, x_new)
            _TestResult = _TestResult and doc_print_test(y_new, [1.0000, 0.5625, 0.2500, 0.0625, 0.0000, 0.0625, 0.2500, 0.5625, 1.0000], "real_vector", 0.0001)
            _TestResult = _TestResult and doc_print_test(d1_new, [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0], "real_vector", 0.0001)
            _TestResult = _TestResult and doc_print_test(d2_new, [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0], "real_vector", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("spline1d_d_convdiff              FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_d_v
    #      Nonlinear least squares optimization using function vector only
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
            #
            #     f0(x0,x1) = 10*(x0+3)^2
            #     f1(x0,x1) = (x1-3)^2
            #
            # using "V" mode of the Levenberg-Marquardt optimizer.
            #
            # Optimization algorithm uses:
            # * function vector f[] = {f1,f2}
            #
            # No other information (Jacobian, gradient, etc.) is needed.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.minlmcreatev(2, x, 0.0001)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_v(state, function1_fvec)
            x, rep = xalglib.minlmresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_d_v                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_d_vj
    #      Nonlinear least squares optimization using function vector and Jacobian
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
            #
            #     f0(x0,x1) = 10*(x0+3)^2
            #     f1(x0,x1) = (x1-3)^2
            #
            # using "VJ" mode of the Levenberg-Marquardt optimizer.
            #
            # Optimization algorithm uses:
            # * function vector f[] = {f1,f2}
            # * Jacobian matrix J = {dfi/dxj}.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.minlmcreatevj(2, x)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_vj(state, function1_fvec, function1_jac)
            x, rep = xalglib.minlmresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_d_vj                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_d_fgh
    #      Nonlinear Hessian-based optimization for general functions
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # This example demonstrates minimization of F(x0,x1) = 100*(x0+3)^4+(x1-3)^4
            # using "FGH" mode of the Levenberg-Marquardt optimizer.
            #
            # F is treated like a monolitic function without internal structure,
            # i.e. we do NOT represent it as a sum of squares.
            #
            # Optimization algorithm uses:
            # * function value F(x0,x1)
            # * gradient G={dF/dxi}
            # * Hessian H={d2F/(dxi*dxj)}
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.minlmcreatefgh(x)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_fgh(state, function1_func, function1_grad, function1_hess)
            x, rep = xalglib.minlmresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_d_fgh                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_d_restarts
    #      Efficient restarts of LM optimizer
    #
    _TestResult = True
    for _spoil_scenario in range(-1,15):
        try:
            #
            # This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
            #
            #     f0(x0,x1) = 10*(x0+3)^2
            #     f1(x0,x1) = (x1-3)^2
            #
            # using several starting points and efficient restarts.
            #
            epsg = 0.0000000001
            if _spoil_scenario==0:
                epsg = float("nan")
            if _spoil_scenario==1:
                epsg = float("+inf")
            if _spoil_scenario==2:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==3:
                epsf = float("nan")
            if _spoil_scenario==4:
                epsf = float("+inf")
            if _spoil_scenario==5:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==6:
                epsx = float("nan")
            if _spoil_scenario==7:
                epsx = float("+inf")
            if _spoil_scenario==8:
                epsx = float("-inf")
            maxits = 0

            #
            # create optimizer using minlmcreatev()
            #
            x = [10,10]
            if _spoil_scenario==9:
                spoil_vec_by_nan(x)
            if _spoil_scenario==10:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==11:
                spoil_vec_by_neginf(x)
            state = xalglib.minlmcreatev(2, x, 0.0001)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_v(state, function1_fvec)
            x, rep = xalglib.minlmresults(state)
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)

            #
            # restart optimizer using minlmrestartfrom()
            #
            # we can use different starting point, different function,
            # different stopping conditions, but problem size
            # must remain unchanged.
            #
            x = [4,4]
            if _spoil_scenario==12:
                spoil_vec_by_nan(x)
            if _spoil_scenario==13:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==14:
                spoil_vec_by_neginf(x)
            xalglib.minlmrestartfrom(state, x)
            xalglib.minlmoptimize_v(state, function2_fvec)
            x, rep = xalglib.minlmresults(state)
            _TestResult = _TestResult and doc_print_test(x, [0,1], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_d_restarts                 FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_t_1
    #      Nonlinear least squares optimization, FJ scheme (obsolete, but supported)
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0
            state = xalglib.minlmcreatefj(2, x)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_fj(state, function1_func, function1_jac)
            x, rep = xalglib.minlmresults(state)
            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_t_1                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minlm_t_2
    #      Nonlinear least squares optimization, FGJ scheme (obsolete, but supported)
    #
    sys.stdout.write("50/78\n")
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0
            state = xalglib.minlmcreatefgj(2, x)
            xalglib.minlmsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minlmoptimize_fgj(state, function1_func, function1_grad, function1_jac)
            x, rep = xalglib.minlmresults(state)
            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,+3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minlm_t_2                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_f
    #      Nonlinear fitting using function value only
    #
    _TestResult = True
    for _spoil_scenario in range(-1,27):
        try:
            #
            # In this example we demonstrate exponential fitting
            # by f(x) = exp(-c*x^2)
            # using function value only.
            #
            # Gradient is estimated using combination of numerical differences
            # and secant updates. diffstep variable stores differentiation step 
            # (we have to tell algorithm what step to use).
            #
            x = [[-1],[-0.8],[-0.6],[-0.4],[-0.2],[0],[0.2],[0.4],[0.6],[0.8],[1.0]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(x)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(x)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(x)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(x)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(x)
            y = [0.223130, 0.382893, 0.582748, 0.786628, 0.941765, 1.000000, 0.941765, 0.786628, 0.582748, 0.382893, 0.223130]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            c = [0.3]
            if _spoil_scenario==10:
                spoil_vec_by_nan(c)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(c)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(c)
            epsf = 0
            if _spoil_scenario==13:
                epsf = float("nan")
            if _spoil_scenario==14:
                epsf = float("+inf")
            if _spoil_scenario==15:
                epsf = float("-inf")
            epsx = 0.000001
            if _spoil_scenario==16:
                epsx = float("nan")
            if _spoil_scenario==17:
                epsx = float("+inf")
            if _spoil_scenario==18:
                epsx = float("-inf")
            maxits = 0
            diffstep = 0.0001
            if _spoil_scenario==19:
                diffstep = float("nan")
            if _spoil_scenario==20:
                diffstep = float("+inf")
            if _spoil_scenario==21:
                diffstep = float("-inf")

            #
            # Fitting without weights
            #
            state = xalglib.lsfitcreatef(x, y, c, diffstep)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_f(state, function_cx_1_func)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)

            #
            # Fitting with weights
            # (you can change weights and see how it changes result)
            #
            w = [1,1,1,1,1,1,1,1,1,1,1]
            if _spoil_scenario==22:
                spoil_vec_by_nan(w)
            if _spoil_scenario==23:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==24:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==25:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==26:
                spoil_vec_by_deleting_element(w)
            state = xalglib.lsfitcreatewf(x, y, w, c, diffstep)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_f(state, function_cx_1_func)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_f                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_fg
    #      Nonlinear fitting using gradient
    #
    _TestResult = True
    for _spoil_scenario in range(-1,24):
        try:
            #
            # In this example we demonstrate exponential fitting
            # by f(x) = exp(-c*x^2)
            # using function value and gradient (with respect to c).
            #
            x = [[-1],[-0.8],[-0.6],[-0.4],[-0.2],[0],[0.2],[0.4],[0.6],[0.8],[1.0]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(x)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(x)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(x)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(x)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(x)
            y = [0.223130, 0.382893, 0.582748, 0.786628, 0.941765, 1.000000, 0.941765, 0.786628, 0.582748, 0.382893, 0.223130]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            c = [0.3]
            if _spoil_scenario==10:
                spoil_vec_by_nan(c)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(c)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(c)
            epsf = 0
            if _spoil_scenario==13:
                epsf = float("nan")
            if _spoil_scenario==14:
                epsf = float("+inf")
            if _spoil_scenario==15:
                epsf = float("-inf")
            epsx = 0.000001
            if _spoil_scenario==16:
                epsx = float("nan")
            if _spoil_scenario==17:
                epsx = float("+inf")
            if _spoil_scenario==18:
                epsx = float("-inf")
            maxits = 0

            #
            # Fitting without weights
            #
            state = xalglib.lsfitcreatefg(x, y, c, True)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_fg(state, function_cx_1_func, function_cx_1_grad)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)

            #
            # Fitting with weights
            # (you can change weights and see how it changes result)
            #
            w = [1,1,1,1,1,1,1,1,1,1,1]
            if _spoil_scenario==19:
                spoil_vec_by_nan(w)
            if _spoil_scenario==20:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==21:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==22:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==23:
                spoil_vec_by_deleting_element(w)
            state = xalglib.lsfitcreatewfg(x, y, w, c, True)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_fg(state, function_cx_1_func, function_cx_1_grad)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_fg                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_fgh
    #      Nonlinear fitting using gradient and Hessian
    #
    _TestResult = True
    for _spoil_scenario in range(-1,24):
        try:
            #
            # In this example we demonstrate exponential fitting
            # by f(x) = exp(-c*x^2)
            # using function value, gradient and Hessian (with respect to c)
            #
            x = [[-1],[-0.8],[-0.6],[-0.4],[-0.2],[0],[0.2],[0.4],[0.6],[0.8],[1.0]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(x)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(x)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(x)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(x)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(x)
            y = [0.223130, 0.382893, 0.582748, 0.786628, 0.941765, 1.000000, 0.941765, 0.786628, 0.582748, 0.382893, 0.223130]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            c = [0.3]
            if _spoil_scenario==10:
                spoil_vec_by_nan(c)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(c)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(c)
            epsf = 0
            if _spoil_scenario==13:
                epsf = float("nan")
            if _spoil_scenario==14:
                epsf = float("+inf")
            if _spoil_scenario==15:
                epsf = float("-inf")
            epsx = 0.000001
            if _spoil_scenario==16:
                epsx = float("nan")
            if _spoil_scenario==17:
                epsx = float("+inf")
            if _spoil_scenario==18:
                epsx = float("-inf")
            maxits = 0

            #
            # Fitting without weights
            #
            state = xalglib.lsfitcreatefgh(x, y, c)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_fgh(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)

            #
            # Fitting with weights
            # (you can change weights and see how it changes result)
            #
            w = [1,1,1,1,1,1,1,1,1,1,1]
            if _spoil_scenario==19:
                spoil_vec_by_nan(w)
            if _spoil_scenario==20:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==21:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==22:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==23:
                spoil_vec_by_deleting_element(w)
            state = xalglib.lsfitcreatewfgh(x, y, w, c)
            xalglib.lsfitsetcond(state, epsf, epsx, maxits)
            xalglib.lsfitfit_fgh(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess)
            info, c, rep = xalglib.lsfitresults(state)
            _TestResult = _TestResult and doc_print_test(info, 2, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.5], "real_vector", 0.05)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_fgh                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_lin
    #      Unconstrained (general) linear least squares fitting with and without weights
    #
    _TestResult = True
    for _spoil_scenario in range(-1,13):
        try:
            #
            # In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
            #
            # We have:
            # * y - vector of experimental data
            # * fmatrix -  matrix of basis functions calculated at sample points
            #              Actually, we have only one basis function F0 = exp(0.5*x).
            #
            fmatrix = [[0.606531],[0.670320],[0.740818],[0.818731],[0.904837],[1.000000],[1.105171],[1.221403],[1.349859],[1.491825],[1.648721]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(fmatrix)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(fmatrix)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(fmatrix)
            y = [1.133719, 1.306522, 1.504604, 1.554663, 1.884638, 2.072436, 2.257285, 2.534068, 2.622017, 2.897713, 3.219371]
            if _spoil_scenario==3:
                spoil_vec_by_nan(y)
            if _spoil_scenario==4:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==5:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==6:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)

            #
            # Linear fitting without weights
            #
            info, c, rep = xalglib.lsfitlinear(y, fmatrix)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.98650], "real_vector", 0.00005)

            #
            # Linear fitting with individual weights.
            # Slightly different result is returned.
            #
            w = [1.414213, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            if _spoil_scenario==8:
                spoil_vec_by_nan(w)
            if _spoil_scenario==9:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==10:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==11:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==12:
                spoil_vec_by_deleting_element(w)
            info, c, rep = xalglib.lsfitlinearw(y, w, fmatrix)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(c, [1.983354], "real_vector", 0.00005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_lin                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_linc
    #      Constrained (general) linear least squares fitting with and without weights
    #
    _TestResult = True
    for _spoil_scenario in range(-1,20):
        try:
            #
            # In this example we demonstrate linear fitting by f(x|a,b) = a*x+b
            # with simple constraint f(0)=0.
            #
            # We have:
            # * y - vector of experimental data
            # * fmatrix -  matrix of basis functions sampled at [0,1] with step 0.2:
            #                  [ 1.0   0.0 ]
            #                  [ 1.0   0.2 ]
            #                  [ 1.0   0.4 ]
            #                  [ 1.0   0.6 ]
            #                  [ 1.0   0.8 ]
            #                  [ 1.0   1.0 ]
            #              first column contains value of first basis function (constant term)
            #              second column contains second basis function (linear term)
            # * cmatrix -  matrix of linear constraints:
            #                  [ 1.0  0.0  0.0 ]
            #              first two columns contain coefficients before basis functions,
            #              last column contains desired value of their sum.
            #              So [1,0,0] means "1*constant_term + 0*linear_term = 0" 
            #
            y = [0.072436,0.246944,0.491263,0.522300,0.714064,0.921929]
            if _spoil_scenario==0:
                spoil_vec_by_nan(y)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(y)
            fmatrix = [[1,0.0],[1,0.2],[1,0.4],[1,0.6],[1,0.8],[1,1.0]]
            if _spoil_scenario==5:
                spoil_mat_by_nan(fmatrix)
            if _spoil_scenario==6:
                spoil_mat_by_posinf(fmatrix)
            if _spoil_scenario==7:
                spoil_mat_by_neginf(fmatrix)
            if _spoil_scenario==8:
                spoil_mat_by_adding_row(fmatrix)
            if _spoil_scenario==9:
                spoil_mat_by_adding_col(fmatrix)
            if _spoil_scenario==10:
                spoil_mat_by_deleting_row(fmatrix)
            if _spoil_scenario==11:
                spoil_mat_by_deleting_col(fmatrix)
            cmatrix = [[1,0,0]]
            if _spoil_scenario==12:
                spoil_mat_by_nan(cmatrix)
            if _spoil_scenario==13:
                spoil_mat_by_posinf(cmatrix)
            if _spoil_scenario==14:
                spoil_mat_by_neginf(cmatrix)

            #
            # Constrained fitting without weights
            #
            info, c, rep = xalglib.lsfitlinearc(y, fmatrix, cmatrix)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(c, [0,0.932933], "real_vector", 0.0005)

            #
            # Constrained fitting with individual weights
            #
            w = [1, 1.414213, 1, 1, 1, 1]
            if _spoil_scenario==15:
                spoil_vec_by_nan(w)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==18:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==19:
                spoil_vec_by_deleting_element(w)
            info, c, rep = xalglib.lsfitlinearwc(y, w, fmatrix, cmatrix)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and doc_print_test(c, [0,0.938322], "real_vector", 0.0005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_linc                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_pol
    #      Unconstrained polynomial fitting
    #
    _TestResult = True
    for _spoil_scenario in range(-1,20):
        try:
            #
            # This example demonstrates polynomial fitting.
            #
            # Fitting is done by two (M=2) functions from polynomial basis:
            #     f0 = 1
            #     f1 = x
            # Basically, it just a linear fit; more complex polynomials may be used
            # (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
            # us to demonstrate polynomialfit() function in action.
            #
            # We have:
            # * x      set of abscissas
            # * y      experimental data
            #
            # Additionally we demonstrate weighted fitting, where second point has
            # more weight than other ones.
            #
            x = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            m = 2
            t = 2
            if _spoil_scenario==10:
                t = float("+inf")
            if _spoil_scenario==11:
                t = float("-inf")

            #
            # Fitting without individual weights
            #
            # NOTE: result is returned as barycentricinterpolant structure.
            #       if you want to get representation in the power basis,
            #       you can use barycentricbar2pow() function to convert
            #       from barycentric to power representation (see docs for 
            #       POLINT subpackage for more info).
            #
            info, p, rep = xalglib.polynomialfit(x, y, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.011, "real", 0.002)

            #
            # Fitting with individual weights
            #
            # NOTE: slightly different result is returned
            #
            w = [1,1.414213562,1,1,1,1,1,1,1,1,1]
            if _spoil_scenario==12:
                spoil_vec_by_nan(w)
            if _spoil_scenario==13:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==14:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==15:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==16:
                spoil_vec_by_deleting_element(w)
            xc = []
            if _spoil_scenario==17:
                spoil_vec_by_adding_element(xc)
            yc = []
            if _spoil_scenario==18:
                spoil_vec_by_adding_element(yc)
            dc = []
            if _spoil_scenario==19:
                spoil_vec_by_adding_element(dc)
            info, p, rep = xalglib.polynomialfitwc(x, y, w, xc, yc, dc, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.023, "real", 0.002)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_pol                      FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_polc
    #      Constrained polynomial fitting
    #
    _TestResult = True
    for _spoil_scenario in range(-1,29):
        try:
            #
            # This example demonstrates polynomial fitting.
            #
            # Fitting is done by two (M=2) functions from polynomial basis:
            #     f0 = 1
            #     f1 = x
            # with simple constraint on function value
            #     f(0) = 0
            # Basically, it just a linear fit; more complex polynomials may be used
            # (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
            # us to demonstrate polynomialfit() function in action.
            #
            # We have:
            # * x      set of abscissas
            # * y      experimental data
            # * xc     points where constraints are placed
            # * yc     constraints on derivatives
            # * dc     derivative indices
            #          (0 means function itself, 1 means first derivative)
            #
            x = [1.0,1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [0.9,1.1]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)
            w = [1,1]
            if _spoil_scenario==10:
                spoil_vec_by_nan(w)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==13:
                spoil_vec_by_adding_element(w)
            if _spoil_scenario==14:
                spoil_vec_by_deleting_element(w)
            xc = [0]
            if _spoil_scenario==15:
                spoil_vec_by_nan(xc)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(xc)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(xc)
            if _spoil_scenario==18:
                spoil_vec_by_adding_element(xc)
            if _spoil_scenario==19:
                spoil_vec_by_deleting_element(xc)
            yc = [0]
            if _spoil_scenario==20:
                spoil_vec_by_nan(yc)
            if _spoil_scenario==21:
                spoil_vec_by_posinf(yc)
            if _spoil_scenario==22:
                spoil_vec_by_neginf(yc)
            if _spoil_scenario==23:
                spoil_vec_by_adding_element(yc)
            if _spoil_scenario==24:
                spoil_vec_by_deleting_element(yc)
            dc = [0]
            if _spoil_scenario==25:
                spoil_vec_by_adding_element(dc)
            if _spoil_scenario==26:
                spoil_vec_by_deleting_element(dc)
            t = 2
            if _spoil_scenario==27:
                t = float("+inf")
            if _spoil_scenario==28:
                t = float("-inf")
            m = 2

            info, p, rep = xalglib.polynomialfitwc(x, y, w, xc, yc, dc, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.000, "real", 0.001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_polc                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_d_spline
    #      Unconstrained fitting by penalized regression spline
    #
    _TestResult = True
    for _spoil_scenario in range(-1,19):
        try:
            #
            # In this example we demonstrate penalized spline fitting of noisy data
            #
            # We have:
            # * x - abscissas
            # * y - vector of experimental data, straight line with small noise
            #
            x = [0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            y = [0.10,0.00,0.30,0.40,0.30,0.40,0.62,0.68,0.75,0.95]
            if _spoil_scenario==5:
                spoil_vec_by_nan(y)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(y)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(y)

            #
            # Fit with VERY small amount of smoothing (rho = -5.0)
            # and large number of basis functions (M=50).
            #
            # With such small regularization penalized spline almost fully reproduces function values
            #
            rho = -5.0
            if _spoil_scenario==10:
                rho = float("nan")
            if _spoil_scenario==11:
                rho = float("+inf")
            if _spoil_scenario==12:
                rho = float("-inf")
            info, s, rep = xalglib.spline1dfitpenalized(x, y, 50, rho)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            v = xalglib.spline1dcalc(s, 0.0)
            _TestResult = _TestResult and doc_print_test(v, 0.10, "real", 0.01)

            #
            # Fit with VERY large amount of smoothing (rho = 10.0)
            # and large number of basis functions (M=50).
            #
            # With such regularization our spline should become close to the straight line fit.
            # We will compare its value in x=1.0 with results obtained from such fit.
            #
            rho = +10.0
            if _spoil_scenario==13:
                rho = float("nan")
            if _spoil_scenario==14:
                rho = float("+inf")
            if _spoil_scenario==15:
                rho = float("-inf")
            info, s, rep = xalglib.spline1dfitpenalized(x, y, 50, rho)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            v = xalglib.spline1dcalc(s, 1.0)
            _TestResult = _TestResult and doc_print_test(v, 0.969, "real", 0.001)

            #
            # In real life applications you may need some moderate degree of fitting,
            # so we try to fit once more with rho=3.0.
            #
            rho = +3.0
            if _spoil_scenario==16:
                rho = float("nan")
            if _spoil_scenario==17:
                rho = float("+inf")
            if _spoil_scenario==18:
                rho = float("-inf")
            info, s, rep = xalglib.spline1dfitpenalized(x, y, 50, rho)
            _TestResult = _TestResult and doc_print_test(info, 1, "int")
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_d_spline                   FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_t_polfit_1
    #      Polynomial fitting, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,10):
        try:
            x = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x)
            y = [0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)
            m = 2
            t = 2
            if _spoil_scenario==8:
                t = float("+inf")
            if _spoil_scenario==9:
                t = float("-inf")
            info, p, rep = xalglib.polynomialfit(x, y, 11, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.011, "real", 0.002)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_t_polfit_1                 FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_t_polfit_2
    #      Polynomial fitting, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,14):
        try:
            x = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x)
            y = [0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)
            w = [1,1.414213562,1,1,1,1,1,1,1,1,1]
            if _spoil_scenario==8:
                spoil_vec_by_nan(w)
            if _spoil_scenario==9:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==10:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==11:
                spoil_vec_by_deleting_element(w)
            xc = []
            yc = []
            dc = []
            m = 2
            t = 2
            if _spoil_scenario==12:
                t = float("+inf")
            if _spoil_scenario==13:
                t = float("-inf")
            info, p, rep = xalglib.polynomialfitwc(x, y, w, 11, xc, yc, dc, 0, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.023, "real", 0.002)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_t_polfit_2                 FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST lsfit_t_polfit_3
    #      Polynomial fitting, full list of parameters.
    #
    _TestResult = True
    for _spoil_scenario in range(-1,23):
        try:
            x = [1.0,1.0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_deleting_element(x)
            y = [0.9,1.1]
            if _spoil_scenario==4:
                spoil_vec_by_nan(y)
            if _spoil_scenario==5:
                spoil_vec_by_posinf(y)
            if _spoil_scenario==6:
                spoil_vec_by_neginf(y)
            if _spoil_scenario==7:
                spoil_vec_by_deleting_element(y)
            w = [1,1]
            if _spoil_scenario==8:
                spoil_vec_by_nan(w)
            if _spoil_scenario==9:
                spoil_vec_by_posinf(w)
            if _spoil_scenario==10:
                spoil_vec_by_neginf(w)
            if _spoil_scenario==11:
                spoil_vec_by_deleting_element(w)
            xc = [0]
            if _spoil_scenario==12:
                spoil_vec_by_nan(xc)
            if _spoil_scenario==13:
                spoil_vec_by_posinf(xc)
            if _spoil_scenario==14:
                spoil_vec_by_neginf(xc)
            if _spoil_scenario==15:
                spoil_vec_by_deleting_element(xc)
            yc = [0]
            if _spoil_scenario==16:
                spoil_vec_by_nan(yc)
            if _spoil_scenario==17:
                spoil_vec_by_posinf(yc)
            if _spoil_scenario==18:
                spoil_vec_by_neginf(yc)
            if _spoil_scenario==19:
                spoil_vec_by_deleting_element(yc)
            dc = [0]
            if _spoil_scenario==20:
                spoil_vec_by_deleting_element(dc)
            m = 2
            t = 2
            if _spoil_scenario==21:
                t = float("+inf")
            if _spoil_scenario==22:
                t = float("-inf")
            info, p, rep = xalglib.polynomialfitwc(x, y, w, 2, xc, yc, dc, 1, m)
            v = xalglib.barycentriccalc(p, t)
            _TestResult = _TestResult and doc_print_test(v, 2.000, "real", 0.001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("lsfit_t_polfit_3                 FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_d_1
    #      Determinant calculation, real matrix, short form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            b = [[1,2],[2,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(b)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(b)
            a = xalglib.rmatrixdet(b)
            _TestResult = _TestResult and doc_print_test(a, -3, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_d_1                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_d_2
    #      Determinant calculation, real matrix, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            b = [[5,4],[4,5]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            a = xalglib.rmatrixdet(b, 2)
            _TestResult = _TestResult and doc_print_test(a, 9, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_d_2                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_d_3
    #      Determinant calculation, complex matrix, short form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            b = [[1+1j,2],[2,1-1j]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(b)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(b)
            a = xalglib.cmatrixdet(b)
            _TestResult = _TestResult and doc_print_test(a, -2, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_d_3                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_d_4
    #      Determinant calculation, complex matrix, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            b = [[5j,4],[4j,5]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            a = xalglib.cmatrixdet(b, 2)
            _TestResult = _TestResult and doc_print_test(a, 9j, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_d_4                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_d_5
    #      Determinant calculation, complex matrix with zero imaginary part, short form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,7):
        try:
            b = [[9,1],[2,1]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(b)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(b)
            a = xalglib.cmatrixdet(b)
            _TestResult = _TestResult and doc_print_test(a, 7, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_d_5                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_0
    #      Determinant calculation, real matrix, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            b = [[3,4],[-4,3]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            a = xalglib.rmatrixdet(b, 2)
            _TestResult = _TestResult and doc_print_test(a, 25, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_0                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_1
    #      Determinant calculation, real matrix, LU, short form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,9):
        try:
            b = [[1,2],[2,5]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(b)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(b)
            p = [1,1]
            if _spoil_scenario==7:
                spoil_vec_by_adding_element(p)
            if _spoil_scenario==8:
                spoil_vec_by_deleting_element(p)
            a = xalglib.rmatrixludet(b, p)
            _TestResult = _TestResult and doc_print_test(a, -5, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_1                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_2
    #      Determinant calculation, real matrix, LU, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            b = [[5,4],[4,5]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            p = [0,1]
            if _spoil_scenario==5:
                spoil_vec_by_deleting_element(p)
            a = xalglib.rmatrixludet(b, p, 2)
            _TestResult = _TestResult and doc_print_test(a, 25, "real", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_2                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_3
    #      Determinant calculation, complex matrix, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,5):
        try:
            b = [[5j,4],[-4,5j]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            a = xalglib.cmatrixdet(b, 2)
            _TestResult = _TestResult and doc_print_test(a, -9, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_3                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_4
    #      Determinant calculation, complex matrix, LU, short form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,9):
        try:
            b = [[1,2],[2,5j]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_adding_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_adding_col(b)
            if _spoil_scenario==5:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_col(b)
            p = [1,1]
            if _spoil_scenario==7:
                spoil_vec_by_adding_element(p)
            if _spoil_scenario==8:
                spoil_vec_by_deleting_element(p)
            a = xalglib.cmatrixludet(b, p)
            _TestResult = _TestResult and doc_print_test(a, -5j, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_4                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST matdet_t_5
    #      Determinant calculation, complex matrix, LU, full form
    #
    _TestResult = True
    for _spoil_scenario in range(-1,6):
        try:
            b = [[5,4j],[4,5]]
            if _spoil_scenario==0:
                spoil_mat_by_nan(b)
            if _spoil_scenario==1:
                spoil_mat_by_posinf(b)
            if _spoil_scenario==2:
                spoil_mat_by_neginf(b)
            if _spoil_scenario==3:
                spoil_mat_by_deleting_row(b)
            if _spoil_scenario==4:
                spoil_mat_by_deleting_col(b)
            p = [0,1]
            if _spoil_scenario==5:
                spoil_vec_by_deleting_element(p)
            a = xalglib.cmatrixludet(b, p, 2)
            _TestResult = _TestResult and doc_print_test(a, 25, "complex", 0.0001)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("matdet_t_5                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minasa_d_1
    #      Nonlinear optimization with bound constraints
    #
    _TestResult = True
    for _spoil_scenario in range(-1,24):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # subject to bound constraints -1<=x0<=+1, -1<=x1<=+1, using ASA.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            if _spoil_scenario==3:
                spoil_vec_by_adding_element(x)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(x)
            bndl = [-1,-1]
            if _spoil_scenario==5:
                spoil_vec_by_nan(bndl)
            if _spoil_scenario==6:
                spoil_vec_by_posinf(bndl)
            if _spoil_scenario==7:
                spoil_vec_by_neginf(bndl)
            if _spoil_scenario==8:
                spoil_vec_by_adding_element(bndl)
            if _spoil_scenario==9:
                spoil_vec_by_deleting_element(bndl)
            bndu = [+1,+1]
            if _spoil_scenario==10:
                spoil_vec_by_nan(bndu)
            if _spoil_scenario==11:
                spoil_vec_by_posinf(bndu)
            if _spoil_scenario==12:
                spoil_vec_by_neginf(bndu)
            if _spoil_scenario==13:
                spoil_vec_by_adding_element(bndu)
            if _spoil_scenario==14:
                spoil_vec_by_deleting_element(bndu)
            epsg = 0.0000000001
            if _spoil_scenario==15:
                epsg = float("nan")
            if _spoil_scenario==16:
                epsg = float("+inf")
            if _spoil_scenario==17:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==18:
                epsf = float("nan")
            if _spoil_scenario==19:
                epsf = float("+inf")
            if _spoil_scenario==20:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==21:
                epsx = float("nan")
            if _spoil_scenario==22:
                epsx = float("+inf")
            if _spoil_scenario==23:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.minasacreate(x, bndl, bndu)
            xalglib.minasasetcond(state, epsg, epsf, epsx, maxits)
            xalglib.minasaoptimize_g(state, function1_grad)
            x, rep = xalglib.minasaresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-1,1], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minasa_d_1                       FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST mincg_d_1
    #      Nonlinear optimization by CG
    #
    _TestResult = True
    for _spoil_scenario in range(-1,12):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # with nonlinear conjugate gradient method.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            maxits = 0

            state = xalglib.mincgcreate(x)
            xalglib.mincgsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.mincgoptimize_g(state, function1_grad)
            x, rep = xalglib.mincgresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("mincg_d_1                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST mincg_d_2
    #      Nonlinear optimization with additional settings and restarts
    #
    _TestResult = True
    for _spoil_scenario in range(-1,18):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # with nonlinear conjugate gradient method.
            #
            # Several advanced techniques are demonstrated:
            # * upper limit on step size
            # * restart from new point
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            epsg = 0.0000000001
            if _spoil_scenario==3:
                epsg = float("nan")
            if _spoil_scenario==4:
                epsg = float("+inf")
            if _spoil_scenario==5:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==6:
                epsf = float("nan")
            if _spoil_scenario==7:
                epsf = float("+inf")
            if _spoil_scenario==8:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==9:
                epsx = float("nan")
            if _spoil_scenario==10:
                epsx = float("+inf")
            if _spoil_scenario==11:
                epsx = float("-inf")
            stpmax = 0.1
            if _spoil_scenario==12:
                stpmax = float("nan")
            if _spoil_scenario==13:
                stpmax = float("+inf")
            if _spoil_scenario==14:
                stpmax = float("-inf")
            maxits = 0

            # first run
            state = xalglib.mincgcreate(x)
            xalglib.mincgsetcond(state, epsg, epsf, epsx, maxits)
            xalglib.mincgsetstpmax(state, stpmax)
            xalglib.mincgoptimize_g(state, function1_grad)
            x, rep = xalglib.mincgresults(state)

            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)

            # second run - algorithm is restarted with mincgrestartfrom()
            x = [10,10]
            if _spoil_scenario==15:
                spoil_vec_by_nan(x)
            if _spoil_scenario==16:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==17:
                spoil_vec_by_neginf(x)
            xalglib.mincgrestartfrom(state, x)
            xalglib.mincgoptimize_g(state, function1_grad)
            x, rep = xalglib.mincgresults(state)

            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-3,3], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("mincg_d_2                        FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minbleic_d_1
    #      Nonlinear optimization with bound constraints
    #
    _TestResult = True
    for _spoil_scenario in range(-1,25):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # subject to bound constraints -1<=x<=+1, -1<=y<=+1, using BLEIC optimizer.
            #
            x = [0,0]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            bndl = [-1,-1]
            if _spoil_scenario==3:
                spoil_vec_by_nan(bndl)
            if _spoil_scenario==4:
                spoil_vec_by_deleting_element(bndl)
            bndu = [+1,+1]
            if _spoil_scenario==5:
                spoil_vec_by_nan(bndu)
            if _spoil_scenario==6:
                spoil_vec_by_deleting_element(bndu)

            #
            # These variables define stopping conditions for the underlying CG algorithm.
            # They should be stringent enough in order to guarantee overall stability
            # of the outer iterations.
            #
            # We use very simple condition - |g|<=epsg
            #
            epsg = 0.000001
            if _spoil_scenario==7:
                epsg = float("nan")
            if _spoil_scenario==8:
                epsg = float("+inf")
            if _spoil_scenario==9:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==10:
                epsf = float("nan")
            if _spoil_scenario==11:
                epsf = float("+inf")
            if _spoil_scenario==12:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==13:
                epsx = float("nan")
            if _spoil_scenario==14:
                epsx = float("+inf")
            if _spoil_scenario==15:
                epsx = float("-inf")

            #
            # These variables define stopping conditions for the outer iterations:
            # * epso controls convergence of outer iterations; algorithm will stop
            #   when difference between solutions of subsequent unconstrained problems
            #   will be less than 0.0001
            # * epsi controls amount of infeasibility allowed in the final solution;
            #   algorithm will stop when infeasibility of the solution will be
            #   less than 0.0001
            #
            epso = 0.00001
            if _spoil_scenario==16:
                epso = float("nan")
            if _spoil_scenario==17:
                epso = float("+inf")
            if _spoil_scenario==18:
                epso = float("-inf")
            epsi = 0.00001
            if _spoil_scenario==19:
                epsi = float("nan")
            if _spoil_scenario==20:
                epsi = float("+inf")
            if _spoil_scenario==21:
                epsi = float("-inf")

            #
            # This variable defines barrier width.
            # 
            # Barrier width must be:
            # * small enough (below some problem-dependent value) in order for algorithm
            #   to  converge.  Necessary  condition  is that the target function must be
            #   well described by linear model in the areas as small as barrier width.
            # * not VERY small (in order to avoid  difficulties  associated  with  rapid
            #   changes in the modified function, ill-conditioning, round-off issues).
            #
            # Here we try small constant barrier width in the hope that this value will
            # be good enough. Another strategy (monotonically decreasing width) is possible,
            # and described in the other examples.
            #
            mu = 0.0001
            if _spoil_scenario==22:
                mu = float("nan")
            if _spoil_scenario==23:
                mu = float("+inf")
            if _spoil_scenario==24:
                mu = float("-inf")

            #
            # Now we are ready to actually optimize something:
            # * first we create optimizer
            # * we add boundary constraints
            # * we tune stopping conditions and barrier width
            # * and, finally, optimize and obtain results...
            #
            state = xalglib.minbleiccreate(x)
            xalglib.minbleicsetbc(state, bndl, bndu)
            xalglib.minbleicsetinnercond(state, epsg, epsf, epsx)
            xalglib.minbleicsetoutercond(state, epso, epsi)
            xalglib.minbleicsetbarrierwidth(state, mu)
            xalglib.minbleicoptimize_g(state, function1_grad)
            x, rep = xalglib.minbleicresults(state)

            #
            # ...and evaluate these results
            #
            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [-1,1], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minbleic_d_1                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    #
    # TEST minbleic_d_2
    #      Nonlinear optimization with linear inequality constraints
    #
    _TestResult = True
    for _spoil_scenario in range(-1,27):
        try:
            #
            # This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
            # subject to inequality constraints:
            # * x>=2 (posed as general linear constraint),
            # * x+y>=6
            # using BLEIC optimizer.
            #
            x = [5,5]
            if _spoil_scenario==0:
                spoil_vec_by_nan(x)
            if _spoil_scenario==1:
                spoil_vec_by_posinf(x)
            if _spoil_scenario==2:
                spoil_vec_by_neginf(x)
            c = [[1,0,2],[1,1,6]]
            if _spoil_scenario==3:
                spoil_mat_by_nan(c)
            if _spoil_scenario==4:
                spoil_mat_by_posinf(c)
            if _spoil_scenario==5:
                spoil_mat_by_neginf(c)
            if _spoil_scenario==6:
                spoil_mat_by_deleting_row(c)
            if _spoil_scenario==7:
                spoil_mat_by_deleting_col(c)
            ct = [1,1]
            if _spoil_scenario==8:
                spoil_vec_by_deleting_element(ct)

            #
            # These variables define stopping conditions for the underlying CG algorithm.
            # They should be stringent enough in order to guarantee overall stability
            # of the outer iterations.
            #
            # We use very simple condition - |g|<=epsg
            #
            epsg = 0.000001
            if _spoil_scenario==9:
                epsg = float("nan")
            if _spoil_scenario==10:
                epsg = float("+inf")
            if _spoil_scenario==11:
                epsg = float("-inf")
            epsf = 0
            if _spoil_scenario==12:
                epsf = float("nan")
            if _spoil_scenario==13:
                epsf = float("+inf")
            if _spoil_scenario==14:
                epsf = float("-inf")
            epsx = 0
            if _spoil_scenario==15:
                epsx = float("nan")
            if _spoil_scenario==16:
                epsx = float("+inf")
            if _spoil_scenario==17:
                epsx = float("-inf")

            #
            # These variables define stopping conditions for the outer iterations:
            # * epso controls convergence of outer iterations; algorithm will stop
            #   when difference between solutions of subsequent unconstrained problems
            #   will be less than 0.0001
            # * epsi controls amount of infeasibility allowed in the final solution;
            #   algorithm will stop when infeasibility of the solution will be
            #   less than 0.0001
            #
            epso = 0.00001
            if _spoil_scenario==18:
                epso = float("nan")
            if _spoil_scenario==19:
                epso = float("+inf")
            if _spoil_scenario==20:
                epso = float("-inf")
            epsi = 0.00001
            if _spoil_scenario==21:
                epsi = float("nan")
            if _spoil_scenario==22:
                epsi = float("+inf")
            if _spoil_scenario==23:
                epsi = float("-inf")

            #
            # This variable defines barrier width.
            # 
            # Barrier width must be:
            # * small enough (below some problem-dependent value) in order for algorithm
            #   to  converge.  Necessary  condition  is that the target function must be
            #   well described by linear model in the areas as small as barrier width.
            # * not VERY small (in order to avoid  difficulties  associated  with  rapid
            #   changes in the modified function, ill-conditioning, round-off issues).
            #
            # Here we try small constant barrier width in the hope that this value will
            # be good enough. Another strategy (monotonically decreasing width) is possible,
            # and described in the other examples.
            #
            mu = 0.0001
            if _spoil_scenario==24:
                mu = float("nan")
            if _spoil_scenario==25:
                mu = float("+inf")
            if _spoil_scenario==26:
                mu = float("-inf")

            #
            # Now we are ready to actually optimize something:
            # * first we create optimizer
            # * we add linear constraints
            # * we tune stopping conditions
            # * and, finally, optimize and obtain results...
            #
            state = xalglib.minbleiccreate(x)
            xalglib.minbleicsetlc(state, c, ct)
            xalglib.minbleicsetinnercond(state, epsg, epsf, epsx)
            xalglib.minbleicsetoutercond(state, epso, epsi)
            xalglib.minbleicsetbarrierwidth(state, mu)
            xalglib.minbleicoptimize_g(state, function1_grad)
            x, rep = xalglib.minbleicresults(state)

            #
            # ...and evaluate these results
            #
            _TestResult = _TestResult and doc_print_test(rep.terminationtype, 4, "int")
            _TestResult = _TestResult and doc_print_test(x, [2,4], "real_vector", 0.005)
            _TestResult = _TestResult and (_spoil_scenario==-1)
        except (RuntimeError, ValueError):
            _TestResult = _TestResult and (_spoil_scenario!=-1)
        except:
            raise
    if not _TestResult:
        sys.stdout.write("minbleic_d_2                     FAILED\n")
        sys.stdout.flush()
    _TotalResult = _TotalResult and _TestResult


    sys.stdout.write("78/78\n")
except:
    sys.stdout.write("Unhandled exception was raised!\n")
    sys.exit(1)
if _TotalResult:
    sys.exit(0)
else:
    sys.exit(1)
