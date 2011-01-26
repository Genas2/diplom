#include <stdio.h>
#include <time.h>
#include "hqrnd.h"
#include "tsort.h"
#include "apserv.h"
#include "nearestneighbor.h"
#include "ablasf.h"
#include "ablas.h"
#include "basicstatops.h"
#include "basestat.h"
#include "bdss.h"
#include "dforest.h"
#include "blas.h"
#include "kmeans.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "lda.h"
#include "gammafunc.h"
#include "normaldistr.h"
#include "igammaf.h"
#include "bdsvd.h"
#include "svd.h"
#include "linreg.h"
#include "mlpbase.h"
#include "xblas.h"
#include "densesolver.h"
#include "logit.h"
#include "fbls.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "mlptrain.h"
#include "mlpe.h"
#include "pca.h"
#include "odesolver.h"
#include "ftbase.h"
#include "fft.h"
#include "conv.h"
#include "corr.h"
#include "fht.h"
#include "gq.h"
#include "gkq.h"
#include "autogk.h"
#include "idwint.h"
#include "ratint.h"
#include "polint.h"
#include "spline1d.h"
#include "minlm.h"
#include "lsfit.h"
#include "pspline.h"
#include "spline2d.h"
#include "matdet.h"
#include "spdgevd.h"
#include "inverseupdate.h"
#include "schur.h"
#include "minasa.h"
#include "mincg.h"
#include "minbleic.h"
#include "nleq.h"
#include "airyf.h"
#include "bessel.h"
#include "betaf.h"
#include "ibetaf.h"
#include "nearunityunit.h"
#include "binomialdistr.h"
#include "chebyshev.h"
#include "chisquaredistr.h"
#include "dawson.h"
#include "elliptic.h"
#include "expintegrals.h"
#include "fdistr.h"
#include "fresnel.h"
#include "hermite.h"
#include "jacobianelliptic.h"
#include "laguerre.h"
#include "legendre.h"
#include "poissondistr.h"
#include "psif.h"
#include "studenttdistr.h"
#include "trigintegrals.h"
#include "correlationtests.h"
#include "jarquebera.h"
#include "mannwhitneyu.h"
#include "stest.h"
#include "studentttests.h"
#include "variancetests.h"
#include "wsr.h"
#include "alglibbasics.h"

#include "testhqrndunit.h"
#include "testtsortunit.h"
#include "testnearestneighborunit.h"
#include "testablasunit.h"
#include "testbasestatunit.h"
#include "testbdssunit.h"
#include "testdforestunit.h"
#include "testblasunit.h"
#include "testkmeansunit.h"
#include "testhblasunit.h"
#include "testreflectionsunit.h"
#include "testcreflectionsunit.h"
#include "testsblasunit.h"
#include "testortfacunit.h"
#include "testevdunit.h"
#include "testmatgenunit.h"
#include "testtrfacunit.h"
#include "testtrlinsolveunit.h"
#include "testsafesolveunit.h"
#include "testrcondunit.h"
#include "testmatinvunit.h"
#include "testldaunit.h"
#include "testgammafuncunit.h"
#include "testbdsvdunit.h"
#include "testsvdunit.h"
#include "testlinregunit.h"
#include "testxblasunit.h"
#include "testdensesolverunit.h"
#include "testfblsunit.h"
#include "testlinminunit.h"
#include "testminlbfgsunit.h"
#include "testmlptrainunit.h"
#include "testmlpeunit.h"
#include "testpcaunit.h"
#include "testodesolverunit.h"
#include "testfftunit.h"
#include "testconvunit.h"
#include "testcorrunit.h"
#include "testfhtunit.h"
#include "testgqunit.h"
#include "testgkqunit.h"
#include "testautogkunit.h"
#include "testidwintunit.h"
#include "testratintunit.h"
#include "testpolintunit.h"
#include "testspline1dunit.h"
#include "testminlmunit.h"
#include "testlsfitunit.h"
#include "testpsplineunit.h"
#include "testspline2dunit.h"
#include "testspdgevdunit.h"
#include "testinverseupdateunit.h"
#include "testschurunit.h"
#include "testminasaunit.h"
#include "testmincgunit.h"
#include "testminbleicunit.h"
#include "testnlequnit.h"
#include "testchebyshevunit.h"
#include "testhermiteunit.h"
#include "testlaguerreunit.h"
#include "testlegendreunit.h"
#include "testalglibbasicsunit.h"


/*$ Body $*/
ae_bool call_unittest(unsigned seed, ae_bool(*testfunc)(ae_bool, ae_state*), int *psticky)
{
#ifndef AE_USE_CPP_ERROR_HANDLING
    ae_state _alglib_env_state;
    ae_frame _frame_block;
    jmp_buf _break_jump;
    ae_bool result;
    
    ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
        *psticky = 1;
        return ae_false;
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    ae_frame_make(&_alglib_env_state, &_frame_block);
    srand(seed);
    result = testfunc(ae_true, &_alglib_env_state);
    ae_state_clear(&_alglib_env_state);
    if( !result )
        *psticky = 1;
    return result;
#else
    try
    {
        ae_state _alglib_env_state;
        ae_frame _frame_block;
        ae_bool result;
    
        ae_state_init(&_alglib_env_state);
        ae_frame_make(&_alglib_env_state, &_frame_block);
        srand(seed);
        result = testfunc(ae_true, &_alglib_env_state);
        ae_state_clear(&_alglib_env_state);
        if( !result )
            *psticky = 1;
        return result;
    }
    catch(...)
    {
        *psticky = 1;
        return ae_false;
    }
#endif
}

int main(int argc, char **argv)
{
    unsigned seed;
    int result;
    double a;
    ae_int32_t *p;
    if( argc==2 )
        seed = (unsigned)atoi(argv[1]);
    else
    {
        time_t t;
        seed = (unsigned)time(&t);
    }
    
    /*
     * TODO: report results from automatic CPU detection
     */
    
    /*
     * Architecture
     */
    if( sizeof(void*)==4 )
        printf("HARDWARE: 32-bit\n");
    else if( sizeof(void*)==8 )
        printf("HARDWARE: 64-bit\n");
    else
        printf("HARDWARE: strange (non-32, non-64)\n");
    
    /* 
     * determine endianness of hardware.
     * 1983 is a good number - non-periodic double representation allow us to
     * easily distinguish between upper and lower halfs and to detect mixed endian hardware.
     */
    a = 1.0/1983.0; 
    p = (ae_int32_t*)(&a);
    if( p[1]==0x3f408642 )
        printf("HARDWARE: little-endian\n");
    else if( p[0]==0x3f408642 )
        printf("HARDWARE: big-endian\n");
    else
        printf("HARDWARE: mixed-endian\n");

    /*
     * now we are ready to test!
     */
    result = 0;
    if( call_unittest(seed, testhqrnd, &result) )
        printf("%-32s OK\n", "hqrnd");
    else
        printf("%-32s FAILED(seed=%u)\n", "hqrnd", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testtsort, &result) )
        printf("%-32s OK\n", "tsort");
    else
        printf("%-32s FAILED(seed=%u)\n", "tsort", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testnearestneighbor, &result) )
        printf("%-32s OK\n", "nearestneighbor");
    else
        printf("%-32s FAILED(seed=%u)\n", "nearestneighbor", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testablas, &result) )
        printf("%-32s OK\n", "ablas");
    else
        printf("%-32s FAILED(seed=%u)\n", "ablas", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testbasestat, &result) )
        printf("%-32s OK\n", "basestat");
    else
        printf("%-32s FAILED(seed=%u)\n", "basestat", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testbdss, &result) )
        printf("%-32s OK\n", "bdss");
    else
        printf("%-32s FAILED(seed=%u)\n", "bdss", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testdforest, &result) )
        printf("%-32s OK\n", "dforest");
    else
        printf("%-32s FAILED(seed=%u)\n", "dforest", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testblas, &result) )
        printf("%-32s OK\n", "blas");
    else
        printf("%-32s FAILED(seed=%u)\n", "blas", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testkmeans, &result) )
        printf("%-32s OK\n", "kmeans");
    else
        printf("%-32s FAILED(seed=%u)\n", "kmeans", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testhblas, &result) )
        printf("%-32s OK\n", "hblas");
    else
        printf("%-32s FAILED(seed=%u)\n", "hblas", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testreflections, &result) )
        printf("%-32s OK\n", "reflections");
    else
        printf("%-32s FAILED(seed=%u)\n", "reflections", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testcreflections, &result) )
        printf("%-32s OK\n", "creflections");
    else
        printf("%-32s FAILED(seed=%u)\n", "creflections", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testsblas, &result) )
        printf("%-32s OK\n", "sblas");
    else
        printf("%-32s FAILED(seed=%u)\n", "sblas", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testortfac, &result) )
        printf("%-32s OK\n", "ortfac");
    else
        printf("%-32s FAILED(seed=%u)\n", "ortfac", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testevd, &result) )
        printf("%-32s OK\n", "evd");
    else
        printf("%-32s FAILED(seed=%u)\n", "evd", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testmatgen, &result) )
        printf("%-32s OK\n", "matgen");
    else
        printf("%-32s FAILED(seed=%u)\n", "matgen", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testtrfac, &result) )
        printf("%-32s OK\n", "trfac");
    else
        printf("%-32s FAILED(seed=%u)\n", "trfac", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testtrlinsolve, &result) )
        printf("%-32s OK\n", "trlinsolve");
    else
        printf("%-32s FAILED(seed=%u)\n", "trlinsolve", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testsafesolve, &result) )
        printf("%-32s OK\n", "safesolve");
    else
        printf("%-32s FAILED(seed=%u)\n", "safesolve", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testrcond, &result) )
        printf("%-32s OK\n", "rcond");
    else
        printf("%-32s FAILED(seed=%u)\n", "rcond", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testmatinv, &result) )
        printf("%-32s OK\n", "matinv");
    else
        printf("%-32s FAILED(seed=%u)\n", "matinv", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlda, &result) )
        printf("%-32s OK\n", "lda");
    else
        printf("%-32s FAILED(seed=%u)\n", "lda", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testgammafunc, &result) )
        printf("%-32s OK\n", "gammafunc");
    else
        printf("%-32s FAILED(seed=%u)\n", "gammafunc", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testbdsvd, &result) )
        printf("%-32s OK\n", "bdsvd");
    else
        printf("%-32s FAILED(seed=%u)\n", "bdsvd", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testsvd, &result) )
        printf("%-32s OK\n", "svd");
    else
        printf("%-32s FAILED(seed=%u)\n", "svd", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlinreg, &result) )
        printf("%-32s OK\n", "linreg");
    else
        printf("%-32s FAILED(seed=%u)\n", "linreg", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testxblas, &result) )
        printf("%-32s OK\n", "xblas");
    else
        printf("%-32s FAILED(seed=%u)\n", "xblas", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testdensesolver, &result) )
        printf("%-32s OK\n", "densesolver");
    else
        printf("%-32s FAILED(seed=%u)\n", "densesolver", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testfbls, &result) )
        printf("%-32s OK\n", "fbls");
    else
        printf("%-32s FAILED(seed=%u)\n", "fbls", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlinmin, &result) )
        printf("%-32s OK\n", "linmin");
    else
        printf("%-32s FAILED(seed=%u)\n", "linmin", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testminlbfgs, &result) )
        printf("%-32s OK\n", "minlbfgs");
    else
        printf("%-32s FAILED(seed=%u)\n", "minlbfgs", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testmlptrain, &result) )
        printf("%-32s OK\n", "mlptrain");
    else
        printf("%-32s FAILED(seed=%u)\n", "mlptrain", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testmlpe, &result) )
        printf("%-32s OK\n", "mlpe");
    else
        printf("%-32s FAILED(seed=%u)\n", "mlpe", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testpca, &result) )
        printf("%-32s OK\n", "pca");
    else
        printf("%-32s FAILED(seed=%u)\n", "pca", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testodesolver, &result) )
        printf("%-32s OK\n", "odesolver");
    else
        printf("%-32s FAILED(seed=%u)\n", "odesolver", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testfft, &result) )
        printf("%-32s OK\n", "fft");
    else
        printf("%-32s FAILED(seed=%u)\n", "fft", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testconv, &result) )
        printf("%-32s OK\n", "conv");
    else
        printf("%-32s FAILED(seed=%u)\n", "conv", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testcorr, &result) )
        printf("%-32s OK\n", "corr");
    else
        printf("%-32s FAILED(seed=%u)\n", "corr", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testfht, &result) )
        printf("%-32s OK\n", "fht");
    else
        printf("%-32s FAILED(seed=%u)\n", "fht", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testgq, &result) )
        printf("%-32s OK\n", "gq");
    else
        printf("%-32s FAILED(seed=%u)\n", "gq", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testgkq, &result) )
        printf("%-32s OK\n", "gkq");
    else
        printf("%-32s FAILED(seed=%u)\n", "gkq", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testautogk, &result) )
        printf("%-32s OK\n", "autogk");
    else
        printf("%-32s FAILED(seed=%u)\n", "autogk", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testidwint, &result) )
        printf("%-32s OK\n", "idwint");
    else
        printf("%-32s FAILED(seed=%u)\n", "idwint", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testratint, &result) )
        printf("%-32s OK\n", "ratint");
    else
        printf("%-32s FAILED(seed=%u)\n", "ratint", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testpolint, &result) )
        printf("%-32s OK\n", "polint");
    else
        printf("%-32s FAILED(seed=%u)\n", "polint", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testspline1d, &result) )
        printf("%-32s OK\n", "spline1d");
    else
        printf("%-32s FAILED(seed=%u)\n", "spline1d", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testminlm, &result) )
        printf("%-32s OK\n", "minlm");
    else
        printf("%-32s FAILED(seed=%u)\n", "minlm", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlsfit, &result) )
        printf("%-32s OK\n", "lsfit");
    else
        printf("%-32s FAILED(seed=%u)\n", "lsfit", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testpspline, &result) )
        printf("%-32s OK\n", "pspline");
    else
        printf("%-32s FAILED(seed=%u)\n", "pspline", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testspline2d, &result) )
        printf("%-32s OK\n", "spline2d");
    else
        printf("%-32s FAILED(seed=%u)\n", "spline2d", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testspdgevd, &result) )
        printf("%-32s OK\n", "spdgevd");
    else
        printf("%-32s FAILED(seed=%u)\n", "spdgevd", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testinverseupdate, &result) )
        printf("%-32s OK\n", "inverseupdate");
    else
        printf("%-32s FAILED(seed=%u)\n", "inverseupdate", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testschur, &result) )
        printf("%-32s OK\n", "schur");
    else
        printf("%-32s FAILED(seed=%u)\n", "schur", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testminasa, &result) )
        printf("%-32s OK\n", "minasa");
    else
        printf("%-32s FAILED(seed=%u)\n", "minasa", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testmincg, &result) )
        printf("%-32s OK\n", "mincg");
    else
        printf("%-32s FAILED(seed=%u)\n", "mincg", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testminbleic, &result) )
        printf("%-32s OK\n", "minbleic");
    else
        printf("%-32s FAILED(seed=%u)\n", "minbleic", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testnleq, &result) )
        printf("%-32s OK\n", "nleq");
    else
        printf("%-32s FAILED(seed=%u)\n", "nleq", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testchebyshev, &result) )
        printf("%-32s OK\n", "chebyshev");
    else
        printf("%-32s FAILED(seed=%u)\n", "chebyshev", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testhermite, &result) )
        printf("%-32s OK\n", "hermite");
    else
        printf("%-32s FAILED(seed=%u)\n", "hermite", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlaguerre, &result) )
        printf("%-32s OK\n", "laguerre");
    else
        printf("%-32s FAILED(seed=%u)\n", "laguerre", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testlegendre, &result) )
        printf("%-32s OK\n", "legendre");
    else
        printf("%-32s FAILED(seed=%u)\n", "legendre", (unsigned int)seed);
    fflush(stdout);
    if( call_unittest(seed, testalglibbasics, &result) )
        printf("%-32s OK\n", "alglibbasics");
    else
        printf("%-32s FAILED(seed=%u)\n", "alglibbasics", (unsigned int)seed);
    fflush(stdout);

    return result;
}
/*$ End $*/
