

#include <stdafx.h>
#include <stdio.h>
#include "testnearestneighborunit.h"


/*$ Declarations $*/
static void testnearestneighborunit_unset2d(/* Complex */ ae_matrix* a,
     ae_state *_state);
static void testnearestneighborunit_unset1d(/* Real    */ ae_vector* a,
     ae_state *_state);
static ae_bool testnearestneighborunit_kdtresultsdifferent(/* Real    */ ae_matrix* refxy,
     ae_int_t ntotal,
     /* Real    */ ae_matrix* qx,
     /* Real    */ ae_matrix* qxy,
     /* Integer */ ae_vector* qt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state);
static double testnearestneighborunit_vnorm(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t normtype,
     ae_state *_state);
static void testnearestneighborunit_testkdtuniform(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_int_t normtype,
     ae_bool* kdterrors,
     ae_state *_state);


/*$ Body $*/


/*************************************************************************
Testing Nearest Neighbor Search
*************************************************************************/
ae_bool testnearestneighbor(ae_bool silent, ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xy;
    ae_int_t i;
    ae_int_t j;
    double v;
    ae_int_t normtype;
    ae_int_t nx;
    ae_int_t ny;
    ae_int_t n;
    ae_int_t smalln;
    ae_int_t largen;
    ae_int_t passcount;
    ae_int_t pass;
    ae_bool waserrors;
    ae_bool kdterrors;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);

    kdterrors = ae_false;
    passcount = 2;
    smalln = 256;
    largen = 2048;
    ny = 3;
    
    /*
     *
     */
    for(pass=1; pass<=passcount; pass++)
    {
        for(normtype=0; normtype<=2; normtype++)
        {
            for(nx=1; nx<=3; nx++)
            {
                
                /*
                 * Test in hypercube
                 */
                ae_matrix_set_length(&xy, largen, nx+ny, _state);
                for(i=0; i<=largen-1; i++)
                {
                    for(j=0; j<=nx+ny-1; j++)
                    {
                        xy.ptr.pp_double[i][j] = 10*ae_randomreal(_state)-5;
                    }
                }
                for(n=1; n<=10; n++)
                {
                    testnearestneighborunit_testkdtuniform(&xy, n, nx, ae_randominteger(ny+1, _state), normtype, &kdterrors, _state);
                }
                testnearestneighborunit_testkdtuniform(&xy, largen, nx, ae_randominteger(ny+1, _state), normtype, &kdterrors, _state);
                
                /*
                 * Test clustered (2*N points, pairs of equal points)
                 */
                ae_matrix_set_length(&xy, 2*smalln, nx+ny, _state);
                for(i=0; i<=smalln-1; i++)
                {
                    for(j=0; j<=nx+ny-1; j++)
                    {
                        xy.ptr.pp_double[2*i+0][j] = 10*ae_randomreal(_state)-5;
                        xy.ptr.pp_double[2*i+1][j] = xy.ptr.pp_double[2*i+0][j];
                    }
                }
                testnearestneighborunit_testkdtuniform(&xy, 2*smalln, nx, ae_randominteger(ny+1, _state), normtype, &kdterrors, _state);
                
                /*
                 * Test degenerate case: all points are same except for one
                 */
                ae_matrix_set_length(&xy, smalln, nx+ny, _state);
                v = ae_randomreal(_state);
                for(i=0; i<=smalln-2; i++)
                {
                    for(j=0; j<=nx+ny-1; j++)
                    {
                        xy.ptr.pp_double[i][j] = v;
                    }
                }
                for(j=0; j<=nx+ny-1; j++)
                {
                    xy.ptr.pp_double[smalln-1][j] = 10*ae_randomreal(_state)-5;
                }
                testnearestneighborunit_testkdtuniform(&xy, smalln, nx, ae_randominteger(ny+1, _state), normtype, &kdterrors, _state);
            }
        }
    }
    
    /*
     * report
     */
    waserrors = kdterrors;
    if( !silent )
    {
        printf("TESTING NEAREST NEIGHBOR SEARCH\n");
        printf("* KD TREES:                              ");
        if( !kdterrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Unsets 2D array.
*************************************************************************/
static void testnearestneighborunit_unset2d(/* Complex */ ae_matrix* a,
     ae_state *_state)
{


    ae_matrix_set_length(a, 0+1, 0+1, _state);
    a->ptr.pp_complex[0][0] = ae_complex_from_d(2*ae_randomreal(_state)-1);
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void testnearestneighborunit_unset1d(/* Real    */ ae_vector* a,
     ae_state *_state)
{


    ae_vector_set_length(a, 0+1, _state);
    a->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
}


/*************************************************************************
Compare results from different queries:
* X     just X-values
* XY    X-values and Y-values
* XT    X-values and tag values
*************************************************************************/
static ae_bool testnearestneighborunit_kdtresultsdifferent(/* Real    */ ae_matrix* refxy,
     ae_int_t ntotal,
     /* Real    */ ae_matrix* qx,
     /* Real    */ ae_matrix* qxy,
     /* Integer */ ae_vector* qt,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    result = ae_false;
    for(i=0; i<=n-1; i++)
    {
        if( qt->ptr.p_int[i]<0||qt->ptr.p_int[i]>=ntotal )
        {
            result = ae_true;
            return result;
        }
        for(j=0; j<=nx-1; j++)
        {
            result = result||ae_fp_neq(qx->ptr.pp_double[i][j],refxy->ptr.pp_double[qt->ptr.p_int[i]][j]);
            result = result||ae_fp_neq(qxy->ptr.pp_double[i][j],refxy->ptr.pp_double[qt->ptr.p_int[i]][j]);
        }
        for(j=0; j<=ny-1; j++)
        {
            result = result||ae_fp_neq(qxy->ptr.pp_double[i][nx+j],refxy->ptr.pp_double[qt->ptr.p_int[i]][nx+j]);
        }
    }
    return result;
}


/*************************************************************************
Returns norm
*************************************************************************/
static double testnearestneighborunit_vnorm(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t normtype,
     ae_state *_state)
{
    ae_int_t i;
    double result;


    result = ae_randomreal(_state);
    if( normtype==0 )
    {
        result = 0;
        for(i=0; i<=n-1; i++)
        {
            result = ae_maxreal(result, ae_fabs(x->ptr.p_double[i], _state), _state);
        }
        return result;
    }
    if( normtype==1 )
    {
        result = 0;
        for(i=0; i<=n-1; i++)
        {
            result = result+ae_fabs(x->ptr.p_double[i], _state);
        }
        return result;
    }
    if( normtype==2 )
    {
        result = 0;
        for(i=0; i<=n-1; i++)
        {
            result = result+ae_sqr(x->ptr.p_double[i], _state);
        }
        result = ae_sqrt(result, _state);
        return result;
    }
    return result;
}


/*************************************************************************
Testing Nearest Neighbor Search on uniformly distributed hypercube

NormType: 0, 1, 2
D: space dimension
N: points count
*************************************************************************/
static void testnearestneighborunit_testkdtuniform(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t nx,
     ae_int_t ny,
     ae_int_t normtype,
     ae_bool* kdterrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    double errtol;
    ae_vector tags;
    ae_vector ptx;
    ae_vector tmpx;
    ae_vector tmpb;
    kdtree treex;
    kdtree treexy;
    kdtree treext;
    ae_matrix qx;
    ae_matrix qxy;
    ae_vector qtags;
    ae_vector qr;
    ae_int_t kx;
    ae_int_t kxy;
    ae_int_t kt;
    double eps;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t task;
    ae_bool isequal;
    double r;
    ae_int_t q;
    ae_int_t qcount;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tags, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ptx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpb, 0, DT_BOOL, _state, ae_true);
    _kdtree_init(&treex, _state, ae_true);
    _kdtree_init(&treexy, _state, ae_true);
    _kdtree_init(&treext, _state, ae_true);
    ae_matrix_init(&qx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&qxy, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&qtags, 0, DT_INT, _state, ae_true);
    ae_vector_init(&qr, 0, DT_REAL, _state, ae_true);

    qcount = 10;
    
    /*
     * Tol - roundoff error tolerance (for '>=' comparisons)
     */
    errtol = 100000*ae_machineepsilon;
    
    /*
     * fill tags
     */
    ae_vector_set_length(&tags, n, _state);
    for(i=0; i<=n-1; i++)
    {
        tags.ptr.p_int[i] = i;
    }
    
    /*
     * build trees
     */
    kdtreebuild(xy, n, nx, 0, normtype, &treex, _state);
    kdtreebuild(xy, n, nx, ny, normtype, &treexy, _state);
    kdtreebuildtagged(xy, &tags, n, nx, 0, normtype, &treext, _state);
    
    /*
     * allocate arrays
     */
    ae_vector_set_length(&tmpx, nx, _state);
    ae_vector_set_length(&tmpb, n, _state);
    ae_matrix_set_length(&qx, n, nx, _state);
    ae_matrix_set_length(&qxy, n, nx+ny, _state);
    ae_vector_set_length(&qtags, n, _state);
    ae_vector_set_length(&qr, n, _state);
    ae_vector_set_length(&ptx, nx, _state);
    
    /*
     * test general K-NN queries (with self-matches):
     * * compare results from different trees (must be equal) and
     *   check that correct (value,tag) pairs are returned
     * * test results from XT tree - let R be radius of query result.
     *   then all points not in result must be not closer than R.
     */
    for(q=1; q<=qcount; q++)
    {
        
        /*
         * Select K: 1..N
         */
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            k = 1+ae_randominteger(n, _state);
        }
        else
        {
            k = 1;
        }
        
        /*
         * Select point (either one of the points, or random)
         */
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            i = ae_randominteger(n, _state);
            ae_v_move(&ptx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
        }
        else
        {
            for(i=0; i<=nx-1; i++)
            {
                ptx.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
        }
        
        /*
         * Test:
         * * consistency of results from different queries
         * * points in query are IN the R-sphere (or at the boundary),
         *   and points not in query are outside of the R-sphere (or at the boundary)
         * * distances are correct and are ordered
         */
        kx = kdtreequeryknn(&treex, &ptx, k, ae_true, _state);
        kxy = kdtreequeryknn(&treexy, &ptx, k, ae_true, _state);
        kt = kdtreequeryknn(&treext, &ptx, k, ae_true, _state);
        if( (kx!=k||kxy!=k)||kt!=k )
        {
            *kdterrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        kdtreequeryresultsxi(&treex, &qx, _state);
        kdtreequeryresultsxyi(&treexy, &qxy, _state);
        kdtreequeryresultstagsi(&treext, &qtags, _state);
        kdtreequeryresultsdistancesi(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
        kdtreequeryresultsx(&treex, &qx, _state);
        kdtreequeryresultsxy(&treexy, &qxy, _state);
        kdtreequeryresultstags(&treext, &qtags, _state);
        kdtreequeryresultsdistances(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
        for(i=0; i<=n-1; i++)
        {
            tmpb.ptr.p_bool[i] = ae_true;
        }
        r = 0;
        for(i=0; i<=k-1; i++)
        {
            tmpb.ptr.p_bool[qtags.ptr.p_int[i]] = ae_false;
            ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
            ae_v_sub(&tmpx.ptr.p_double[0], 1, &qx.ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
            r = ae_maxreal(r, testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state), _state);
        }
        for(i=0; i<=n-1; i++)
        {
            if( tmpb.ptr.p_bool[i] )
            {
                ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
                ae_v_sub(&tmpx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
                *kdterrors = *kdterrors||ae_fp_less(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state),r*(1-errtol));
            }
        }
        for(i=0; i<=k-2; i++)
        {
            *kdterrors = *kdterrors||ae_fp_greater(qr.ptr.p_double[i],qr.ptr.p_double[i+1]);
        }
        for(i=0; i<=k-1; i++)
        {
            ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
            ae_v_sub(&tmpx.ptr.p_double[0], 1, &xy->ptr.pp_double[qtags.ptr.p_int[i]][0], 1, ae_v_len(0,nx-1));
            *kdterrors = *kdterrors||ae_fp_greater(ae_fabs(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state)-qr.ptr.p_double[i], _state),errtol);
        }
        
        /*
         * Test reallocation properties: buffered functions must automatically
         * resize array which is too small, but leave unchanged array which is
         * too large.
         */
        if( n>=2 )
        {
            
            /*
             * First step: array is too small, two elements are required
             */
            k = 2;
            kx = kdtreequeryknn(&treex, &ptx, k, ae_true, _state);
            kxy = kdtreequeryknn(&treexy, &ptx, k, ae_true, _state);
            kt = kdtreequeryknn(&treext, &ptx, k, ae_true, _state);
            if( (kx!=k||kxy!=k)||kt!=k )
            {
                *kdterrors = ae_true;
                ae_frame_leave(_state);
                return;
            }
            ae_matrix_set_length(&qx, 1, 1, _state);
            ae_matrix_set_length(&qxy, 1, 1, _state);
            ae_vector_set_length(&qtags, 1, _state);
            ae_vector_set_length(&qr, 1, _state);
            kdtreequeryresultsx(&treex, &qx, _state);
            kdtreequeryresultsxy(&treexy, &qxy, _state);
            kdtreequeryresultstags(&treext, &qtags, _state);
            kdtreequeryresultsdistances(&treext, &qr, _state);
            *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
            
            /*
             * Second step: array is one row larger than needed, so only first
             * row is overwritten. Test it.
             */
            k = 1;
            kx = kdtreequeryknn(&treex, &ptx, k, ae_true, _state);
            kxy = kdtreequeryknn(&treexy, &ptx, k, ae_true, _state);
            kt = kdtreequeryknn(&treext, &ptx, k, ae_true, _state);
            if( (kx!=k||kxy!=k)||kt!=k )
            {
                *kdterrors = ae_true;
                ae_frame_leave(_state);
                return;
            }
            for(i=0; i<=nx-1; i++)
            {
                qx.ptr.pp_double[1][i] = _state->v_nan;
            }
            for(i=0; i<=nx+ny-1; i++)
            {
                qxy.ptr.pp_double[1][i] = _state->v_nan;
            }
            qtags.ptr.p_int[1] = 999;
            qr.ptr.p_double[1] = _state->v_nan;
            kdtreequeryresultsx(&treex, &qx, _state);
            kdtreequeryresultsxy(&treexy, &qxy, _state);
            kdtreequeryresultstags(&treext, &qtags, _state);
            kdtreequeryresultsdistances(&treext, &qr, _state);
            *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
            for(i=0; i<=nx-1; i++)
            {
                *kdterrors = *kdterrors||!ae_isnan(qx.ptr.pp_double[1][i], _state);
            }
            for(i=0; i<=nx+ny-1; i++)
            {
                *kdterrors = *kdterrors||!ae_isnan(qxy.ptr.pp_double[1][i], _state);
            }
            *kdterrors = *kdterrors||!(qtags.ptr.p_int[1]==999);
            *kdterrors = *kdterrors||!ae_isnan(qr.ptr.p_double[1], _state);
        }
        
        /*
         * Test reallocation properties: 'interactive' functions must allocate
         * new array on each call.
         */
        if( n>=2 )
        {
            
            /*
             * On input array is either too small or too large
             */
            for(k=1; k<=2; k++)
            {
                ae_assert(k==1||k==2, "KNN: internal error (unexpected K)!", _state);
                kx = kdtreequeryknn(&treex, &ptx, k, ae_true, _state);
                kxy = kdtreequeryknn(&treexy, &ptx, k, ae_true, _state);
                kt = kdtreequeryknn(&treext, &ptx, k, ae_true, _state);
                if( (kx!=k||kxy!=k)||kt!=k )
                {
                    *kdterrors = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                ae_matrix_set_length(&qx, 3-k, 3-k, _state);
                ae_matrix_set_length(&qxy, 3-k, 3-k, _state);
                ae_vector_set_length(&qtags, 3-k, _state);
                ae_vector_set_length(&qr, 3-k, _state);
                kdtreequeryresultsxi(&treex, &qx, _state);
                kdtreequeryresultsxyi(&treexy, &qxy, _state);
                kdtreequeryresultstagsi(&treext, &qtags, _state);
                kdtreequeryresultsdistancesi(&treext, &qr, _state);
                *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
                *kdterrors = (*kdterrors||qx.rows!=k)||qx.cols!=nx;
                *kdterrors = (*kdterrors||qxy.rows!=k)||qxy.cols!=nx+ny;
                *kdterrors = *kdterrors||qtags.cnt!=k;
                *kdterrors = *kdterrors||qr.cnt!=k;
            }
        }
    }
    
    /*
     * test general approximate K-NN queries (with self-matches):
     * * compare results from different trees (must be equal) and
     *   check that correct (value,tag) pairs are returned
     * * test results from XT tree - let R be radius of query result.
     *   then all points not in result must be not closer than R/(1+Eps).
     */
    for(q=1; q<=qcount; q++)
    {
        
        /*
         * Select K: 1..N
         */
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            k = 1+ae_randominteger(n, _state);
        }
        else
        {
            k = 1;
        }
        
        /*
         * Select Eps
         */
        eps = 0.5+ae_randomreal(_state);
        
        /*
         * Select point (either one of the points, or random)
         */
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            i = ae_randominteger(n, _state);
            ae_v_move(&ptx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
        }
        else
        {
            for(i=0; i<=nx-1; i++)
            {
                ptx.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
        }
        
        /*
         * Test:
         * * consistency of results from different queries
         * * points in query are IN the R-sphere (or at the boundary),
         *   and points not in query are outside of the R-sphere (or at the boundary)
         * * distances are correct and are ordered
         */
        kx = kdtreequeryaknn(&treex, &ptx, k, ae_true, eps, _state);
        kxy = kdtreequeryaknn(&treexy, &ptx, k, ae_true, eps, _state);
        kt = kdtreequeryaknn(&treext, &ptx, k, ae_true, eps, _state);
        if( (kx!=k||kxy!=k)||kt!=k )
        {
            *kdterrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        kdtreequeryresultsxi(&treex, &qx, _state);
        kdtreequeryresultsxyi(&treexy, &qxy, _state);
        kdtreequeryresultstagsi(&treext, &qtags, _state);
        kdtreequeryresultsdistancesi(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
        kdtreequeryresultsx(&treex, &qx, _state);
        kdtreequeryresultsxy(&treexy, &qxy, _state);
        kdtreequeryresultstags(&treext, &qtags, _state);
        kdtreequeryresultsdistances(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, k, nx, ny, _state);
        for(i=0; i<=n-1; i++)
        {
            tmpb.ptr.p_bool[i] = ae_true;
        }
        r = 0;
        for(i=0; i<=k-1; i++)
        {
            tmpb.ptr.p_bool[qtags.ptr.p_int[i]] = ae_false;
            ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
            ae_v_sub(&tmpx.ptr.p_double[0], 1, &qx.ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
            r = ae_maxreal(r, testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state), _state);
        }
        for(i=0; i<=n-1; i++)
        {
            if( tmpb.ptr.p_bool[i] )
            {
                ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
                ae_v_sub(&tmpx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
                *kdterrors = *kdterrors||ae_fp_less(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state),r*(1-errtol)/(1+eps));
            }
        }
        for(i=0; i<=k-2; i++)
        {
            *kdterrors = *kdterrors||ae_fp_greater(qr.ptr.p_double[i],qr.ptr.p_double[i+1]);
        }
        for(i=0; i<=k-1; i++)
        {
            ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
            ae_v_sub(&tmpx.ptr.p_double[0], 1, &xy->ptr.pp_double[qtags.ptr.p_int[i]][0], 1, ae_v_len(0,nx-1));
            *kdterrors = *kdterrors||ae_fp_greater(ae_fabs(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state)-qr.ptr.p_double[i], _state),errtol);
        }
    }
    
    /*
     * test general R-NN queries  (with self-matches):
     * * compare results from different trees (must be equal) and
     *   check that correct (value,tag) pairs are returned
     * * test results from XT tree - let R be radius of query result.
     *   then all points not in result must be not closer than R.
     */
    for(q=1; q<=qcount; q++)
    {
        
        /*
         * Select R
         */
        if( ae_fp_greater(ae_randomreal(_state),0.3) )
        {
            r = ae_maxreal(ae_randomreal(_state), ae_machineepsilon, _state);
        }
        else
        {
            r = ae_machineepsilon;
        }
        
        /*
         * Select point (either one of the points, or random)
         */
        if( ae_fp_greater(ae_randomreal(_state),0.5) )
        {
            i = ae_randominteger(n, _state);
            ae_v_move(&ptx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
        }
        else
        {
            for(i=0; i<=nx-1; i++)
            {
                ptx.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
        }
        
        /*
         * Test:
         * * consistency of results from different queries
         * * points in query are IN the R-sphere (or at the boundary),
         *   and points not in query are outside of the R-sphere (or at the boundary)
         * * distances are correct and are ordered
         */
        kx = kdtreequeryrnn(&treex, &ptx, r, ae_true, _state);
        kxy = kdtreequeryrnn(&treexy, &ptx, r, ae_true, _state);
        kt = kdtreequeryrnn(&treext, &ptx, r, ae_true, _state);
        if( kxy!=kx||kt!=kx )
        {
            *kdterrors = ae_true;
            ae_frame_leave(_state);
            return;
        }
        kdtreequeryresultsxi(&treex, &qx, _state);
        kdtreequeryresultsxyi(&treexy, &qxy, _state);
        kdtreequeryresultstagsi(&treext, &qtags, _state);
        kdtreequeryresultsdistancesi(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, kx, nx, ny, _state);
        kdtreequeryresultsx(&treex, &qx, _state);
        kdtreequeryresultsxy(&treexy, &qxy, _state);
        kdtreequeryresultstags(&treext, &qtags, _state);
        kdtreequeryresultsdistances(&treext, &qr, _state);
        *kdterrors = *kdterrors||testnearestneighborunit_kdtresultsdifferent(xy, n, &qx, &qxy, &qtags, kx, nx, ny, _state);
        for(i=0; i<=n-1; i++)
        {
            tmpb.ptr.p_bool[i] = ae_true;
        }
        for(i=0; i<=kx-1; i++)
        {
            tmpb.ptr.p_bool[qtags.ptr.p_int[i]] = ae_false;
        }
        for(i=0; i<=n-1; i++)
        {
            ae_v_move(&tmpx.ptr.p_double[0], 1, &ptx.ptr.p_double[0], 1, ae_v_len(0,nx-1));
            ae_v_sub(&tmpx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
            if( tmpb.ptr.p_bool[i] )
            {
                *kdterrors = *kdterrors||ae_fp_less(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state),r*(1-errtol));
            }
            else
            {
                *kdterrors = *kdterrors||ae_fp_greater(testnearestneighborunit_vnorm(&tmpx, nx, normtype, _state),r*(1+errtol));
            }
        }
        for(i=0; i<=kx-2; i++)
        {
            *kdterrors = *kdterrors||ae_fp_greater(qr.ptr.p_double[i],qr.ptr.p_double[i+1]);
        }
    }
    
    /*
     * Test self-matching:
     * * self-match - nearest neighbor of each point in XY is the point itself
     * * no self-match - nearest neighbor is NOT the point itself
     */
    if( n>1 )
    {
        
        /*
         * test for N=1 have non-general form, but it is not really needed
         */
        for(task=0; task<=1; task++)
        {
            for(i=0; i<=n-1; i++)
            {
                ae_v_move(&ptx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nx-1));
                kx = kdtreequeryknn(&treex, &ptx, 1, task==0, _state);
                kdtreequeryresultsxi(&treex, &qx, _state);
                if( kx!=1 )
                {
                    *kdterrors = ae_true;
                    ae_frame_leave(_state);
                    return;
                }
                isequal = ae_true;
                for(j=0; j<=nx-1; j++)
                {
                    isequal = isequal&&ae_fp_eq(qx.ptr.pp_double[0][j],ptx.ptr.p_double[j]);
                }
                if( task==0 )
                {
                    *kdterrors = *kdterrors||!isequal;
                }
                else
                {
                    *kdterrors = *kdterrors||isequal;
                }
            }
        }
    }
    ae_frame_leave(_state);
}


/*$ End $*/
