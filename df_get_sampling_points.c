#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "mex.h"
#include "blas.h"
#include "lapack.h"

/* Eelke Visser, 11 Aug 2011 (eelke.visser@donders.ru.nl) */

void getSamplingPoints(
        double *xyz,
        unsigned int n,
        unsigned int nsl,
        double *TT,
        double *iPM,
        double *T2iMD2RPM,
        double *beta,
        double *txyz,
        mxLogical invert)
{
    char notrans = 'N';
    double zero = 0;
    double one = 1;
    ptrdiff_t pd1 = 1;
    ptrdiff_t pd4 = 4;
    
    for (unsigned int i = 0; i < n; i++, xyz += 4, txyz += 4)
    {
        double tmp[4];

        dgemv_(&notrans, &pd4, &pd4, &one, TT, &pd4, xyz, &pd1, &zero, tmp, &pd1);

        double T1[16] =
            {1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1};
        
        // NB: Everything is column-major!
        // (2, 1) -> 1 + 0 * 4 -> T1[1]
        // (2, 2) -> 1 + 1 * 4 -> T1[5]
        // (2, 4) -> 1 + 3 * 4 -> T1[13]

        int ti = tmp[2] - 1;
             
        if (ti <= 0)
        {
            T1[1] = beta[6];
            T1[5] = 1 + beta[6 + nsl];
            T1[13] = beta[6 + 2 * nsl];
        }
        else if (ti >= nsl - 1)
        {
            T1[1] = beta[6 + nsl - 1];
            T1[5] = 1 + beta[6 + 2 * nsl - 1];
            T1[13] = beta[6 + 3 * nsl - 1];
        }
        else
        {
            double ip;
            double fp = modf(ti, &ip);
            int ipi = (int)ip;
            
            T1[1] = (1 - fp) * beta[6 + ipi] + fp * beta[6 + ipi + 1];
            T1[5] = 1 + (1 - fp) * beta[6 + ipi + nsl] + fp * beta[6 + ipi + 1 + nsl];
            T1[13] = (1 - fp) * beta[6 + ipi + 2 * nsl] + fp * beta[6 + ipi + 1 + 2 * nsl];
        }

        /*
        double tmp2[4];
        dgemv_(&notrans, &pd4, &pd4, &one, T2iMD2RPM, &pd4, xyz, &pd1, &zero, tmp, &pd1);
        dgemv_(&notrans, &pd4, &pd4, &one, T1, &pd4, tmp, &pd1, &zero, tmp2, &pd1);
        dgemv_(&notrans, &pd4, &pd4, &one, iPM, &pd4, tmp2, &pd1, &zero, txyz, &pd1);
        */
        
        double tmpmat[16], final[16];
        dgemm_(&notrans, &notrans, &pd4, &pd4, &pd4, &one, T1, &pd4, T2iMD2RPM, &pd4, &zero, tmpmat, &pd4);
        dgemm_(&notrans, &notrans, &pd4, &pd4, &pd4, &one, iPM, &pd4, tmpmat, &pd4, &zero, final, &pd4);
        
        if (invert)
        {
            ptrdiff_t piv[4];
            ptrdiff_t info;
            memcpy(txyz, xyz, 4 * sizeof(double));
            dgesv_(&pd4, &pd1, final, &pd4, piv, txyz, &pd4, &info);
            if (info)
                mexWarnMsgTxt("df_get_sampling_points: Warning - inversion failed");
        } else
            dgemv_(&notrans, &pd4, &pd4, &one, final, &pd4, xyz, &pd1, &zero, txyz, &pd1);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*
     * [txyz] = df_get_sampling_points(xyz, TT, iPM, T2iMD2RPM, beta, invert)
     *
     * Everything is column-major in here, as both matlab and blas use that convention ..
     */
    
    if (nrhs != 6)
        mexErrMsgTxt("Wrong number of input arguments");
    
    if (mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetM(prhs[1]) != 4 || mxGetN(prhs[1]) != 4
            || mxGetNumberOfDimensions(prhs[2]) != 2 || mxGetM(prhs[2]) != 4 || mxGetN(prhs[2]) != 4
            || mxGetNumberOfDimensions(prhs[3]) != 2 || mxGetM(prhs[3]) != 4 || mxGetN(prhs[3]) != 4
            || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
            || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
        mexErrMsgTxt("Transformations should be 4 x 4 matrices of real doubles");
    
    if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetM(prhs[0]) != 4 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("xyz should be a 4 x N matrix of real doubles");
    
    if (mxGetNumberOfDimensions(prhs[4]) != 2 || mxGetM(prhs[4]) != 1 || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])
            || mxGetN(prhs[4]) % 3 || mxGetN(prhs[4]) < 9)
        mexErrMsgTxt("beta should be a row vector of real doubles and of length 3 N + 6, where N is the number of slices");
    
    if (mxGetNumberOfDimensions(prhs[5]) != 2 || mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 1 || !mxIsLogical(prhs[5]))
        mexErrMsgTxt("invert should be a logical");
    
    unsigned int n = mxGetN(prhs[0]);
    unsigned int nsl = mxGetN(prhs[4]) / 3 - 2;
    double *xyz = mxGetPr(prhs[0]);
    double *TT = mxGetPr(prhs[1]);
    double *iPM = mxGetPr(prhs[2]);
    double *T2iMD2RPM = mxGetPr(prhs[3]);
    double *beta = mxGetPr(prhs[4]);
    mxLogical invert = mxGetLogicals(prhs[5]);
    
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    plhs[0] = mxCreateNumericMatrix(4, n, mxDOUBLE_CLASS, mxREAL);
    double *txyz = (double*)mxGetPr(plhs[0]);
    
    getSamplingPoints(xyz, n, nsl, TT, iPM, T2iMD2RPM, beta, txyz, invert);
}
