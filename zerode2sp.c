#include "mex.h"
#include <math.h>
#include <stdlib.h>

void erode(float *A, float *B, size_t nx, size_t ny, char *w,
        size_t nwx, size_t nwy)
{
    // nwx and nwy are assumed to be odd!
    size_t hnwx = nwx / 2;
    size_t hnwy = nwy / 2;
    
    for (size_t j = 0; j < ny; j++)
    {
        for (size_t i = 0; i < nx; i++)
        {
            size_t wi1 = i < hnwx ? hnwx - i : 0;
            size_t wi2 = i >= nx - hnwx ? hnwx - i + nx - 1 : nwx - 1;
            size_t wj1 = j < hnwy ? hnwy - j : 0;
            size_t wj2 = j >= ny - hnwy ? hnwy - j + ny - 1 : nwy - 1;
            
            float min = INFINITY;
            
            for (size_t wj = wj1; wj <= wj2; wj++)
            {
                for (size_t wi = wi1; wi <= wi2; wi++)
                {
                    if (w[wj * nwx + wi] > 0)
                    {
                        float val = A[(j - hnwy + wj) * nx + i - hnwx + wi];
                        if (val < min)
                            min = val;
                    }
                }
            }
            
            *B++ = min;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* arguments:
     *
     * img - image (m-by-n-by-k float, i.e. float image[k][n][m])
     * w - window (m-by-n char, i.e. char w[n][m])
     */
    
    if (nrhs != 2)
        mexErrMsgTxt("Wrong number of arguments (expected: image, window).");

    if (!mxIsSingle(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 3 ||
            mxIsComplex(prhs[0]))
        mexErrMsgTxt("Image should be a 3D array of real singles.");
    
    if  (!mxIsInt8(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 ||
            mxIsComplex(prhs[1]))
        mexErrMsgTxt("Window should be a 2D array of real int8s.");
    
    const mwSize *imgdims = mxGetDimensions(prhs[0]);
    const mwSize imgx = imgdims[0];
    const mwSize imgy = imgdims[1];
    const mwSize imgz = imgdims[2];
    float *img = (float*)mxGetPr(prhs[0]);
    
    size_t wx = mxGetM(prhs[1]);
    size_t wy = mxGetN(prhs[1]);
    char *w = (char*)mxGetPr(prhs[1]);
    
    if (!(wx % 2) || !(wy % 2))
        mexErrMsgTxt("Window size should be odd-by-odd.");

    plhs[0] = mxCreateNumericArray(3, imgdims, mxSINGLE_CLASS, mxREAL);
    float *newimg = (float*)mxGetPr(plhs[0]);

    #pragma omp parallel for schedule(static) num_threads(4)
    for (size_t k = 0; k < imgz; k++)
    {
        size_t offs = imgx * imgy * k;
        erode(img + offs, newimg + offs, imgx, imgy, w, wx, wy);
    }
}
