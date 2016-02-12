//
// A class to compute the intersection of two simplices
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "SimplexIntersect.hpp"
#include "ESMCI_Exception.h"

// Fortran
#define FC_GLOBAL(name,NAME) name##_

#ifdef FC_GLOBAL
#define _GESV_ FC_GLOBAL(dgesv,DGESV)
#else
#define _GESV_ dgesv_
#endif

extern "C" void _GESV_(int*, int*, double*, int*, int*, double*, int*, int*);

bool SimplexIntersect::solve(const Simplex& smplx1, const Simplex& smplx2) {

    // make sure the dimensions agree
    if (smplx1.getNumberOfDimensions() != mndims || 
        smplx2.getNumberOfDimensions() != mndims) {
        throw ESMCI::Ex() << "SimplexIntersect::solve mismatch in the number of space dimensions";
    }
    if (smplx1.getNumberOfParametricDimensions() + 
        smplx2.getNumberOfParametricDimensions() != mndims) {
        throw ESMCI::Ex() << 
          "SimplexIntersect::solve the two simplices must be complementary (e.g. face and edge in 3d)";
    }

    //
    // build matrix system
    //

    // first simplex
    const std::vector< std::vector<double> >& p1s = smplx1.getPoints();
    mndims1 = p1s.size() - 1;
    for (size_t i = 0; i < mndims; ++i) {
        for (size_t j = 0; j < mndims1; ++j) {
            size_t k = i + mndims*j;
            // slightly perturb the vertex to avoid degeneracy issues
            mmat[k] = p1s[j+1][i] - p1s[0][i] + (i + 1)*moffset;
        }
        msol[i] = -p1s[0][i];
    }
    // second simplex
    const std::vector< std::vector<double> >& p2s = smplx2.getPoints();
    for (size_t i = 0; i < mndims; ++i) {
        for (size_t j = 0; j < mndims - mndims1; ++j) {
            size_t k = i + mndims*(j + mndims1);
            mmat[k] = p2s[0][i] - p2s[j+1][i];
        }
        msol[i] += p2s[0][i];
    }

    int info;
    int m = (int) mndims;
    int nrhs = 1;
    std::vector<int> ipiv(m);

    // solve
    _GESV_(&m, &nrhs, &mmat.front(), &m, &ipiv.front(), &msol.front(), &m, &info);

    return (info == 0);
}   
    
bool SimplexIntersect::intersect(double tol) {

    bool inSmplx1 = true;
    double sum1 = 0;
    for (size_t i =0; i < mndims1; ++i) {
        double val = msol[i];
        sum1 += val;
        inSmplx1 &= (val >= 0.0 - tol) && (val <= 1.0 + tol);
    }
    bool inSmplx2 = true;
    double sum2 = 0;
    for (size_t i =0; i < mndims - mndims1; ++i) {
        double val = msol[i + mndims1];
        sum2 += val;
        inSmplx2 &= (val >= 0.0 - tol) && (val <= 1.0 + tol);
    }

    return inSmplx1 && inSmplx2 && sum1 <= 1.0 + tol && sum2 <= 1.0 + tol;
}
