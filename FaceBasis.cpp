// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "FaceBasis.hpp"
#include "ESMCI_Exception.h"

void FaceBasis::computeWeights(const Simplex& smplx) {

    this->mweights.clear();

    if (smplx.getNumberOfParametricDimensions() != 2) {
        throw ESMCI::Ex() << "FaceBasis::computeWeights can only compute weights for a 2-simplex (triangle)";
    }

    // the coordinates of the triangle we wish to interpolate over
    const std::vector<double>& xiA = smplx.getPoint(0); // target is a triangle
    const std::vector<double>& xiB = smplx.getPoint(1);
    const std::vector<double>& xiC = smplx.getPoint(2);

    // sigmas[0] = 1 - xiA[0] - xiA[1] - ...
    // dXiB = length in parameter space between points B and A
    // dXiC = length in parameter space between points C and A
    // sigmas[1] = - sum of dXiB's components
    // sigmas[2] = - sum of dXiC's components
    std::vector<double> sigmas(3, 0.0);
    std::vector<double> dXiB(xiB);
    std::vector<double> dXiC(xiC);
    sigmas[0] = 1.0;
    for (size_t i = 0; i < this->mnpdims; ++i) {
        dXiB[i] -= xiA[i];
        dXiC[i] -= xiA[i];
        sigmas[0] -= xiA[i];
        sigmas[1] -= dXiB[i];
        sigmas[2] -= dXiC[i];
    }

    // compute the d phi_j ^ d phi_i determinants for all combinations
    // of j and i = j+ 1...npdims. Because the basis functions are linear
    // d phi_j ^ d phi_i is constant inside the cell, this factor can be 
    // taken out of the integrals.
    double det;
    std::map< std::vector<size_t>, double > dets;

    // j and i
    std::vector<size_t> surfInds(2);

    // special case for j = 0
    surfInds[0] = 0; 
    for (size_t i = 1; i < this->mnpdims + 1; ++i) {
        det = sigmas[1]*dXiC[i-1] - sigmas[2]*dXiB[i-1];
        surfInds[1] = i;
        dets.insert( std::pair<std::vector<size_t>, double>(surfInds, det) );
    }

    for (size_t j = 1; j < this->mnpdims; ++j) {
        surfInds[0] = j;
        for (size_t i = j + 1; i < this->mnpdims + 1; ++i) {
            det = dXiB[j-1]*dXiC[i-1] - dXiC[j-1]*dXiB[i-1];
            surfInds[1] = i;
            dets.insert( std::pair<std::vector<size_t>, double>(surfInds, det) );
        }
    }

    // integrated basis functions
    std::vector<double> phiIntegrals(this->mnpdims + 1);
    phiIntegrals[0] = sigmas[0]/2.0 + sigmas[1]/6.0 + sigmas[2]/6.0;
    for (size_t i = 1; i < this->mnpdims + 1; ++i) {
        size_t im1 = i - 1;
        phiIntegrals[i] = xiA[im1]/2.0 + dXiB[im1]/6.0 + dXiC[im1]/6.0;
    }

    // iterate over all the faces
    this->msmplxIter.begin();
    bool go = true;
    while (go) {

        // the set of indices that uniquely identify the face
        const std::vector<size_t>& inds = this->msmplxIter();

        size_t i = inds[0];
        size_t j = inds[1];
        size_t k = inds[2]; 
        std::vector<size_t> surfInds(2);
        double integral = 0.0;
        double sign = 1.0;

        // assume i, j, k to be in ascending order
        surfInds[0] = j; 
        surfInds[1] = k;
        integral += 2.0*phiIntegrals[i]*dets[surfInds];

        // term: j, k, i
        sign = 1.0;
        surfInds[0] = k;
        surfInds[1] = i;
        if (k > i) {
            // swap
            surfInds[0] = i;
            surfInds[1] = k;
            sign = -1.0;
        }
        integral += 2.0*sign*phiIntegrals[j]*dets[surfInds];

        // term k, i, j
        sign = 1.0;
        surfInds[0] = i;
        surfInds[1] = j;
        if (i > j) {
            // swap 
            surfInds[0] = j;
            surfInds[1] = i;
            sign = -1.0;
        }
        integral += 2.0*sign*phiIntegrals[k]*dets[surfInds];

        std::pair< std::vector<size_t>, double > p(inds, integral);
        this->mweights.insert(p);
        
        // next face
        go = this->msmplxIter.next();
    }
}


