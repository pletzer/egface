// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "EdgeBasis.hpp"
#include "ESMCI_Exception.h"


void EdgeBasis::computeWeights(const Simplex& smplx) {

    this->mweights.clear();

    if (smplx.getNumberOfParametricDimensions() != 1) {
        // error
        throw ESMCI::Ex() << "EdgeBasis::computeWeights can only compute weights for a 1-simplex (line)";
    }
    const std::vector<double>& xiA = smplx.getPoint(0); // target is a line
    const std::vector<double>& xiB = smplx.getPoint(1);
    std::vector<double> deltaXi(this->mnpdims);
    double sumXiAs = 0;
    double sumDeltaXis = 0;
    for (size_t i = 0; i < this->mnpdims; ++i) {
        deltaXi[i] = xiB[i] - xiA[i];
        sumXiAs += xiA[i];
        sumDeltaXis += deltaXi[i];
    }

    // iterate over the edges of the source cell
    this->msmplxIter.begin();
    bool go = true;
    while (go) {
        const std::vector<size_t>& inds = this->msmplxIter();
        // i, j are in ascending order
        // i >= 0, j != i, so j > 0
        size_t i = inds[0];
        size_t j = inds[1];
        size_t jm1 = j - 1; // j >= 1

        // compute the integral
        double integral;
        if (i > 0) {
            size_t im1 = i - 1;
            integral = deltaXi[jm1]*(xiA[im1] + 0.5*deltaXi[im1]) 
                     - deltaXi[im1]*(xiA[jm1] + 0.5*deltaXi[jm1]);
        }
        else {
            integral = deltaXi[jm1]*(1.0 - sumXiAs - 0.5*sumDeltaXis)
                     + sumDeltaXis*(xiA[jm1] + 0.5*deltaXi[jm1]);
        }

        std::pair< std::vector<size_t>, double > p(inds, integral);
        this->mweights.insert(p);

        go = this->msmplxIter.next();
    }
}


