// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "NodalBasis.hpp"
#include "ESMCI_Exception.h"

void NodalBasis::computeWeights(const Simplex& smplx) {

    this->mweights.clear();

    if (smplx.getNumberOfParametricDimensions() != 0) {
        // error
        throw ESMCI::Ex() << "NodalBasis::computeWeights can only compute weights for a 0-simplex (point)";
    }
    const std::vector<double>& xiCoords = smplx.getPoint(0); // target is a point
    double sumXi = 0;
    for (size_t i = 0; i < xiCoords.size(); ++i) {
        sumXi += xiCoords[i];
    }

    // iterate over the nodes of the source cell
    this->msmplxIter.begin();
    bool go = true;
    while (go) {
        const std::vector<size_t>& inds = this->msmplxIter();
        size_t i = inds[0];

        // evaluate the basis function at xiCoords
        double val = 0;
        if (i > 0) {
            val = xiCoords[i-1];
        }
        else {
            val = 1.0 - sumXi;
        }

        std::pair< std::vector<size_t>, double > p(inds, val);
        this->mweights.insert(p);

        go = this->msmplxIter.next();
    }
}


