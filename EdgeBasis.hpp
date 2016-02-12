//
// A class for linearly interpolating edge fields on simplices
//
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "Basis.hpp" // base class

#ifndef EDGE_BASIS_HPP
#define EDGE_BASIS_HPP

class EdgeBasis : public Basis {

public:

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    */
    EdgeBasis(size_t npdims) : Basis(npdims) {
        this->msmplxIter.setup(npdims, 1);
        this->mnumElements = this->msmplxIter.getNumberOfElements();
    }

    virtual ~EdgeBasis(){}

   /**
    * Compute the interpolation weights
    * @param smplx a simplex whose vertices are expressed in parametric coordinates
    * @note assumes smplx is fully contained 
    */
    void computeWeights(const Simplex& smplx);

private:

};

#endif // EDGE_BASIS_HPP
