//
// A class for linearly interpolating nodal fields on simplices
//
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "Factorial.hpp"
#include "Basis.hpp" // base class

#ifndef NODAL_BASIS_HPP
#define NODAL_BASIS_HPP

class NodalBasis : public Basis {

public:

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    */
    NodalBasis(size_t npdims) : Basis(npdims) {
      this->msmplxIter.setup(npdims, 0);
      this->mnumElements = this->msmplxIter.getNumberOfElements();
    }

    virtual ~NodalBasis(){};

   /**
    * Compute interpolation weights
    * @param smplx a target simplex whose vertices are expressed in parametric coordinates
    * @note assumes smplx is fully contained 
    */
    void computeWeights(const Simplex& smplx);

};

#endif // NODAL_BASIS_HPP