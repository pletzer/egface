//
// A class for linearly interpolating face fields on simplices
//
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "Basis.hpp" // base class

#ifndef FACE_BASIS_HPP
#define FACE_BASIS_HPP

class FaceBasis : public Basis {

public:

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    */
    FaceBasis(size_t npdims) : Basis(npdims) {
        this->msmplxIter.setup(npdims, 2);
        this->mnumElements = this->msmplxIter.getNumberOfElements();
    }

    virtual ~FaceBasis(){}

   /**
    * Compute interpolation weights
    * @param smplx a target simplex whose vertices are expressed in parametric coordinates
    * @note assumes smplx is fully contained 
    */
    void computeWeights(const Simplex& smplx);

private:

};

#endif // FACE_BASIS_HPP