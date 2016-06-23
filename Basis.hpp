//
// Base class for linearly interpolating fields on simplices
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <map>
#include <vector>
#include <string> // for size_t

#include "Simplex.hpp"
#include "SimplexIter.hpp"

#ifndef BASIS_HPP
#define BASIS_HPP

class Basis {

public:

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    */
    Basis(size_t npdims) {
        mnpdims = npdims;
    }

    virtual ~Basis(){};

   /**
    * Get the list of elements (nodes, edges, faces, ...) and corresponding 
    * interpolation weights
    * @return map element -> weight
    */
    const std::map< std::vector<size_t>, double >& getWeights() const {
         return mweights;
    }

   /**
    * Compute the interpolation weights
    * @param smplx a target simplex whose vertices are expressed in parametric coordinates
    * @return value
    */
    virtual void computeWeights(const Simplex& smplx) = 0;

   /**
    * Get the number of elements 
    * @return number
    */
    size_t getNumberOfElements() {
        return mnumElements;
    }

protected:

    size_t mnpdims;
    size_t mnumElements;
    std::map< std::vector<size_t>, double > mweights;
    SimplexIter msmplxIter;

};

#endif // NODAL_BASIS_HPP
