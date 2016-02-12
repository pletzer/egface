//
// A class to compute the intersection of two simplices
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <limits>
#include <vector>
#include <string> // for size_t
#include "Simplex.hpp"

#ifndef SIMPLEX_INTERSECT_HPP
#define SIMPLEX_INTERSECT_HPP

class SimplexIntersect {

public:

   /**
    * Constructor
    * @param ndims number of dimensions
    */
    SimplexIntersect(size_t ndims) {
        
        mndims = ndims;
        mmat.resize(ndims * ndims);
        msol.resize(ndims);

        // very small number use to displace the first simplex's coordinates
        // to avoid a singular system
        moffset = 12.345678901234 * std::numeric_limits<double>::epsilon();
    }

    virtual ~SimplexIntersect(){};

   /**
    * Solve the system of equations
    * @param smplx1 first simplex
    * @param smplx2 second simplex
    * @return true if a solution was found
    */
    bool solve(const Simplex& smplx1, const Simplex& smplx2);
    
   /**
    * Get the solution
    * @return array of parametric coordinates for the first simplex (host), followed
    *         by the parametric coordinates of the second simplex (target)
    */
    std::vector<double> getSolution() const {
         return msol;
    }

   /**
    * Does the target simplex intersect with the host simplex?
    * @param tol tolerance
    * @return true if the simplices intersect, false otherwise
    */
    bool intersect(double tol);

private:

    size_t mndims;
    size_t mndims1;
    std::vector<double> msol;
    std::vector<double> mmat;
    double moffset;
};

#endif // SIMPLEX_INTERSECT_HPP
