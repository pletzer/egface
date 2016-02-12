//
// A class for representing simplices (line segments, triangles, tets, etc)
//
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <vector>
#include <string> // for size_t

#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

class Simplex {

public:

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    * @param ndims number of space dimensions
    * @note this will construct a default simplex with vertices (0, ...), 
    *       (1, 0, ...), (0, 1, ...)
    */
    Simplex(size_t npdims, size_t ndims);

    virtual ~Simplex(){};

   /**
    * Set vertex
    * @param index index of the vertex
    * @param pt point
    */
    void setPoint(size_t index, const std::vector<double>& pt);

   /**
    * Get vertices
    * @return vertex
    */
    const std::vector< std::vector<double> >& getPoints() const {
        return mvertices;
    }



   /**
    * Get the simplex parametric coordinates
    * @param position point in physical space
    * @param xiCoords parametric coordinates (output)
    * @param tol tolerance in xiCoords for determining whether point is inside simplex
    * @return return true if point is inside simplex, false otherwise
    */
    bool getParametricCoordinates(const std::vector<double>& position, 
                                  double tol, std::vector<double>& xiCoords);

   /**
    * Get the volume (length, area, ...) of the simplex
    * @return positive definite number 
    */
    double getVolume() const;

   /**
    * Get vertex 
    * @param index index in the range 0...npdims (incl.)
    * @return vertex
    */
    const std::vector<double>& getPoint(size_t index) const {
        return mvertices[index];
    }

   /**
    * Get number of parametric dimensions
    * @return number
    */
    size_t getNumberOfParametricDimensions() const {
        return mnpdims;
    }

   /**
    * Get number of space dimensions
    * @return number
    */
    size_t getNumberOfDimensions() const {
        return mndims;
    }    

private:

    double getDeterminant(const std::vector<double>& mat) const;

    int solve(std::vector<double>& rhs) const; // rhs is solution on output

    size_t mnpdims;
    size_t mndims;
    std::vector< std::vector<double> > mvertices;
};

#endif // SIMPLEX_HPP