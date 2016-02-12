//
// A class for iterating over the elements (nodes, edges, faces, ...) of a simplex 
//
// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "Factorial.hpp"
#include <vector>
#include <string> // for size_t

#ifndef SIMPLEX_ITER_HPP
#define SIMPLEX_ITER_HPP

class SimplexIter {

public:

   /**
    * No argument constructor
    */
    SimplexIter() {}

   /**
    * Constructor
    * @param npdims number of parametric dimensions
    * @param order (0 for nodes, 1 for edges, 2 for faces, ...)
    */
    SimplexIter(size_t npdims, size_t order);

    virtual ~SimplexIter(){};

    /**
    * Setup
    * @param npdims number of parametric dimensions
    * @param order (0 for nodes, 1 for edges, 2 for faces)
    */
    void setup(size_t npdims, size_t order);

   /**
    * Set the iterator to the beginning
    */
    void begin();

   /**
    * Increment
    * @return true if new element is valid, false otherwise
    */
    bool next();

   /**
    * Accessor
    * @return array of node indices
    */
    const std::vector<size_t>& operator()() const {
        return mcurrentElem;
    }

    size_t getNumberOfElements() const {
        // binomial coefficient
        return factorial(mnpdims)/(factorial(morder)*factorial(mnpdims-morder));
    }


private:

    size_t mnpdims;
    size_t morder;
    size_t mindx;
    std::vector<size_t> mcurrentElem;

};

#endif // SIMPLEX_ITER_HPP