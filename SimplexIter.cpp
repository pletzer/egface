// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include "SimplexIter.hpp"

SimplexIter::SimplexIter(size_t npdims, size_t order) {
    this->setup(npdims, order);
}

void 
SimplexIter::setup(size_t npdims, size_t order) {
    mnpdims = npdims;
    morder = order;
    mindx = order;
    mcurrentElem.resize(order + 1);
    // starting element
    this->begin();
}

void
SimplexIter::begin() {
    for (size_t i = 0; i < mcurrentElem.size(); ++i) {
        mcurrentElem[i] = i;
    }
}

bool
SimplexIter::next() {
    size_t val = mcurrentElem[mindx] + 1;
    if (val <= mnpdims) {
        // increment the last index
        mcurrentElem[mindx] = val;
        return true;
    }
    else {
        // step back 
        int newIndx = (int) mindx - 1;
        while (newIndx >= 0) {
            val = mcurrentElem[newIndx] + 1;
            if (val <= mnpdims + newIndx - morder) {
                // reset 
                mcurrentElem[newIndx] = val;
                for (size_t j = newIndx + 1; j < morder + 1; ++j) {
                    mcurrentElem[j] = mcurrentElem[j- 1] + 1;
                }
                mindx = morder;
                return true;
            }
            else {
                newIndx--;
            }
        }
        // cannot increment, reached the end
        return false;
    }
}
