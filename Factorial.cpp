// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <string> // for size_t

size_t factorial(size_t n) {
    size_t res = 1;
    if (n <= 1) return res;
    for (size_t i = 1; i < n + 1; ++i) {
        res *= i;
    }
    return res;
}
