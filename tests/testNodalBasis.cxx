// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <iostream>
#include <vector>
#include <cmath>
#include <Simplex.hpp>
#include <NodalBasis.hpp>
#include <ESMCI_Exception.h>

const double eps = 1.e-10;

bool test1DBadTarget() {
    const size_t npdims = 1;
    NodalBasis bss(npdims);
    Simplex targetSmplx(1, npdims);
    // this should throw an exception
    try {
        bss.computeWeights(targetSmplx);
        return false;
    }
    catch (ESMCI::Ex e) {
        std::cout << "testing exception: " << e.what() << '\n';
    }
    return true;
}

bool test1D() {
    const size_t order = 0;
    const size_t npdims = 1;
    std::map< std::vector<size_t>, double > weights;
    std::map< std::vector<size_t>, double > values;
    std::vector<double> xi(npdims);
    double val, exact;    

    // create the basis function
    NodalBasis bss(npdims);

    // set the nodal field values
    std::vector<size_t> elem(order + 1);
    elem[0] = 0; 
    values[elem] = 0.2;
    elem[0] = 1;
    values[elem] = 0.5;

    // create the target
    Simplex targetSmplx(order, npdims);
    xi[0] = 0.4;
    targetSmplx.setPoint(0, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = (1.0 - xi[0])*0.2 + xi[0]*0.5;
    if (std::abs(val - exact) < eps) {
        return true;
    }
    return false;
}

bool test2D() {
    const size_t order = 0;
    const size_t npdims = 2;
    std::map< std::vector<size_t>, double > weights;
    std::map< std::vector<size_t>, double > values;
    std::vector<double> xi(npdims);
    double val, exact;    

    // create the basis function
    NodalBasis bss(npdims);

    // set the nodal field values
    std::vector<size_t> elem(order + 1);
    elem[0] = 0; 
    values[elem] = 0.2;
    elem[0] = 1;
    values[elem] = 0.3;
    elem[0] = 2;
    values[elem] = 0.4;

    // create the target
    Simplex targetSmplx(order, npdims);
    xi[0] = 0.5; xi[1] = 0.4;
    targetSmplx.setPoint(0, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = (1.0 - xi[0] - xi[1])*0.2 + xi[0]*0.3 + xi[1]*0.4;
    if (std::abs(val - exact) < eps) {
        return true;
    }
    return false;
}

bool test3D() {
    const size_t order = 0;
    const size_t npdims = 3;
    std::map< std::vector<size_t>, double > weights;
    std::map< std::vector<size_t>, double > values;
    std::vector<double> xi(npdims);
    double val, exact;    

    // create the basis function
    NodalBasis bss(npdims);

    // set the nodal field values
    std::vector<size_t> elem(order + 1);
    elem[0] = 0; 
    values[elem] = 0.2;
    elem[0] = 1;
    values[elem] = 0.3;
    elem[0] = 2;
    values[elem] = 0.4;
    elem[0] = 3;
    values[elem] = 0.5;

    // create the target
    Simplex targetSmplx(order, npdims);
    xi[0] = 0.5; xi[1] = 0.4; xi[2] = 0.3;
    targetSmplx.setPoint(0, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = (1.0 - xi[0] - xi[1] - xi[2])*0.2 + xi[0]*0.3 + xi[1]*0.4 + xi[2]*0.5;
    if (std::abs(val - exact) < eps) {
        return true;
    }
    return false;
}

int main() {

    if (!test1DBadTarget()) {
        std::cout << "test1DBadTarget FAILED\n";
        return 1;        
    }

    if (!test1D()) {
        std::cout << "test1D FAILED\n";
        return 1;
    }
    if (!test2D()) {
        std::cout << "test2D FAILED\n";
        return 1;
    }
    if (!test3D()) {
        std::cout << "test3D FAILED\n";
        return 1;
    }

    std::cout << "PASSED\n";
    return 0;
}
