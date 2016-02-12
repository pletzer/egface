// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <iostream>
#include <vector>
#include <cmath>
#include <Simplex.hpp>
#include <FaceBasis.hpp>
#include <ESMCI_Exception.h>

const double eps = 1.e-10;

bool test2DBadTarget() {
    const size_t npdims = 2;
    FaceBasis bss(npdims);
    Simplex targetSmplx(0, npdims);
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

bool test2D() {
    const size_t order = 2;
    const size_t npdims = 2;
    std::map< std::vector<size_t>, double > weights;
    std::map< std::vector<size_t>, double > values;
    std::vector<double> xi(npdims);
    double val, exact;

    // create the basis function
    FaceBasis bss(npdims);

    // set the face field values
    std::vector<size_t> elem(order + 1);
    elem[0] = 0; elem[1] = 1; elem[2] = 2;
    values[elem] = 0.2;

    // create the target
    Simplex targetSmplx(order, npdims);

    // face 0 -> 1 -> 2
    xi[0] = 0.0; xi[1] = 0.0;
    targetSmplx.setPoint(0, xi);
    xi[0] = 1.0; xi[1] = 0.0;
    targetSmplx.setPoint(1, xi);
    xi[0] = 0.0; xi[1] = 1.0;
    targetSmplx.setPoint(2, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = 0.2;
    if (std::abs(val - exact) > eps) {
        return false;
    }

    return true;
}


bool test3D() {
    const size_t order = 2;
    const size_t npdims = 3;
    std::map< std::vector<size_t>, double > weights;
    std::map< std::vector<size_t>, double > values;
    std::vector<double> xi(npdims);
    double val, exact;

    // create the basis function
    FaceBasis bss(npdims);

    // set the face field values
    std::vector<size_t> elem(order + 1);
    elem[0] = 0; elem[1] = 1; elem[2] = 2;
    values[elem] = 0.2;

    elem[0] = 0; elem[1] = 1; elem[2] = 3;
    values[elem] = 0.3;

    elem[0] = 0; elem[1] = 2; elem[2] = 3;
    values[elem] = 0.4;

    elem[0] = 1; elem[1] = 2; elem[2] = 3;
    values[elem] = 0.5;

    // create the target
    Simplex targetSmplx(order, npdims);

    // face 0 -> 1 -> 2
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(0, xi);
    xi[0] = 1.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(1, xi);
    xi[0] = 0.0; xi[1] = 1.0; xi[2] = 0.0;
    targetSmplx.setPoint(2, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = 0.2;
    if (std::abs(val - exact) > eps) {
        return false;
    }

    // face 0 -> 1 -> 3
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(0, xi);
    xi[0] = 1.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(1, xi);
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 1.0;
    targetSmplx.setPoint(2, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = 0.3;
    if (std::abs(val - exact) > eps) {
        return false;
    }

    // face 0 -> 2 -> 3
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(0, xi);
    xi[0] = 0.0; xi[1] = 1.0; xi[2] = 0.0;
    targetSmplx.setPoint(1, xi);
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 1.0;
    targetSmplx.setPoint(2, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = 0.4;
    if (std::abs(val - exact) > eps) {
        return false;
    }

    // face 1 -> 2 -> 3
    xi[0] = 1.0; xi[1] = 0.0; xi[2] = 0.0;
    targetSmplx.setPoint(0, xi);
    xi[0] = 0.0; xi[1] = 1.0; xi[2] = 0.0;
    targetSmplx.setPoint(1, xi);
    xi[0] = 0.0; xi[1] = 0.0; xi[2] = 1.0;
    targetSmplx.setPoint(2, xi);

    bss.computeWeights(targetSmplx);
    weights = bss.getWeights();

    // interpolate
    val = 0;
    for (std::map< std::vector<size_t>, double >::const_iterator it = values.begin();
        it != values.end(); ++it) {
        val += it->second * weights[it->first];
    }

    // check
    exact = 0.5;
    if (std::abs(val - exact) > eps) {
        return false;
    }

    return true;
}

int main() {

    if (!test2DBadTarget()) {
        std::cout << "test2DBadTarget FAILED\n";
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
