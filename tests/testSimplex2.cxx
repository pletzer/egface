#include <iostream>
#include <vector>
#include <cmath>
#include <Simplex.hpp>

const double eps = 1.e-10;

bool testLineIn1D() {
    const size_t npdims = 1;
    const size_t ndims = 1;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && std::abs(xi[0] - 0.2) < eps) {
        return true;
    }
    return false;
}

bool testLineIn2D() {
    const size_t npdims = 1;
    const size_t ndims = 2;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2; pt[1] = 0.1;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && std::abs(xi[0] - 0.2) < eps) {
        return true;
    }
    return false;
}

bool testTriangleIn2D() {
    const size_t npdims = 2;
    const size_t ndims = 2;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2; pt[1] = 0.1;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && 
        std::abs(xi[0] - 0.2) < eps &&
        std::abs(xi[1] - 0.1) < eps) {
        return true;
    }
    return false;
}

bool testLineIn3D() {
    const size_t npdims = 1;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2; pt[1] = 0.0; pt[2] = 0.0;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && 
        std::abs(xi[0] - 0.2) < eps) {
        return true;
    }
    return false;
}

bool testTriangleIn3D() {
    const size_t npdims = 2;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2; pt[1] = 0.1; pt[2] = 0.0;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && 
        std::abs(xi[0] - 0.2) < eps &&
        std::abs(xi[1] - 0.1) < eps) {
        return true;
    }
    return false;
}

bool testTetrahedronIn3D() {
    const size_t npdims = 3;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    pt[0] = 0.2; pt[1] = 0.1; pt[2] = 0.3;
    std::vector<double> xi(npdims);
    bool inside = smplx.getParametricCoordinates(pt, eps, xi);
    if (inside && 
        std::abs(xi[0] - 0.2) < eps &&
        std::abs(xi[1] - 0.1) < eps &&
        std::abs(xi[2] - 0.3) < eps) {
        return true;
    }
    return false;
}

int main() {

    if (!testLineIn1D()) {
        std::cout << "testLineIn1D FAILED\n";
        return 1;
    }

    if (!testLineIn2D()) {
        std::cout << "testLineIn2D FAILED\n";
        return 1;
    }
    if (!testTriangleIn2D()) {
        std::cout << "testTriangleIn2D FAILED\n";
        return 1;
    }

    if (!testLineIn3D()) {
        std::cout << "testLineIn3D FAILED\n";
        return 1;
    }
    if (!testTriangleIn3D()) {
        std::cout << "testTriangleIn3D FAILED\n";
        return 1;
    }
    if (!testTetrahedronIn3D()) {
        std::cout << "testTetrahedronIn3D FAILED\n";
        return 1;
    }
    
    std::cout << "PASSED\n";
    return 0;
}
