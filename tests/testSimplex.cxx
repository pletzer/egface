// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <iostream>
#include <vector>
#include <cmath>
#include <Simplex.hpp>

const double eps = 1.e-10;

bool testPointIn0D() {
    const size_t npdims = 0;
    const size_t ndims = 0;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 0.0) < eps) {
        return true;
    }
    return false;
}


bool testPointIn1D() {
    const size_t npdims = 0;
    const size_t ndims = 1;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 0.0) < eps) {
        return true;
    }
    return false;
}

bool testLineIn1D() {
    const size_t npdims = 1;
    const size_t ndims = 1;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0;
    smplx.setPoint(1, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 1.0) < eps) {
        return true;
    }
    return false;
}

bool testPointIn2D() {
    const size_t npdims = 0;
    const size_t ndims = 2;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 0.0) < eps) {
        return true;
    }
    return false;
}

bool testLineIn2D() {
    const size_t npdims = 1;
    const size_t ndims = 2;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0;
    smplx.setPoint(1, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 1.0) < eps) {
        return true;
    }
    return false;
}

bool testTriangleIn2D() {
    const size_t npdims = 2;
    const size_t ndims = 2;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0;
    smplx.setPoint(1, pt);
    pt[0] = 0.0; pt[1] = 1.0;
    smplx.setPoint(2, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 0.5) < eps) {
        return true;
    }
    return false;
}

bool testPointIn3D() {
    const size_t npdims = 0;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 0.0) < eps) {
        return true;
    }
    return false;
}

bool testLineIn3D() {
    const size_t npdims = 1;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx.setPoint(1, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 1.0) < eps) {
        return true;
    }
    return false;
}

bool testTriangleIn3D() {
    const size_t npdims = 2;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx.setPoint(1, pt);
    pt[0] = 0.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx.setPoint(2, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 1.0/2.0) < eps) {
        return true;
    }
    return false;
}

bool testTetrahedronIn3D() {
    const size_t npdims = 3;
    const size_t ndims = 3;
    Simplex smplx(npdims, ndims);
    std::vector<double> pt(ndims, 0.0);
    smplx.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx.setPoint(1, pt);
    pt[0] = 0.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx.setPoint(2, pt);
    pt[0] = 0.0; pt[1] = 0.0; pt[2] = 1.0;
    smplx.setPoint(3, pt);
    double vol = smplx.getVolume();
    if (std::abs(vol - 1.0/6.0) < eps) {
        return true;
    }
    return false;
}

int main() {

    if (!testPointIn0D()) {
        std::cout << "testPointIn0D FAILED\n";
        return 1;
    }

    if (!testPointIn1D()) {
        std::cout << "testPointIn1D FAILED\n";
        return 1;
    }
    if (!testLineIn1D()) {
        std::cout << "testLineIn1D FAILED\n";
        return 1;
    }

    if (!testPointIn2D()) {
        std::cout << "testPointIn2D FAILED\n";
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

    if (!testPointIn3D()) {
        std::cout << "testPointIn3D FAILED\n";
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