// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#include <iostream>
#include <vector>
#include <cmath>
#include <Simplex.hpp>
#include <SimplexIntersect.hpp>
#include <ESMCI_Exception.h>

const double eps = 1.e-10;

bool test1DBadTarget() {
    const size_t ndims = 1;
    SimplexIntersect si(ndims);
    Simplex smplx1(0, ndims);
    Simplex smplx2(0, ndims);
    // this should throw an exception
    try {
        si.solve(smplx1, smplx2);
        return false;
    }
    catch (ESMCI::Ex e) {
        std::cout << "testing exception: " << e.what() << '\n';
    }
    return true;
}

bool test1D() {
    const size_t ndims = 1;
    std::vector<double> pt(ndims);

    SimplexIntersect si(ndims);

    // inside test

    Simplex smplx1(0, ndims);
    pt[0] = 0.3;
    smplx1.setPoint(0, pt);

    Simplex smplx2(1, ndims);
    pt[0] = 0.0;
    smplx2.setPoint(0, pt);
    pt[0] = 1.0;
    smplx2.setPoint(1, pt);

    bool success = si.solve(smplx1, smplx2);
    if (!success) return false;

    bool ok = si.intersect(1.e-10);
    if (!ok) return false;

    // outside test

    pt[0] = -0.1;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (ok) return false;

    // on the boundary test

    pt[0] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;


    return true;
}

bool test2D() {
    const size_t ndims = 2;
    std::vector<double> pt(ndims);

    SimplexIntersect si(ndims);

    // inside test

    Simplex smplx1(1, ndims);
    pt[0] = 0.3; pt[1] = 0.4;
    smplx1.setPoint(0, pt);
    pt[0] = 0.0; pt[1] = 0.0;
    smplx1.setPoint(1, pt);

    Simplex smplx2(1, ndims);
    pt[0] = -1.0; pt[1] = 0.0;
    smplx2.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0;
    smplx2.setPoint(1, pt);

    bool success = si.solve(smplx1, smplx2);
    if (!success) return false;

    bool ok = si.intersect(1.e-10);
    if (!ok) return false;

    // outside test

    pt[0] = 0.1; pt[1] = 0.1;
    smplx1.setPoint(1, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (ok) return false;

    // on the boundary test

    pt[0] = -1.0; pt[1] = 0.0;
    smplx1.setPoint(1, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    return true;
}

bool test3D_point() {
    const size_t ndims = 3;
    std::vector<double> pt(ndims);

    SimplexIntersect si(ndims);

    // inside test

    Simplex smplx1(0, ndims);
    pt[0] = 0.3; pt[1] = 0.4; pt[1] = 0.5;
    smplx1.setPoint(0, pt);

    Simplex smplx2(3, ndims);
    pt[0] = 0.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx2.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx2.setPoint(1, pt);
    pt[0] = 0.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx2.setPoint(2, pt);
    pt[0] = 0.0; pt[1] = 0.0; pt[2] = 1.0;
    smplx2.setPoint(3, pt);

    bool success = si.solve(smplx1, smplx2);
    if (!success) return false;

    bool ok = si.intersect(1.e-10);
    if (!ok) return false;

    // outside test

    pt[0] = 1.1; pt[1] = 0.1; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (ok) return false;

    // on a node test

    pt[0] = 0.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    // on an edge test

    pt[0] = 0.5; pt[1] = 0.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    // on a face test

    pt[0] = 0.5; pt[1] = 0.5; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    return true;
}

bool test3D_line() {
    const size_t ndims = 3;
    std::vector<double> pt(ndims);

    SimplexIntersect si(ndims);

    // inside test

    Simplex smplx1(1, ndims);
    pt[0] = 0.3; pt[1] = 0.4; pt[1] = 0.5;
    smplx1.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 1.0;
    smplx1.setPoint(1, pt);

    Simplex smplx2(2, ndims);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx2.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx2.setPoint(1, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 1.0;
    smplx2.setPoint(2, pt);

    bool success = si.solve(smplx1, smplx2);
    if (!success) return false;

    bool ok = si.intersect(1.e-10);
    if (!ok) return false;

    // outside test

    pt[0] = 0.0; pt[1] = -1.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);
    pt[0] = 0.9; pt[1] = 0.0; pt[2] = 0.0;
    smplx1.setPoint(1, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (ok) return false;

    // on a node test

    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    // on an edge test

    pt[0] = 1.0; pt[1] = 0.5; pt[2] = 0.0;
    smplx1.setPoint(0, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    return true;
}

bool test3D_line_tangent() {
    const size_t ndims = 3;
    std::vector<double> pt(ndims);

    SimplexIntersect si(ndims);

    // barely inside test

    Simplex smplx1(1, ndims);
    pt[0] = 1.0; pt[1] = 0.0; pt[1] = 0.0;
    smplx1.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 1.0;
    smplx1.setPoint(1, pt);

    Simplex smplx2(2, ndims);
    pt[0] = 1.0; pt[1] = 0.0; pt[2] = 0.0;
    smplx2.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 0.0;
    smplx2.setPoint(1, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 1.0;
    smplx2.setPoint(2, pt);

    bool success = si.solve(smplx1, smplx2);
    if (!success) return false;

    bool ok = si.intersect(1.e-6);
    if (!ok) return false;

    // outside test

    pt[0] = 1.01; pt[1] = 0.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);
    pt[0] = 1.01; pt[1] = 1.0; pt[2] = 1.0;
    smplx1.setPoint(1, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (ok) return false;

    // on a node test

    pt[0] = 1.01; pt[1] = 0.0; pt[2] = 0.0;
    smplx1.setPoint(0, pt);
    pt[0] = 1.0; pt[1] = 1.0; pt[2] = 1.0;
    smplx1.setPoint(1, pt);    

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    // on an edge test

    pt[0] = 0.0; pt[1] = 0.5; pt[2] = 0.0;
    smplx1.setPoint(0, pt);
    pt[0] = 1.01; pt[1] = 0.5; pt[2] = 0.0;
    smplx1.setPoint(1, pt);

    success = si.solve(smplx1, smplx2);
    if (!success) return false;

    ok = si.intersect(1.e-10);
    if (!ok) return false;

    return true;
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

    if (!test3D_point()) {
        std::cout << "test3D_point FAILED\n";
        return 1;
    }

    if (!test3D_line()) {
        std::cout << "test3D_line FAILED\n";
        return 1;
    }

    if (!test3D_line_tangent()) {
        std::cout << "test3D_line_tangent FAILED\n";
        return 1;
    }

    std::cout << "PASSED\n";
    return 0;
}
