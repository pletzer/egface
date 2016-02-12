#include <iostream>
#include <vector>
#include <SimplexIter.hpp>

void printArray(const std::vector<size_t>& arr) {
    for (size_t i = 0; i < arr.size(); ++i) {
        std::cout << arr[i] << ", ";
    }
    std::cout << '\n';
}

void test(size_t npdims, size_t order) {
    SimplexIter sit(npdims, order);
    bool go = true;
    while (go) {
        printArray(sit());
        go = sit.next();
    }
}

int main() {

    std::cout << "0d/order = 0 (points)\n";
    test(0, 0);

    std::cout << "1d/order = 0 (points)\n";
    test(1, 0);
    std::cout << "1d/order = 1 (edges)\n";
    test(1, 1);

    std::cout << "2d/order = 0 (points)\n";
    test(2, 0);
    std::cout << "2d/order = 1 (edges)\n";
    test(2, 1);
    std::cout << "2d/order = 2 (faces)\n";
    test(2, 2);

    std::cout << "3d/order = 0 (points)\n";
    test(3, 0);
    std::cout << "3d/order = 1 (edges)\n";
    test(3, 1);
    std::cout << "3d/order = 2 (faces) \n";
    test(3, 2);
    std::cout << "3d/order = 3 (cells) \n";
    test(3, 3);

    std::cout << "4d/order = 0 (points)\n";
    test(4, 0);
    std::cout << "4d/order = 1 (edges)\n";
    test(4, 1);
    std::cout << "4d/order = 2 (faces)\n";
    test(4, 2);
    std::cout << "4d/order = 3 (3-cells)\n";
    test(4, 3);
    std::cout << "4d/order = 4 (4-cells)\n";
    test(4, 4);
    
    std::cout << "PASSED\n";
    return 0;
}