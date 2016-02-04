/** 
 * Test point intersector
 */

#include <iostream>
#include <egfPointIntersectorType.hpp>

int main(int argc, char**argv) {
    
    if (argc < 2) {
        std::cout << "ERROR: Must provide path to vtk file\n";
        return 1;
    }
    char* filename = argv[1];

    // Load a grid from file
    egfGridType* grid = NULL;
    int ier;
    ier = egfGrid_new(&grid);
    ier = egfGrid_loadFromFile(&grid, filename);

    egfPointIntersectorType* intrsctr = NULL;

    // Constructor
    ier = egfPointIntersector_new(&intrsctr);

    // Set the tolerance
    double tol = 1.e-12;
    ier = egfPointIntersector_setTolerance(&intrsctr, tol);

    // Set the number of cells per bucket
    ier = egfPointIntersector_setNumberOfCellsPerBucket(&intrsctr, 100);

    // Set the grid, should be called after setNumberOfCellsPerBucket
    ier = egfPointIntersector_setGrid(&intrsctr, grid);

    // Compute all the intersection points
    const double p0[] = {1, 0.5, 0};
    const double p1[] = {1, 0.5, 6};
    std::cerr << "*** entering gridWithLine\n";
    ier = egfPointIntersector_gridWithLine(&intrsctr, p0, p1);
    std::cerr << "*** exiting gridWithLine\n";

    // Print the object
    ier = egfPointIntersector_print(&intrsctr);

    // Destructor
    ier = egfPointIntersector_del(&intrsctr);

    // Delete grid
    ier = egfGrid_del(&grid);
    
    return 0;
}
