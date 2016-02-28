/** 
 * Test point intersector
 */

#include <iostream>
#include <egfPointIntersectorType.hpp>

int main(int argc, char**argv) {

    const int ADD_END_POINTS = 1;
    
    if (argc < 2) {
        std::cout << "ERROR: Must provide path to vtk file\n";
        return 1;
    }
    char* filename = argv[1];

    // Load a grid from file
    egfUnstructuredGrid3dType* grid = NULL;
    int ier;
    ier = egfUnstructuredGrid3d_new(&grid);
    ier = egfUnstructuredGrid3d_loadFromFile(&grid, filename);

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
    ier = egfPointIntersector_gridWithLine(&intrsctr, p0, p1, ADD_END_POINTS);

    // Print the object
    ier = egfPointIntersector_print(&intrsctr);

    // Destructor
    ier = egfPointIntersector_del(&intrsctr);

    // Delete grid
    ier = egfUnstructuredGrid3d_del(&grid);
    
    return 0;
}
