/**
 * Grid object
 */

#ifndef EGF_GRID
#define EGF_GRID


#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

struct egfGridType {
    vtkUnstructuredGridReader* reader;
    vtkUnstructuredGrid* ugrid;
};

extern "C" {

/**
 * Constructor
 * @param self handle
 * @return 0 upon success
 */
int egfGrid_new(egfGridType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfGrid_del(egfGridType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfGrid_print(egfGridType** self);

/**
 * Load from file
 * @param self handle
 * @param filename file name
 * @return 0 upon success
 */
int egfGrid_loadFromFile(egfGridType** self, const char* filename);

/**
 * Save to file
 * @param self handle
 * @param filename file name
 * @return 0 upon success
 */
int egfGrid_saveToFile(egfGridType** self, const char* filename);

/**
 * Get number of cells
 * @param self handle
 * @param numCells (output)
 * @return 0 upon success
 */
 int egfGrid_getNumberOfCells(egfGridType** self, int* numCells);

/**
 * Get number of cells
 * @param self handle
 * @param numPoints (output)
 * @return 0 upon success
 */
 int egfGrid_getNumberOfPoints(egfGridType** self, int* numPoints);

}

#endif // EGF_GRID