/**
 * Grid object
 */

#ifndef EGF_GRID
#define EGF_GRID


#include <vtkUnstructuredGrid.h>

struct egfGridType {
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

/**
 * Get the domain bounds
 * @param self handle
 * @param posMin lower box point (output)
 * @param posMin upper box point (output)
 * @return 0 upon success
 */
 int egfGrid_getDomainBounds(egfGridType** self, double* posMin, double* posMax);

/**
 * Get the number of cell elements of type "order"
 * @param self handle
 * @param order 0=nodal, 1=edge, 2=face, 3=cell
 * @param numElems number of points, edges, faces, etc (output)
 * @return 0 upon success
 */
 int egfGrid_getNumberOfElements(egfGridType** self, int order, int* numElems);

/**
 * Get the point connectivity of points, edges, faces, cells for each cell
 * @param self handle
 * @param order 0=nodal, 1=edge, 2=face, 3=cell
 * @param ptConnect array of size numCells * numElems * (order + 1) (output)
 * @return 0 upon success
 * @note caller owns ptConnect
 */
 int egfGrid_getElementConnectivity(egfGridType** self, int order, int ptConnect[]);

}

#endif // EGF_GRID