/**
 * Grid object
 */

#ifndef EGF_UNSTRUCTURED_GRID_3D
#define EGF_UNSTRUCTURED_GRID_3D


#include <vtkUnstructuredGrid.h>

struct egfUnstructuredGrid3dType {
    vtkUnstructuredGrid* ugrid;
};

extern "C" {

/**
 * Constructor
 * @param self handle
 * @return 0 upon success
 */
int egfUnstructuredGrid3d_new(egfUnstructuredGrid3dType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfUnstructuredGrid3d_del(egfUnstructuredGrid3dType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfUnstructuredGrid3d_print(egfUnstructuredGrid3dType** self);

/**
 * Load from file
 * @param self handle
 * @param filename file name
 * @return 0 upon success
 */
int egfUnstructuredGrid3d_loadFromFile(egfUnstructuredGrid3dType** self, const char* filename);

/**
 * Save to file
 * @param self handle
 * @param filename file name
 * @return 0 upon success
 */
int egfUnstructuredGrid3d_saveToFile(egfUnstructuredGrid3dType** self, const char* filename);

/**
 * Get number of cells
 * @param self handle
 * @param numCells (output)
 * @return 0 upon success
 */
 int egfUnstructuredGrid3d_getNumberOfCells(egfUnstructuredGrid3dType** self, int* numCells);

/**
 * Get number of cells
 * @param self handle
 * @param numPoints (output)
 * @return 0 upon success
 */
 int egfUnstructuredGrid3d_getNumberOfPoints(egfUnstructuredGrid3dType** self, int* numPoints);

/**
 * Get the domain bounds
 * @param self handle
 * @param posMin lower box point (output)
 * @param posMin upper box point (output)
 * @return 0 upon success
 */
 int egfUnstructuredGrid3d_getDomainBounds(egfUnstructuredGrid3dType** self,
 	                                       double* posMin, double* posMax);

}

#endif // EGF_UNSTRUCTURED_GRID_3D