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
 * @return 0 upon success
 */
int egfGrid_new(egfGridType** self);

/**
 * Destructor
 * @return 0 upon success
 */
int egfGrid_del(egfGridType** self);

/**
 * Load from file
 * @param filename file name
 * @return 0 upon success
 */
int egfGrid_loadFromFile(egfGridType** self, const char* filename);

/**
 * Save to file
 * @param filename file name
 * @return 0 upon success
 */
int egfGrid_saveToFile(egfGridType** self, const char* filename);
}

#endif // EGF_GRID