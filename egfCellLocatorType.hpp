/**
 * Cell locator object
 */

#ifndef EGF_CELL_LOCATOR
#define EGF_CELL_LOCATOR

#include <vtkCellLocator.h>
#include <vector>

struct egfCellLocatorType {
    vtkCellLocator* cellLocator;
    std::vector<int> cellIds;
    bool treeHasBeenBuilt;
};

extern "C" {

/**
 * Constructor
 * @param self handle
 * @return 0 upon success
 */
int egfCellLocator_new(egfCellLocatorType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfCellLocator_del(egfCellLocatorType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfCellLocator_print(egfCellLocatorType** self);

/**
 * Load from file
 * @param self handle
 * @param grid unstructured grid instance
 * @return 0 upon success
 */
int egfCellLocator_setGrid(egfCellLocatorType** self,
	                       egfUnstructuredGrid3dType* grid);

/**
 * Set the number of cells per bucket used to accelerate the search for cells
 * @param self handle
 * @param numCellsPerBucket
 * @return 0 upon success
 */
 int egfCellLocator_setNumberOfCellsPerBucket(egfCellLocatorType** self, int numCellsPerBucket);

/**
 * Find cells in box
 * @param self handle
 * @param posMin lower box position
 * @param posMax upper box position
 * @return 0 upon success
 */
 int egfCellLocator_findCellsInBox(egfCellLocatorType** self, const double posMin[], const double posMax[]); 

/**
 * Get number of cells found
 * @param self handle
 * @param numCells (output)
 * @return 0 upon success
 */
 int egfCellLocator_getNumberOfCells(egfCellLocatorType** self, int* numCells);

 /**
  * Fill in cell IDs
  * @param self handle
  * @param cell ID array (output)
  * @return 0 upon success
  */
 int egfCellLocator_fillInCellIds(egfCellLocatorType** self, int* cellIds);  

}

#endif // EGF_CELL_LOCATOR