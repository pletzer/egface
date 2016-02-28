/**
 * Field object
 */

#ifndef EGF_FIELD
#define EGF_FIELD


#include "egfUnstructuredGrid3dType.hpp"
#include "SimplexIter.hpp"
#include <vtkDoubleArray.h>
#include <vector>
#include <map>
#include <string>

struct egfFieldType {
    vtkUnstructuredGrid* ugrid;
    std::vector<vtkDoubleArray*> dataArrays;
    std::map<int, std::vector<size_t> > elems;
    SimplexIter smplxIt;
    std::string type;
    int order;
};

extern "C" {

/**
 * Constructor
 * @param self handle
 * @return 0 upon success
 */
int egfField_new(egfFieldType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfField_del(egfFieldType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfField_print(egfFieldType** self);
    
/** 
 * Set field order
 * @param order (0 = nodal, 1 = edge, 2 = face, 3 = cell)
 * @return 0 upon success
 */
 int egfField_setOrder(egfFieldType** self, int order);
    
 /**
  * Set grid
  * @param grid grid
  * @return 0 upon success
  */
 int egfField_setGrid(egfFieldType** self, egfUnstructuredGrid3dType* grid);
    
 /**
  * Get number of elements (nodes, edges, faces) in a cell
  * @param numElems number of elements (output)
  * @return 0 upon success
  */
 int egfField_getNumberOfElements(egfFieldType** self, int* numElems);
    
 /**
  * Get element
  * @param elem cell element index in the range 0...numElems - 1
  * @param inds array of size order + 1 for local cell point indices (output)
  * @return 0 upon success
  */
 int egfField_getElement(egfFieldType** self, int elem, int inds[]);
    
  /**
   * Set the field values for a given cell element
   * @param elem cell element index in the range 0...numElems - 1
   * @param vals arrays of size number of cells
   * @return 0 upon success
   */
 int egfField_setValues(egfFieldType** self, int elem, const double vals[]);

}

#endif // EGF_FIELD
