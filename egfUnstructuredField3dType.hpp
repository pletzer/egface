/**
 * Field object
 */

#ifndef EGF_UNSTRUCTURED_FIELD_3D
#define EGF_UNSTRUCTURED_FIELD_3D


#include "egfUnstructuredGrid3dType.hpp"
#include "SimplexIter.hpp"
#include <vtkDoubleArray.h>
#include <vector>
#include <map>
#include <string>

struct egfUnstructuredField3dType {
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
int egfUnstructuredField3d_new(egfUnstructuredField3dType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfUnstructuredField3d_del(egfUnstructuredField3dType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfUnstructuredField3d_print(egfUnstructuredField3dType** self);
    
/** 
 * Set field order
 * @param order (0 = nodal, 1 = edge, 2 = face, 3 = cell)
 * @return 0 upon success
 */
 int egfUnstructuredField3d_setOrder(egfUnstructuredField3dType** self, int order);
    
 /**
  * Set grid
  * @param grid grid
  * @return 0 upon success
  */
 int egfUnstructuredField3d_setGrid(egfUnstructuredField3dType** self,
                                    egfUnstructuredGrid3dType* grid);
    
 /**
  * Get number of elements (nodes, edges, faces) in a cell
  * @param numElems number of elements (output)
  * @return 0 upon success
  */
 int egfUnstructuredField3d_getNumberOfElements(egfUnstructuredField3dType** self,
                                                int* numElems);
    
 /**
  * Get element
  * @param elem cell element index in the range 0...numElems - 1
  * @param inds array of size order + 1 for local cell point indices (output)
  * @return 0 upon success
  */
 int egfUnstructuredField3d_getElement(egfUnstructuredField3dType** self,
                                       int elem, int inds[]);
    
  /**
   * Set the field values for a given cell element
   * @param elem cell element index in the range 0...numElems - 1
   * @param vals arrays of size number of cells
   * @return 0 upon success
   */
 int egfUnstructuredField3d_setValues(egfUnstructuredField3dType** self,
                                      int elem, const double vals[]);

}

#endif // EGF_UNSTRUCTURED_FIELD_3D
