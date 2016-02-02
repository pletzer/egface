/**
 * Point intersector object
 */

#ifndef EGF_POINT_INTERSECTOR
#define EGF_POINT_INTERSECTOR

#include "egfGridType.hpp"
#include <vtkCellLocator.h>
#include <vector>

struct egfPointIntersectorType {
	vtkUnstructuredGrid* ugrid;
    vtkCellLocator* cellLocator;
    std::vector< std::vector<double> > intersectPoints;
    double tol;
};

extern "C" {

/**
 * Constructor
 * @param self handle
 * @return 0 upon success
 */
int egfPointIntersector_new(egfPointIntersectorType** self);

/**
 * Destructor
 * @param self handle
 * @return 0 upon success
 */
int egfPointIntersector_del(egfPointIntersectorType** self);

/**
 * Print 
 * @param self handle
 * @return 0 upon success
 */
 int egfPointIntersector_print(egfPointIntersectorType** self);

/**
 * Set grid
 * @param self handle
 * @param grid unstructured grid instance
 * @return 0 upon success
 */
int egfPointIntersector_setGrid(egfPointIntersectorType** self, egfGridType* grid);

/**
 * Set tolerance
 * @param self handle
 * @param tol tolerance
 * @return 0 upon success
 */
int egfPointIntersector_setTolerance(egfPointIntersectorType** self, double tol);

/**
 * Find grid points intersecting with line segment
 * @param self handle
 * @param p0 first tetrahedron vertex
 * @param p1 second tetrahedron vertex
 * @param p2 third tetrahedron vertex
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithLine(egfPointIntersectorType** self, 
 	                                  const double p0[],
 	                                  const double p1[]); 

/**
 * Find grid points intersecting with triangle
 * @param self handle
 * @param p0 first tetrahedron vertex
 * @param p1 second tetrahedron vertex
 * @param p2 third tetrahedron vertex
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithTriangle(egfPointIntersectorType** self, 
 	                                      const double p0[],
 	                                      const double p1[],
 	                                      const double p2[]); 

/**
 * Find grid points intersecting with tetrahedron
 * @param self handle
 * @param p0 first tetrahedron vertex
 * @param p1 second tetrahedron vertex
 * @param p2 third tetrahedron vertex
 * @param p3 fourth tetrahedron vertex
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithTetrahedron(egfPointIntersectorType** self, 
 	                                     	 const double p0[],
 	                                     	 const double p1[],
 	                                     	 const double p2[],
 	                                     	 const double p3[]); 

/**
 * Get number of points found
 * @param self handle
 * @param numPoints (output)
 * @return 0 upon success
 */
 int egfPointIntersector_getNumberOfPoints(egfPointIntersectorType** self, int* numPoints);

 /**
  * Fill in points
  * @param self handle
  * @param point array (output)
  * @return 0 upon success
  */
 int egfPointIntersector_fillInPoints(egfPointIntersectorType** self, double* points);  

}

#endif // EGF_POINT