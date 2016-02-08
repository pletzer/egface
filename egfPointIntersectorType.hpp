/**
 * Point intersector object
 */

#ifndef EGF_POINT_INTERSECTOR
#define EGF_POINT_INTERSECTOR

#include "egfGridType.hpp"
#include <vtkCellLocator.h>
#include <vector>
#include <map>
#include <set>

struct egfPointIntersectorType {

    // source grid
	vtkUnstructuredGrid* ugrid;

    // cell finder
    vtkCellLocator* cellLocator;

    // collection of intersection points for each source grid cell
    std::map<vtkIdType, std::set<std::vector<double> > > intersectPoints;
    double tol;

    // Add a single point
    void addEntry(vtkIdType cellId, const std::vector<double>& point) {

        // Is there an entry for cellId?
        std::map<vtkIdType, std::set<std::vector<double> > >::iterator
            it = this->intersectPoints.find(cellId);

        if (it == this->intersectPoints.end()) {
            std::set<std::vector<double> > pointsInCell;
            pointsInCell.insert(point);
            std::pair<vtkIdType, std::set<std::vector<double> > > cp(cellId, pointsInCell);
            this->intersectPoints.insert(cp);
        }
        else {
            it->second.insert(point);
        }
    }

    // Add a collection of points to the map
    void addEntry(vtkIdType cellId, const std::set<std::vector<double> >& pointsInCell) {
            
        if (pointsInCell.size() == 0) {
            // Nothing to do
            return;
        }

        // Is there an entry for cellId?
        std::map<vtkIdType, std::set<std::vector<double> > >::iterator
            it = this->intersectPoints.find(cellId);

        if (it == this->intersectPoints.end()) {
            // No, create a new entry and insert
            std::pair<vtkIdType, std::set<std::vector<double> > > cp(cellId, pointsInCell);
            this->intersectPoints.insert(cp);
        }
            else {
            // Yes, insert
            it->second.insert(pointsInCell.begin(), pointsInCell.end());
        }

    }
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
 * Set number of cells per bucket
 * @param self handle
 * @param numCells number of cells
 * @return 0 upon success
 */
int egfPointIntersector_setNumberOfCellsPerBucket(egfPointIntersectorType** self, int numCells);

/**
 * Find source grid cell containing point
 * @param self handle
 * @param p0 vertex
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithPoint(egfPointIntersectorType** self, 
                                      const double p0[]); 

/**
 * Find grid points intersecting with line segment
 * @param self handle
 * @param p0 first line vertex
 * @param p1 second line vertex
 * @param addEndPoints whether (1) or not (0) the end points should be added
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithLine(egfPointIntersectorType** self, 
 	                                  const double p0[],
 	                                  const double p1[],
                                      int addEndPoints); 

/**
 * Find grid points intersecting with triangle
 * @param self handle
 * @param p0 first triangle vertex
 * @param p1 second triangle vertex
 * @param p2 third triangle vertex
 * @param addEndPoints whether (1) or not (0) the end points should be added
 * @return 0 upon success
 */
 int egfPointIntersector_gridWithTriangle(egfPointIntersectorType** self, 
 	                                      const double p0[],
 	                                      const double p1[],
 	                                      const double p2[],
                                          int addEndPoints); 

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
 * Get number of intersected cells 
 * @param self handle
 * @param numCells number of cells (output)
 * @return 0 upon success
 */
 int egfPointIntersector_getNumberOfCells(egfPointIntersectorType** self, int* numCells);

/**
 * Fill in the intersected cell Ids 
 * @param self handle
 * @param cellIds cell ids (output)
 * @return 0 upon success
 * @note caller allocates memory of cellids, size numCells
 */
 int egfPointIntersector_fillInCellIds(egfPointIntersectorType** self, int* cellIds);

/**
 * Get number of intersection points in cell 
 * @param self handle
 * @param cellId cell Id
 * @param numPoints number of points (output)
 * @return 0 upon success
 */
 int egfPointIntersector_getNumberOfPointsInCell(egfPointIntersectorType** self, int cellId, int* numPoints);

/**
 * Fill in the intersection points belonging to a cell 
 * @param self handle
 * @param cellId cell Id
 * @param points as a flat array (output)
 * @return 0 upon success
 * @note caller allocates memory of points, size 3*numPoints
 */
 int egfPointIntersector_fillInPointsInCell(egfPointIntersectorType** self, int cellId, double* points);

}

#endif // EGF_POINT