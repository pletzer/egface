#include "egfPointIntersectorType.hpp"
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <iostream>
#include <algorithm>

//#define DEBUG

extern "C" {

int egfPointIntersector_new(egfPointIntersectorType** self) {
	(*self) = new egfPointIntersectorType();
    (*self)->ugrid = NULL;
	(*self)->cellLocator = vtkCellLocator::New();
    (*self)->cellLocator->SetNumberOfCellsPerBucket(20);
	(*self)->tol = 1.e-10;
	return 0;
}

int egfPointIntersector_del(egfPointIntersectorType** self) {
    (*self)->ugrid = NULL;
	(*self)->cellLocator->Delete();
	delete *self;
    return 0;
}

int egfPointIntersector_print(egfPointIntersectorType** self) {
	std::cout << "egfPointIntersector:\n";
	std::cout << "found " << (*self)->intersectPoints.size() << " cells\n";
    for (std::map<vtkIdType, std::set<std::vector<double> > >::const_iterator
      it = (*self)->intersectPoints.begin(); it != (*self)->intersectPoints.end(); ++it) { 
        std::cout << "\tcell " << it->first << " has " << it->second.size() << " points\n";
        size_t i = 0;
        for (std::set<std::vector<double> >::const_iterator
            it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            std::cout << "\t\tpoint " << i << ": ";
            for (size_t j = 0; j < (*it2).size(); ++j) {
                std::cout << (*it2)[j] << ", ";
            }
            std::cout << '\n';
            i++;
        }
	}
    return 0;
}

int egfPointIntersector_setGrid(egfPointIntersectorType** self,
                                egfUnstructuredGrid3dType* grid) {
    (*self)->ugrid = grid->ugrid;
    (*self)->cellLocator->SetDataSet(grid->ugrid);
    (*self)->cellLocator->BuildLocator();
    return 0;
}

int egfPointIntersector_getNumberOfCells(egfPointIntersectorType** self, int* numCells) {
    *numCells = (int) (*self)->intersectPoints.size();
    return 0;
}

int egfPointIntersector_fillInCellIds(egfPointIntersectorType** self, int* cellIds) {
    size_t i = 0;
    for (std::map<vtkIdType, std::set<std::vector<double> > >::const_iterator 
      it = (*self)->intersectPoints.begin(); it != (*self)->intersectPoints.end(); ++it) {
        cellIds[i] = (int) it->first;
        i++;
    }
    return 0;
}

int egfPointIntersector_getNumberOfPointsInCell(egfPointIntersectorType** self, int cellId, int* numPoints) {
    *numPoints = 0;
    std::map<vtkIdType, std::set<std::vector<double> > >::const_iterator
      it = (*self)->intersectPoints.find((vtkIdType) cellId);
    if (it != (*self)->intersectPoints.end()) {
        *numPoints = (int) it->second.size();
    }
    return 0;
}

int egfPointIntersector_fillInPointsInCell(egfPointIntersectorType** self, int cellId, double* points) {
    std::map<vtkIdType, std::set<std::vector<double> > >::const_iterator
      it = (*self)->intersectPoints.find((vtkIdType) cellId);
    if (it != (*self)->intersectPoints.end()) {
        size_t i = 0;
        for (std::set<std::vector<double> >::const_iterator it2 = it->second.begin(); 
          it2 != it->second.end(); ++it2) {
            for (size_t j = 0; j < 3; ++j) {
                points[i*3 + j] = (*it2)[j];
            }
            i++;
        }
    }    
    return 0;
}

int egfPointIntersector_setTolerance(egfPointIntersectorType** self, double tol) {
    (*self)->tol = tol;
    return 0;
}

int egfPointIntersector_setNumberOfCellsPerBucket(egfPointIntersectorType** self, int numCells) {
    (*self)->cellLocator->SetNumberOfCellsPerBucket(numCells);
    return 0;
}

int egfPointIntersector_gridWithPoint(egfPointIntersectorType** self, 
                                      const double p0[]) {

    double dist2;
    std::vector<double> point(p0, p0 + 3);
    std::vector<double> closestPoint(3);
    double tol2 = (*self)->tol * (*self)->tol;
    int subId;
    double pcoords[] = {0., 0., 0.};
    double weights[] = {0., 0., 0., 0., 0., 0., 0., 0.};
    vtkIdType cellId = (*self)->ugrid->FindCell(&point[0], NULL, 0, 
        tol2, subId, pcoords, weights);
    if (cellId >= 0) {
        (*self)->addEntry(cellId, point);

        // Also search among neighbors since the point my be shared
        vtkCell* cell = (*self)->ugrid->GetCell(cellId);
        vtkIdList* ptIds = cell->GetPointIds();
        vtkIdType numPoints = ptIds->GetNumberOfIds();
        vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> sharedPtIds = vtkSmartPointer<vtkIdList>::New();
        sharedPtIds->SetNumberOfIds(1);
        for (vtkIdType i = 0; i < numPoints; ++i) {
            sharedPtIds->SetId(0, ptIds->GetId(i));
            (*self)->ugrid->GetCellNeighbors(cellId, sharedPtIds, neighborCellIds);
            vtkIdType numNeighborCells = neighborCellIds->GetNumberOfIds();
            for (vtkIdType j = 0; j < numNeighborCells; ++j) {
                vtkIdType neighborCellId = neighborCellIds->GetId(j);
                vtkCell* neighborCell = (*self)->ugrid->GetCell(neighborCellId);
                vtkIdType cellId2 = (*self)->ugrid->FindCell(&point[0], 
                    neighborCell, neighborCellId, tol2, subId, pcoords, weights);
                if (cellId2 >= 0) {
                    (*self)->addEntry(cellId2, point);
                }
            }
        }
    }

    return 0;

}

int egfPointIntersector_gridWithLine(egfPointIntersectorType** self, 
                                     const double p0[], 
                                     const double p1[], 
                                     int addEndPoints) {

    std::vector<double> pa(p0, p0 + 3);
    std::vector<double> pb(p1, p1 + 3);
    double bounds[6];

    // Truncate the segment if it goes beyond the grid bounds
    (*self)->ugrid->GetBounds(bounds);
    const double xmins[] = {bounds[0], bounds[2], bounds[4]};
    const double xmaxs[] = {bounds[1], bounds[3], bounds[5]};
    for (size_t j = 0; j < 3; ++j) {
        pa[j] = (pa[j] < xmins[j]? xmins[j]: pa[j]);
        pa[j] = (pa[j] > xmaxs[j]? xmaxs[j]: pa[j]);
        pb[j] = (pb[j] < xmins[j]? xmins[j]: pb[j]);
        pb[j] = (pb[j] > xmaxs[j]? xmaxs[j]: pb[j]);
    }

    // Build a one cell unstructured grid representing the line segment
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(2);
    points->SetPoint(0, &pa[0]);
    points->SetPoint(1, &pb[0]);

    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints(points);
    ug->Allocate(1, 1);

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds->SetNumberOfIds(2);
    for (vtkIdType i = 0; i < 2; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_LINE, ptIds);

    // Find the cells along the line
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    (*self)->cellLocator->FindCellsAlongLine(&pa[0], &pb[0],
                                             (*self)->tol, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection
        return 0;
    }

    // Compute the intersection between each grid cell face with the segment
    double t; // parametric position along the line
    std::vector<double> pt(3); // intersection point
    double pcoords[] = {0, 0, 0}; // tetrahedron parametric coordinates
    int subId; // not used
    // Iterate over the grid cells along the segment
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {

        vtkIdType cellId = cellIds->GetId(j);

        // Collect all the intersection points in this cell
        std::set<std::vector<double> > pointsInCell;

        vtkCell* cell = (*self)->ugrid->GetCell(cellId);

        // Iterate over faces of the cell
        int numFaces = cell->GetNumberOfFaces();
        for (int k = 0; k < numFaces; ++k) {
            vtkCell* face = cell->GetFace(k);
            int res = face->IntersectWithLine(&pa[0], &pb[0],
                                              (*self)->tol, t, &pt[0], 
                                              pcoords, subId);
            if (res) {
                pointsInCell.insert(pt);
            }
        }

        (*self)->addEntry(cellId, pointsInCell);

    }

    // Add the segment's vertices
    if (addEndPoints) {
        // Add the triangle's vertices
        egfPointIntersector_gridWithPoint(self, p0);
        egfPointIntersector_gridWithPoint(self, p1);
    }

    return 0;
}

int egfPointIntersector_gridWithTriangle(egfPointIntersectorType** self, 
                                         const double p0[], 
                                         const double p1[], 
                                         const double p2[], 
                                         int addEndPoints) {

    // Copy the points so we can pass double* arrays, which VTK expects
    std::vector<double> pa(p0, p0 + 3);
    std::vector<double> pb(p1, p1 + 3);
    std::vector<double> pc(p2, p2 + 3);

    // Construct triangle by building a one cell unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(3);
    points->SetPoint(0, &pa[0]);
    points->SetPoint(1, &pb[0]);
    points->SetPoint(2, &pc[0]);
    ug->SetPoints(points);
    ug->Allocate(1, 1);
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds->SetNumberOfIds(3);
    for (vtkIdType i = 0; i < 3; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_TRIANGLE, ptIds);

    // Compute the bounding box of the triangle
    double* bbox = ug->GetBounds();

    // Find all the grid cells within the triangle's bounds
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection possible
        return 0;
    }

    double t; // parametric position along the line
    std::vector<double> pt(3); // intersection point
    double pcoords[] = {0, 0, 0}; // triangle parametric coordinates
    int subId; // not used
    std::vector<double> pBeg(3); // start point on the line
    std::vector<double> pEnd(3); // end point of the line
    vtkCell* tri = ug->GetCell(0);

    // Compute the intersection between each cell edge with the triangle
    // Iterate over the grid cells in the box
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkIdType cellId = cellIds->GetId(j);
        std::set<std::vector<double> > pointsInCell;
        vtkCell* cell = (*self)->ugrid->GetCell(cellId);
        int numEdges = cell->GetNumberOfEdges();
        for (int k = 0; k < numEdges; ++k) {
            vtkCell* edge = cell->GetEdge(k);
            vtkIdType iBeg = edge->GetPointId(0);
            vtkIdType iEnd = edge->GetPointId(1);
            (*self)->ugrid->GetPoints()->GetPoint(iBeg, &pBeg[0]);
            (*self)->ugrid->GetPoints()->GetPoint(iEnd, &pEnd[0]);
            int res = tri->IntersectWithLine(&pBeg[0], &pEnd[0], 
                (*self)->tol, t, &pt[0], pcoords, subId);
            if (res) {
                pointsInCell.insert(pt);
            }                
        }

        (*self)->addEntry(cellId, pointsInCell);

    }

    // Compute the intersection between each grid cell face with 
    // the triangle's edges
    const int NOT_ADD_VERTEX_POINTS = 0;
    egfPointIntersector_gridWithLine(self, p0, p1, NOT_ADD_VERTEX_POINTS);
    egfPointIntersector_gridWithLine(self, p1, p2, NOT_ADD_VERTEX_POINTS);
    egfPointIntersector_gridWithLine(self, p2, p0, NOT_ADD_VERTEX_POINTS);

    if (addEndPoints) {
        // Add the triangle's vertices
        egfPointIntersector_gridWithPoint(self, p0);
        egfPointIntersector_gridWithPoint(self, p1);
        egfPointIntersector_gridWithPoint(self, p2);
    }

    return 0;
}

int egfPointIntersector_gridWithTetrahedron(egfPointIntersectorType** self, 
                                            const double p0[], 
                                            const double p1[], 
                                            const double p2[], 
                                            const double p3[]) {

    std::vector<double> pa(p0, p0 + 3);
    std::vector<double> pb(p1, p1 + 3);
    std::vector<double> pc(p2, p2 + 3);
    std::vector<double> pd(p3, p3 + 3);

    // Construct tetrahedron by building a one cell unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(4);
    points->SetPoint(0, &pa[0]);
    points->SetPoint(1, &pb[0]);
    points->SetPoint(2, &pc[0]);
    points->SetPoint(3, &pd[0]);
    ug->SetPoints(points);
    ug->Allocate(1, 1);
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_TETRA, ptIds);

    // Compute the bounding box of the tetrahedron
    double* bbox = ug->GetBounds();

    // Find all the grid cells in the bounding box
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection
        return 0;
    }

    double pcoords[] = {0, 0, 0}; // tetrahedron parametric coordinates
    int subId; // not used
    vtkCell* tet = ug->GetCell(0);
    double weights[8]; // size is max number of nodes per cell

    // Add all the grid cell vertices that are inside this tet
    std::vector<double> point(3);
    std::vector<double> closestPoint(3);
    double dist2;
    // Iterate over the grid cells
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkIdType cellId = cellIds->GetId(j);
        std::set<std::vector<double> > pointsInCell;
        vtkCell* cell = (*self)->ugrid->GetCell(cellId);
        vtkIdType numPoints = cell->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; ++i) {
            (*self)->ugrid->GetPoints()->GetPoint(cell->GetPointId(i), &point[0]);
            int res = tet->EvaluatePosition(&point[0], &closestPoint[0], 
                subId, pcoords, dist2, weights);
            bool inside = true;
            double sum = 0;
            for (size_t k = 0; k < 3; ++k) {
                inside &= pcoords[k] > -(*self)->tol;
                inside &= pcoords[k] < 1 + (*self)->tol;
                sum += pcoords[k];
            }
            if (inside && sum < 1 + (*self)->tol) {
                pointsInCell.insert(point);
            }
        }

        (*self)->addEntry(cellId, pointsInCell);

    }

    // Compute the intersection between each grid cell edge and the tet's faces
    egfPointIntersector_gridWithTriangle(self, p0, p1, p2, 0);
    egfPointIntersector_gridWithTriangle(self, p1, p2, p3, 0);
    egfPointIntersector_gridWithTriangle(self, p2, p0, p3, 0);
    egfPointIntersector_gridWithTriangle(self, p0, p1, p3, 0);

    // Compute the intersection between each grid cell face with the tet's edges
    egfPointIntersector_gridWithLine(self, p0, p1, 0);
    egfPointIntersector_gridWithLine(self, p1, p2, 0);
    egfPointIntersector_gridWithLine(self, p2, p0, 0);
    egfPointIntersector_gridWithLine(self, p0, p3, 0);
    egfPointIntersector_gridWithLine(self, p1, p3, 0);
    egfPointIntersector_gridWithLine(self, p2, p3, 0);


    // Add the tet's vertices
    egfPointIntersector_gridWithPoint(self, p0);
    egfPointIntersector_gridWithPoint(self, p1);
    egfPointIntersector_gridWithPoint(self, p2);
    egfPointIntersector_gridWithPoint(self, p3);

    return 0;
}

 
} // extern "C"
