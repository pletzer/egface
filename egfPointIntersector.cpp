#include "egfPointIntersectorType.hpp"
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <iostream>
#include <algorithm>

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
	std::cout << "found " << (*self)->intersectPoints.size() << " points\n";
	for (std::set< std::vector<double> >::const_iterator p = (*self)->intersectPoints.begin(); 
         p != (*self)->intersectPoints.end(); ++p) {
        std::cout << "point ";
        for (size_t j = 0; j < (*p).size(); ++j) {
		  std::cout << (*p)[j] << ", ";
        }
        std::cout << '\n';
	}
    return 0;
}

int egfPointIntersector_setGrid(egfPointIntersectorType** self, egfGridType* grid) {
    (*self)->ugrid = grid->ugrid;
    (*self)->cellLocator->SetDataSet(grid->ugrid);
    (*self)->cellLocator->BuildLocator();
    return 0;
}

int egfPointIntersector_getNumberOfPoints(egfPointIntersectorType** self, int* numPoints) {
    *numPoints = (int) (*self)->intersectPoints.size();
    return 0;
}

int egfPointIntersector_fillInPoints(egfPointIntersectorType** self, double* points) {
    size_t i = 0;
    for (std::set< std::vector<double> >::const_iterator p = (*self)->intersectPoints.begin(); 
         p != (*self)->intersectPoints.end(); ++p) {
        for (size_t j = 0; j < (*p).size(); ++j) {
            points[i*3 + j] = (*p)[j];
        }
        i++;
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


int egfPointIntersector_gridWithLine(egfPointIntersectorType** self, 
                                     const double p0[], 
                                     const double p1[]) {

    // VTK does nto enforce const'ness
    std::vector<double> pa(p0, p0 + 3);
    std::vector<double> pb(p1, p1 + 3);

    // Construct triangle by building a one cell unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(2);
    points->SetPoint(0, p0);
    points->SetPoint(1, p1);
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
    (*self)->cellLocator->FindCellsAlongLine(&pa[0], &pb[0], (*self)->tol, cellIds);
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
        vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
        int numFaces = cell->GetNumberOfFaces();
        for (int k = 0; k < numFaces; ++k) {
            vtkCell* face = cell->GetFace(k);
            int res = face->IntersectWithLine(&pa[0], &pb[0], (*self)->tol, t, &pt.front(), pcoords, subId);
            if (res) {
                (*self)->intersectPoints.insert(pt);
            }
        }
    }

    // Add the segment's vertices
    (*self)->intersectPoints.insert(std::vector<double>(p0, p0 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p1, p1 + 3));

    return 0;
}

int egfPointIntersector_gridWithTriangle(egfPointIntersectorType** self, 
                                         const double p0[], 
                                         const double p1[], 
                                         const double p2[]) {

    // Construct triangle by building a one cell unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(3);
    points->SetPoint(0, p0);
    points->SetPoint(1, p1);
    points->SetPoint(2, p2);
    ug->SetPoints(points);
    ug->Allocate(1, 1);
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds->SetNumberOfIds(3);
    for (vtkIdType i = 0; i < 3; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_TRIANGLE, ptIds);

    // Compute the bounding box
    double* bbox = ug->GetBounds();

    // Find all the grid cells in the bounding box
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection
        return 0;
    }

    // Compute the intersection between each grid cell face with the triangle's edges
    double t; // parametric position along the line
    std::vector<double> pt(3); // intersection point
    double pcoords[] = {0, 0, 0}; // tetrahedron parametric coordinates
    int subId; // not used
    double* pa; // start point on the line
    double* pb; // end point of the line
    vtkCell* tri = ug->GetCell(0);
    int numEdges = tri->GetNumberOfEdges();
    for (int i = 0; i < numEdges; ++i) {
        vtkCell* edge = tri->GetEdge(i);
        pa = points->GetPoint(edge->GetPointId(0));
        pb = points->GetPoint(edge->GetPointId(1));
        // Iterate over the grid cells in the box
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
            int numFaces = cell->GetNumberOfFaces();
            for (int k = 0; k < numFaces; ++k) {
                vtkCell* face = cell->GetFace(k);
                int res = face->IntersectWithLine(pa, pb, (*self)->tol, t, &pt.front(), pcoords, subId);
                if (res) {
                    (*self)->intersectPoints.insert(pt);
                }
            }
        }
    }

    // Compute the intersection between each grid cell edge with the triangle
    // Iterate over the grid cells in the box
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
        int numEdges = cell->GetNumberOfEdges();
        for (int k = 0; k < numEdges; ++k) {
            vtkCell* edge = cell->GetEdge(k);
            pa = points->GetPoint(edge->GetPointId(0));
            pb = points->GetPoint(edge->GetPointId(1));
            int res = tri->IntersectWithLine(pa, pb, (*self)->tol, t, &pt.front(), pcoords, subId);
            if (res) {
                (*self)->intersectPoints.insert(pt);
            }                
        }
    }

    // Add the triangle's vertices
    (*self)->intersectPoints.insert(std::vector<double>(p0, p0 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p1, p1 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p2, p2 + 3));

    // Add all the grid cell vertices that are inside this triangle
    double* closestPoint;
    double dist2;
    double weights[] = {0, 0};
    // Iterate over the grid cells
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
        vtkIdType numPoints = cell->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double* p = points->GetPoint(cell->GetPointId(i));
            int res = tri->EvaluatePosition(p, closestPoint, subId, pcoords, dist2, weights);
            bool inside = true;
            double sum = 0;
            for (size_t k = 0; k < 2; ++k) {
                inside &= pcoords[k] > -(*self)->tol;
                inside &= pcoords[k] < 1 + (*self)->tol;
                sum += pcoords[k];
            }
            if (inside && sum < 1 + (*self)->tol) {
                (*self)->intersectPoints.insert(std::vector<double>(p, p + 3));
            }
        }
    }

    return 0;
}

int egfPointIntersector_gridWithTetrahedron(egfPointIntersectorType** self, 
                                            const double p0[], 
                                            const double p1[], 
                                            const double p2[], 
                                            const double p3[]) {

    // Construct tetrahedron by building a one cell unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(4);
    points->SetPoint(0, p0);
    points->SetPoint(1, p1);
    points->SetPoint(2, p2);
    points->SetPoint(3, p3);
    ug->SetPoints(points);
    ug->Allocate(1, 1);
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_TETRA, ptIds);

    // Compute the bounding box
    double* bbox = ug->GetBounds();

    // Find all the grid cells in the bounding box
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection
        return 0;
    }

    // Compute the intersection between each grid cell face with the tet's edges
    double t; // parametric position along the line
    std::vector<double> pt(3); // intersection point
    double pcoords[] = {0, 0}; // tetrahedron parametric coordinates
    int subId; // not used
    double* pa; // start point on the line
    double* pb; // end point of the line
    vtkCell* tet = ug->GetCell(0);
    int numEdges = tet->GetNumberOfEdges();
    for (int i = 0; i < numEdges; ++i) {
        vtkCell* edge = tet->GetEdge(i);
        pa = points->GetPoint(edge->GetPointId(0));
        pb = points->GetPoint(edge->GetPointId(1));
        // Iterate over the grid cells in the box
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
            int numFaces = cell->GetNumberOfFaces();
            for (int k = 0; k < numFaces; ++k) {
                vtkCell* face = cell->GetFace(k);
                int res = face->IntersectWithLine(pa, pb, (*self)->tol, t, &pt.front(), pcoords, subId);
                if (res) {
                    (*self)->intersectPoints.insert(pt);
                }
            }
        }
    }

    // Compute the intersection between each grid cell edge and the tet's faces
    int numFaces = tet->GetNumberOfFaces();
    for (int i = 0; i < numFaces; ++i) {
        vtkCell* face = tet->GetFace(i);
        // Iterate over the grid cells in the box
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
            vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
            int numEdges = cell->GetNumberOfEdges();
            for (int k = 0; k < numEdges; ++k) {
                vtkCell* edge = cell->GetEdge(k);
                pa = points->GetPoint(edge->GetPointId(0));
                pb = points->GetPoint(edge->GetPointId(1));
                int res = face->IntersectWithLine(pa, pb, (*self)->tol, t, &pt.front(), pcoords, subId);
                if (res) {
                    (*self)->intersectPoints.insert(pt);
                }                
            }
        }
    }

    // Add the tet's vertices
    (*self)->intersectPoints.insert(std::vector<double>(p0, p0 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p1, p1 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p2, p2 + 3));
    (*self)->intersectPoints.insert(std::vector<double>(p3, p3 + 3));

    // Add all the grid cell vertices that are inside this tet
    double* closestPoint;
    double dist2;
    double weights[] = {0, 0, 0};
    // Iterate over the grid cells
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
        vtkIdType numPoints = cell->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double* p = points->GetPoint(cell->GetPointId(i));
            int res = tet->EvaluatePosition(p, closestPoint, subId, pcoords, dist2, weights);
            bool inside = true;
            double sum = 0;
            for (size_t k = 0; k < 3; ++k) {
                inside &= pcoords[k] > -(*self)->tol;
                inside &= pcoords[k] < 1 + (*self)->tol;
                sum += pcoords[k];
            }
            if (inside && sum < 1 + (*self)->tol) {
                (*self)->intersectPoints.insert(std::vector<double>(p, p + 3));
            }
        }
    }

    return 0;
}

 
} // extern "C"
