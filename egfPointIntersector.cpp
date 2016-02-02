#include "egfPointIntersectorType.hpp"
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <iostream>
#include <algorithm>

extern "C" {

int egfPointIntersector_new(egfPointIntersectorType** self) {
	(*self) = new egfPointIntersectorType();
    (*self)->ugrid = NULL;
	(*self)->cellLocator = vtkCellLocator::New();
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
	std::cout << "found points: ";
	for (size_t i = 0; i < (*self)->intersectPoints.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
		  std::cout << (*self)->intersectPoints[i][j] << ", ";
        }
        std::cout << '\n';
	}
    return 0;
}

int egfPointIntersector_setGrid(egfPointIntersectorType** self, egfGridType* grid) {
    (*self)->ugrid = grid->ugrid;
    (*self)->cellLocator->SetDataSet(grid->ugrid);
    return 0;
}

 int egfPointIntersector_getNumberOfPoints(egfPointIntersectorType** self, int* numPoints) {
    *numPoints = (int) (*self)->intersectPoints.size();
    return 0;
 }

 int egfPointIntersector_fillInPoints(egfPointIntersectorType** self, double* points) {
    for (size_t i = 0; i < (*self)->intersectPoints.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            points[i*3 + j] = (*self)->intersectPoints[i][j];
        }
    }
    return 0;
 } 

int egfPointIntersector_gridWithLine(egfPointIntersectorType** self, 
                                     const double p0[], 
                                     const double p1[]) {
    return 0;
}

int egfPointIntersector_gridWithTriangle(egfPointIntersectorType** self, 
                                         const double p0[], 
                                         const double p1[], 
                                         const double p2[]) {
    return 0;
}

int egfPointIntersector_gridWithTetrahedron(egfPointIntersectorType** self, 
                                            const double p0[], 
                                            const double p1[], 
                                            const double p2[], 
                                            const double p3[]) {

    // Construct tetrahedron by building a one cell unstructured grid
    vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(4);
    points->SetPoint(0, p0);
    points->SetPoint(1, p1);
    points->SetPoint(2, p2);
    points->SetPoint(3, p3);
    ug->SetPoints(points);
    ug->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    ug->InsertNextCell(VTK_TETRA, ptIds);
    ptIds->Delete();

    // Compute the bounding box
    double* bbox = ug->GetBounds();

    // Find all the grid cells in the bounding box
    vtkIdList* cellIds = vtkIdList::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIds);
    if (cellIds->GetNumberOfIds() == 0) {
        // No intersection
        return 0;
    }

    // Compute the intersection between each grid cell face with the tet's edges
    double t; // parametric position along the line
    std::vector<double> pt(3); // intersection point
    double pcoords[] = {0, 0, 0}; // tetrahedron parametric coordinates
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
                    (*self)->intersectPoints.push_back(pt);
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
                    (*self)->intersectPoints.push_back(pt);
                }                
            }
        }
    }

    // Add the tet's vertices
    (*self)->intersectPoints.push_back(std::vector<double>(p0, p0 + 4));
    (*self)->intersectPoints.push_back(std::vector<double>(p1, p1 + 4));
    (*self)->intersectPoints.push_back(std::vector<double>(p2, p2 + 4));
    (*self)->intersectPoints.push_back(std::vector<double>(p3, p3 + 4));

    // Add all the grid cell vertices that are inside this tet
    double* closestPoint;
    double dist2;
    double weights[] = {0, 0, 0};
    // Itereate over the grid cells
    for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
        vtkCell* cell = (*self)->ugrid->GetCell(cellIds->GetId(j));
        vtkIdType numPoints = cell->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double* p = points->GetPoint(cell->GetPointId(i));
            int res = cell->EvaluatePosition(p, closestPoint, subId, pcoords, dist2, weights);
            if (res) {
                (*self)->intersectPoints.push_back(std::vector<double>(p, p + 4));
            }
        }
    }

    // Clean up
    cellIds->Delete();
    points->Delete();
    ug->Delete();

    return 0;
}

 
} // extern "C"
