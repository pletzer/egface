#include "egfUnstructuredGrid3dType.hpp"
#include "SimplexIter.hpp"
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <iostream>

extern "C" {

int egfUnstructuredGrid3d_new(egfUnstructuredGrid3dType** self) {
	(*self) = new egfUnstructuredGrid3dType();
	(*self)->ugrid = vtkUnstructuredGrid::New();
	return 0;
}

int egfUnstructuredGrid3d_del(egfUnstructuredGrid3dType** self) {
	(*self)->ugrid->Delete();
	delete *self;
    return 0;
}

int egfUnstructuredGrid3d_print(egfUnstructuredGrid3dType** self) {
	std::cout << "egfUnstructuredGrid3d:\n";
	vtkIdType numPoints = (*self)->ugrid->GetNumberOfPoints();
	double point[] = {0, 0, 0};
	for (vtkIdType i = 0; i < numPoints; ++i) {
		(*self)->ugrid->GetPoint(i, point);
		std::cout << i << " point " << point[0] << ", " << point[1] << ", " << point[2] << '\n';
	}
	vtkIdType numCells = (*self)->ugrid->GetNumberOfCells();
	vtkIdType* ptIds;
	vtkIdType npts;
	for (vtkIdType i = 0; i < numCells; ++i) {
		int cellType = (*self)->ugrid->GetCellType(i);
		(*self)->ugrid->GetCellPoints(i, npts, ptIds);
		std::cout << i << " cell type " << cellType << " point ids ";
		for (vtkIdType j = 0; j < npts; ++j) {
			std::cout << ptIds[j] << ' ';
		}
		std::cout << '\n';
	}

    return 0;
}

int egfUnstructuredGrid3d_loadFromFile(egfUnstructuredGrid3dType** self, const char* filename) {
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    (*self)->ugrid->DeepCopy(reader->GetOutput());
	return 0;
}

int egfUnstructuredGrid3d_saveToFile(egfUnstructuredGrid3dType** self, const char* filename) {
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(filename);
#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput((*self)->ugrid);
#else
    writer->SetInputData((*self)->ugrid);
#endif
    writer->Update();
	return 0;
}

 int egfUnstructuredGrid3d_getNumberOfCells(egfUnstructuredGrid3dType** self, int* numCells) {
 	*numCells = (int) (*self)->ugrid->GetNumberOfCells();
 	return 0;
 }

 int egfUnstructuredGrid3d_getNumberOfPoints(egfUnstructuredGrid3dType** self, int* numPoints) {
 	*numPoints = (int) (*self)->ugrid->GetNumberOfPoints();
 	return 0; 	
 }

 int egfUnstructuredGrid3d_getDomainBounds(egfUnstructuredGrid3dType** self, double* posMin, double* posMax) {
    (*self)->ugrid->GetPoints()->ComputeBounds();
    double* bbox = (*self)->ugrid->GetPoints()->GetBounds();
    posMin[0] = bbox[0];
    posMax[0] = bbox[1];
    posMin[1] = bbox[2];
    posMax[1] = bbox[3];
    posMin[2] = bbox[4];
    posMax[2] = bbox[5];
    return 0;
 }

 int egfUnstructuredGrid3d_getNumberOfElements(egfUnstructuredGrid3dType** self, int order, int* numElems) {
 	SimplexIter smplxIt(3, order);
 	*numElems = smplxIt.getNumberOfElements();
 	return 0;
 }

 int egfUnstructuredGrid3d_getElementConnectivity(egfUnstructuredGrid3dType** self, int order, int ptConnect[]) {
 	int res = 0;
 	SimplexIter smplxIt(3, order);
 	size_t numElems = smplxIt.getNumberOfElements();
	vtkIdType* ptIds;
	vtkIdType npts;
 	vtkIdType numCells = (*self)->ugrid->GetNumberOfCells();
 	for (vtkIdType i = 0; i < numCells; ++i) {
 		(*self)->ugrid->GetCellPoints(i, npts, ptIds);
 		if (npts == 4) {
 			// Iterate over the cell elements
 			smplxIt.begin();
 			for (size_t j = 0; j < numElems; ++j) {
 				const std::vector<size_t>& inds = smplxIt();
 				for (int k = 0; k < order + 1; ++k) {
 					size_t index = (order + 1)*(i*numElems + j) + k;
 					ptConnect[index] = ptIds[inds[k]];
 				}
 				smplxIt.next();
 			}
 		}
 		else {
 			std::cerr << "ERROR: cell is not a tetrahedron!\n";
 			res = 1;
 		}
 	}
 	return res;
 }


} // extern "C"
