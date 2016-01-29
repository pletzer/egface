#include "egfGridType.hpp"
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <iostream>

extern "C" {

int egfGrid_new(egfGridType** self) {
	(*self) = new egfGridType();
	(*self)->reader = vtkUnstructuredGridReader::New();
	(*self)->ugrid = vtkUnstructuredGrid::New();
	return 0;
}

int egfGrid_del(egfGridType** self) {
	(*self)->ugrid->Delete();
	(*self)->reader->Delete();
	delete *self;
    return 0;
}

int egfGrid_print(egfGridType** self) {
	std::cout << "egfGrid:\n";
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

int egfGrid_loadFromFile(egfGridType** self, const char* filename) {
    (*self)->reader->SetFileName(filename);
    (*self)->reader->Update();
#if (VTK_MAJOR_VERSION < 6)
    (*self)->ugrid = (*self)->reader->GetOutput();
#else
    (*self)->ugrid = (*self)->reader->GetOutputData();
#endif
	return 0;
}

int egfGrid_saveToFile(egfGridType** self, const char* filename) {
	vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    writer->SetFileName(filename);
    writer->Update();
#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput((*self)->ugrid);
#else
    writer->SetInputData((*self)->ugrid);
#endif
	writer->Delete();
	return 0;
}

 int egfGrid_getNumberOfCells(egfGridType** self, int* numCells) {
 	*numCells = (int) (*self)->ugrid->GetNumberOfCells();
 	return 0;
 }

 int egfGrid_getNumberOfPoints(egfGridType** self, int* numPoints) {
 	*numPoints = (int) (*self)->ugrid->GetNumberOfPoints();
 	return 0; 	
 }

} // extern "C"
