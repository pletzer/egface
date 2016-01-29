#include <egfGrid.hpp>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

extern "C" {

int egfGrid_new(vtkUnstructuredGrid** self) {
	(*self) = new egfGridType();
}

int egfGrid_del(vtkUnstructuredGrid** self) {
	delete *self);
}

int egfGrid_loadFromFile(vtkUnstructuredGrid** self, const char* filename) {
	vtkUnstructuredGridReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(filename);
    reader->Update();
    (*self)->ugrid = reader->GetOutputData();
	reader->Delete();
}

int egfGrid_saveToFile(vtkUnstructuredGrid** self, const char* filename);
	vtkUnstructuredGridWriter* writer = vtkPolyDataWriter::New();
    writer->SetFileName(filename);
    writer->SetInputData((*self)->ugrid);
    writer->Update();
	writer->Delete();
}
