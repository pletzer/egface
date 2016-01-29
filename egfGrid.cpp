#include "egfGridType.hpp"
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

extern "C" {

int egfGrid_new(egfGridType** self) {
	(*self) = new egfGridType();
	return 0;
}

int egfGrid_del(egfGridType** self) {
	delete *self;
    return 0;
}

int egfGrid_loadFromFile(egfGridType** self, const char* filename) {
	vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
    reader->SetFileName(filename);
    reader->Update();
#if (VTK_MAJOR_VERSION < 6)
    (*self)->ugrid = reader->GetOutput();
#else
    (*self)->ugrid = reader->GetOutputData();
#endif
	reader->Delete();
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

} // extern "C"
