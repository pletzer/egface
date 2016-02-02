#include "egfGridType.hpp"
#include "egfCellLocatorType.hpp"
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <iostream>

extern "C" {

int egfCellLocator_new(egfCellLocatorType** self) {
	(*self) = new egfCellLocatorType();
	(*self)->cellLocator = vtkCellLocator::New();
	(*self)->treeHasBeenBuilt = false;
	return 0;
}

int egfCellLocator_del(egfCellLocatorType** self) {
	(*self)->cellLocator->Delete();
	delete *self;
    return 0;
}

int egfCellLocator_print(egfCellLocatorType** self) {
	std::cout << "egfCellLocator:\n";
	std::cout << "found cell Ids: ";
	for (size_t i = 0; i < (*self)->cellIds.size(); ++i) {
		std::cout << (*self)->cellIds[i] << ' ';
	}
    std::cout << '\n';
    return 0;
}

int egfCellLocator_setGrid(egfCellLocatorType** self, egfGridType* grid) {
    (*self)->cellLocator->SetDataSet(grid->ugrid);
    return 0;
}

 int egfCellLocator_setNumberOfCellsPerBucket(egfCellLocatorType** self, int numCellsPerBucket) {
    (*self)->cellLocator->SetNumberOfCellsPerBucket(numCellsPerBucket);
    return 0;
 }

 int egfCellLocator_findCellsInBox(egfCellLocatorType** self, const double posMin[], const double posMax[]) {
    if (! (*self)->treeHasBeenBuilt) {
        (*self)->cellLocator->BuildLocator();
        (*self)->treeHasBeenBuilt = true;
    }
    double bbox[] = {posMin[0], posMax[0], posMin[1], posMax[1], posMin[2], posMax[2]};
    vtkIdList* cellIdList = vtkIdList::New();
    (*self)->cellLocator->FindCellsWithinBounds(bbox, cellIdList);
    size_t numCells = (size_t) cellIdList->GetNumberOfIds();
    (*self)->cellIds.resize((size_t) numCells);
    for (size_t i = 0; i < (size_t) numCells; ++i) {
        (*self)->cellIds[i] = (int) cellIdList->GetId(i);
    }
    cellIdList->Delete();
    return 0;
 }

 int egfCellLocator_getNumberOfCells(egfCellLocatorType** self, int* numCells) {
    *numCells = (int) (*self)->cellIds.size();
    return 0;
 }

 int egfCellLocator_fillInCellIds(egfCellLocatorType** self, int* cellIds) {
    for (size_t i = 0; i < (*self)->cellIds.size(); ++i) {
        cellIds[i] = (*self)->cellIds[i];
    }
    return 0;
 } 

 
} // extern "C"
