
#include "egfFieldType.hpp"
#include <vtkFieldData.h>
#include <iostream>
#include <sstream>

extern "C" {

int egfField_new(egfFieldType** self) {
    *self = new egfFieldType();
    (*self)->ugrid = NULL;
    (*self)->order = -1;
    (*self)->type = "UNKNOWN";
    return 0;
}

int egfField_del(egfFieldType** self) {
    for (size_t i = 0; i < (*self)->dataArrays.size(); ++i) {
        (*self)->dataArrays[i]->Delete();
    }
    delete *self;
    return 0;
}

int egfField_print(egfFieldType** self) {
    std::cout << "egfField\n";
    if ((*self)->ugrid && (*self)->order >= 0) {
        size_t numElems = (*self)->elems.size();
        std::cout << "type: " << (*self)->type << " order: " << (*self)->order 
        << " number of elements: " << numElems << '\n';
        for (size_t i = 0; i < numElems; ++i) {
            std::vector<size_t>& elm = (*self)->elems[i];
            std::cout << "element " << i << ": ";
            for (size_t j = 0; j < (*self)->order + 1; ++j) {
                std::cout << elm[j] << ' ';
            }
            std::cout << '\n';
        }
    }
    return 0;
}
   
int egfField_setOrder(egfFieldType** self, int order) {
    (*self)->order = order;
    if ((*self)->order == 0) {
        (*self)->type = "nodal";
    }
    else if ((*self)->order == 1) {
        (*self)->type = "edge";
    }
    else if ((*self)->order == 2) {
        (*self)->type = "face";
    }
    else if ((*self)->order == 3) {
        (*self)->type = "cell";
    }

    (*self)->smplxIt.setup(3, order); // Always in 3D
    size_t numElems = (*self)->smplxIt.getNumberOfElements();
    (*self)->smplxIt.begin();
    for (int i = 0; i < numElems; ++i) {
        const std::vector<size_t>& es = (*self)->smplxIt();
        (*self)->elems.insert(std::pair<int, std::vector<size_t> >(i, es));
        (*self)->smplxIt.next();
    }
    (*self)->ugrid->GetFieldData()->AllocateArrays(numElems);
    (*self)->dataArrays.resize(numElems);

    return 0;
}

int egfField_setGrid(egfFieldType** self, egfGridType* grid) {
    (*self)->ugrid = grid->ugrid;
    return 0;
}
    
int egfField_getNumberOfElements(egfFieldType** self, int* numElems) {
    *numElems = (int) (*self)->smplxIt.getNumberOfElements();
    return 0;
}
    
int egfField_getElement(egfFieldType** self, int elem, int inds[]) {
    for (size_t j = 0; j < (*self)->elems[elem].size(); ++j) {
        inds[j] = (*self)->elems[elem][j];
    }
    return 0;
}
    
int egfField_setValues(egfFieldType** self, int elem, const double vals[]) {
    vtkDoubleArray* data = vtkDoubleArray::New();
    vtkIdType numCells = (*self)->ugrid->GetNumberOfCells();
    size_t numElems = (*self)->smplxIt.getNumberOfElements();
    data->SetNumberOfComponents(numElems);
    std::map<int, std::vector<size_t> >::const_iterator it = (*self)->elems.find(elem);
    std::ostringstream ss;
    ss << elem << ": ";
    for (size_t i = 0; i < it->second.size(); ++i) {
        ss << it->second[i] << ' ';
    }
    data->SetName(ss.str().c_str());

    // Caller owns data, vals' life expectancy must be longer than this instance's
    data->SetVoidArray((double*) vals, numCells * numElems, 1);

    (*self)->dataArrays[elem] = data;

    (*self)->ugrid->GetFieldData()->AddArray(data);
    return 0;
}

} // extern "C"
