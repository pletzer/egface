
#include "egfUnstructuredField3dType.hpp"
#include <vtkFieldData.h>
#include <iostream>
#include <sstream>

extern "C" {

int egfUnstructuredField3d_new(egfUnstructuredField3dType** self) {
    *self = new egfUnstructuredField3dType();
    (*self)->ugrid = NULL;
    (*self)->order = -1;
    (*self)->type = "UNKNOWN";
    return 0;
}

int egfUnstructuredField3d_del(egfUnstructuredField3dType** self) {
    for (size_t i = 0; i < (*self)->dataArrays.size(); ++i) {
        if ((*self)->dataArrays[i]) {
            (*self)->dataArrays[i]->Delete();
        }
    }
    delete *self;
    return 0;
}

int egfUnstructuredField3d_print(egfUnstructuredField3dType** self) {
    std::cout << "egfUnstructuredField3d\n";
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
   
int egfUnstructuredField3d_setOrder(egfUnstructuredField3dType** self, int order) {
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
    (*self)->dataArrays.resize(numElems, NULL);

    return 0;
}

int egfUnstructuredField3d_setGrid(egfUnstructuredField3dType** self, egfUnstructuredGrid3dType* grid) {
    (*self)->ugrid = grid->ugrid;
    return 0;
}
    
int egfUnstructuredField3d_getNumberOfElements(egfUnstructuredField3dType** self, int* numElems) {
    *numElems = (int) (*self)->smplxIt.getNumberOfElements();
    return 0;
}
    
int egfUnstructuredField3d_getElement(egfUnstructuredField3dType** self, int elem, int inds[]) {
    for (size_t j = 0; j < (*self)->elems[elem].size(); ++j) {
        inds[j] = (*self)->elems[elem][j];
    }
    return 0;
}
    
int egfUnstructuredField3d_setValues(egfUnstructuredField3dType** self, int elem, const double vals[]) {
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
