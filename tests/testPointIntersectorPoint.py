import argparse
import os
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int
import numpy
import sys

parser = argparse.ArgumentParser(description='Test point intersector.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
parser.add_argument('--input', dest='input', default='',
                    help='Specify input file')
parser.add_argument('--p0', dest='p0', type=str, default='0,0,0',
                    help='Specify segment starting point')
args = parser.parse_args()

suffix = 'so'
if os.uname()[0] == 'Darwin':
    suffix = 'dylib'
elif os.uname()[0] == 'Windows':
    suffix = 'dll'

libName = 'libegface'
# Open the shared library
lib = cdll.LoadLibrary(args.build_dir + '/' + libName + '.' + suffix)

# Load a grid from file
grid = c_void_p(0)
ier = lib.egfUnstructuredGrid3d_new(byref(grid))
ier = lib.egfUnstructuredGrid3d_loadFromFile(byref(grid), args.input)

# Opaque handle
intrsctr = c_void_p(0)

# Constructor
ier = lib.egfPointIntersector_new(byref(intrsctr))
assert ier == 0

# Set the tolerance
tol = c_double(1.e-12)
ier = lib.egfPointIntersector_setTolerance(byref(intrsctr), tol)
assert ier == 0

# Set the number of cells per bucket
ier = lib.egfPointIntersector_setNumberOfCellsPerBucket(byref(intrsctr), 100)
assert ier == 0

# Set the grid, should be called after setNumberOfCellsPerBucket
ier = lib.egfPointIntersector_setGrid(byref(intrsctr), grid)
assert ier == 0

# Compute all the intersection points
p0 = numpy.array(eval(args.p0), numpy.float64)
print 'p0 = ', p0
ier = lib.egfPointIntersector_gridWithPoint(byref(intrsctr),
                                           p0.ctypes.data_as(POINTER(c_double)))
assert ier == 0

# Print the object
ier = lib.egfPointIntersector_print(byref(intrsctr))
assert ier == 0

# Run checks
numCells = c_int()
ier = lib.egfPointIntersector_getNumberOfCells(byref(intrsctr), byref(numCells))
assert ier == 0

cellIds = numpy.zeros((numCells.value,), numpy.int32)
ier = lib.egfPointIntersector_fillInCellIds(byref(intrsctr), cellIds.ctypes.data_as(POINTER(c_int)))
assert ier == 0

numPoints = c_int()
for cellId in cellIds:
    ci = c_int(cellId)
    ier = lib.egfPointIntersector_getNumberOfPointsInCell(byref(intrsctr), ci, byref(numPoints))
    assert ier == 0
    print 'cell ', cellId, ' has ', numPoints.value, 'intersection points'
    points = numpy.zeros((numPoints.value*3,), numpy.float64)
    ier = lib.egfPointIntersector_fillInPointsInCell(byref(intrsctr), ci, 
                                                     points.ctypes.data_as(POINTER(c_double)))
    assert ier == 0

# Destructor
ier = lib.egfPointIntersector_del(byref(intrsctr))
assert ier == 0

# Delete grid
ier = lib.egfUnstructuredGrid3d_del(byref(grid))
