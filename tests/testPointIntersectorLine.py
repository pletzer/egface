import argparse
import os
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int
import numpy

parser = argparse.ArgumentParser(description='Test point intersector.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
parser.add_argument('--input', dest='input', default='',
                    help='Specify input file')
parser.add_argument('--p0', dest='p0', type=str, default='0,0,0',
                    help='Specify segment starting point')
parser.add_argument('--p1', dest='p1', type=str, default='1,0,0',
                    help='Specify segment ending point')
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
ier = lib.egfGrid_new(byref(grid))
ier = lib.egfGrid_loadFromFile(byref(grid), args.input)

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
p1 = numpy.array(eval(args.p1), numpy.float64)
print 'p0 = ', p0, ' p1 = ', p1
ier = lib.egfPointIntersector_gridWithLine(byref(intrsctr),
                                           p0.ctypes.data_as(POINTER(c_double)),
                                           p1.ctypes.data_as(POINTER(c_double)))
assert ier == 0

# Print the object
ier = lib.egfPointIntersector_print(byref(intrsctr))
assert ier == 0

# Check
numPoints = c_int()
ier = lib.egfPointIntersector_getNumberOfPoints(byref(intrsctr), byref(numPoints))
assert ier == 0
print 'number of intersection points: ', numPoints.value

# Fill in the intersection points
points = numpy.zeros((3*numPoints.value,), numpy.float64)
ier = lib.egfPointIntersector_fillInPoints(byref(intrsctr), points.ctypes.data_as(POINTER(c_double)))
assert ier == 0
print 'intersection points: \n', points.reshape((numPoints.value, 3))

# Destructor
ier = lib.egfCellLocator_del(byref(intrsctr))
assert ier == 0

# Delete grid
ier = lib.egfGrid_del(byref(grid))