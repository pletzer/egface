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
locator = c_void_p(0)

# Constructor
ier = lib.egfPointIntersector_new(byref(locator))
assert ier == 0

# Set the grid
ier = lib.egfPointIntersector_setGrid(byref(locator), grid)
assert ier == 0

# Set the tolerance
lib.egfPointIntersector_setTolerance.argtypes = [c_void_p, c_double]
ier = lib.egfPointIntersector_setTolerance(byref(locator), 1.e-12)
assert ier == 0

# Compute all the intersection points
p0 = numpy.array(eval(args.p0))
p1 = numpy.array(eval(args.p1))
ier = lib.egfPointIntersector_gridWithLine(byref(locator),
                                           p0.ctypes.data_as(POINTER(c_double)),
                                           p1.ctypes.data_as(POINTER(c_double)))
assert ier == 0

# Print the object
ier = lib.egfPointIntersector_print(byref(locator))
assert ier == 0

# Destructor
ier = lib.egfCellLocator_del(byref(locator))
assert ier == 0

# Delete grid
ier = lib.egfGrid_del(byref(grid))
