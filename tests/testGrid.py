import argparse
import os
import sys
import numpy
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int

parser = argparse.ArgumentParser(description='Test grid.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
parser.add_argument('--input', dest='input', default='',
                    help='Specify input file')
parser.add_argument('--output', dest='output', default='',
                    help='Specify output file')
args = parser.parse_args()

suffix = 'so'
if os.uname()[0] == 'Darwin':
    suffix = 'dylib'
elif os.uname()[0] == 'Windows':
    suffix = 'dll'

libName = 'libegface'
# Open the shared library
lib = cdll.LoadLibrary(args.build_dir + '/' + libName + '.' + suffix)

# Opaque handle
handle = c_void_p(0)

# Constructor
ier = lib.egfGrid_new(byref(handle))
assert ier == 0

# Load from file
ier = lib.egfGrid_loadFromFile(byref(handle), args.input)
assert ier == 0

# Check number of points
numPoints = c_int()
ier = lib.egfGrid_getNumberOfPoints(byref(handle), byref(numPoints))
assert ier == 0
print 'numPoints = ', numPoints.value

# Check number of cells
numCells = c_int()
ier = lib.egfGrid_getNumberOfCells(byref(handle), byref(numCells))
assert ier == 0
print 'numCells = ', numCells.value

# Get the connectivity
for order in range(0, 4):
	numElems = c_int(-1)
	ier = lib.egfGrid_getNumberOfElements(byref(handle), order, byref(numElems))
	assert ier == 0
	print >> sys.stderr, 'order = ', order, ' numElems = ', numElems.value
	ptConnect = numpy.zeros((numCells.value, numElems.value, order + 1), numpy.int32)
	ier = lib.egfGrid_getElementConnectivity(byref(handle), order, 
		ptConnect.ctypes.data_as(POINTER(c_double)))
	assert ier == 0

# Print
ier = lib.egfGrid_print(byref(handle));
assert ier == 0

# Save to file
ier = lib.egfGrid_saveToFile(byref(handle), args.output)
assert ier == 0

# Destructor
ier = lib.egfGrid_del(byref(handle))
assert ier == 0