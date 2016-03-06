import argparse
import os
import sys
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int

parser = argparse.ArgumentParser(description='Test field.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
parser.add_argument('--input', dest='input', default='',
                    help='Specify input grid file')
parser.add_argument('--order', dest='order', default=0, type=int,
	                help='Specify order (0=nodal, 1=edge, 2=face, 3=cell)')
parser.add_argument('--expression', dest='expression',
	default='dx ^ dy + 2 * dy ^ dz',
	help='Specify form as an expression of x, y, and z, eg dx ^ dy + 2 * dy ^ dz for a 2-form')
args = parser.parse_args()

suffix = 'so'
if os.uname()[0] == 'Darwin':
    suffix = 'dylib'
elif os.uname()[0] == 'Windows':
    suffix = 'dll'

libName = 'libegface'
# Open the shared library
lib = cdll.LoadLibrary(args.build_dir + '/' + libName + '.' + suffix)

# Load the grid from file
grid = c_void_p(0)
ier = lib.egfUnstructuredGrid3d_new(byref(grid))
ier = lib.egfUnstructuredGrid3d_loadFromFile(byref(grid), args.input)
numCells = c_int()
ier = lib.egfUnstructuredGrid3d_getNumberOfCells(byref(grid), byref(numCells))
assert ier == 0
print 'number of grid cells: ', numCells.value

# Opaque handle
field = c_void_p(0)

# Constructor
ier = lib.egfUnstructuredField3d_new(byref(field))
assert ier == 0

# Set the grid
ier = lib.egfUnstructuredField3d_setGrid(byref(field), grid)
assert ier == 0

# Set the staggering
ier = lib.egfUnstructuredField3d_setOrder(byref(field), args.order)
assert ier == 0

# Get the number of elements per cell
numElems = c_int()
ier = lib.egfUnstructuredField3d_getNumberOfElements(byref(field), byref(numElems))
assert ier == 0
print 'number of elements per cell: ', numElems.value

np1 = args.order + 1
inds = (c_int * np1)()
for i in range(numElems.value):
	ier = lib.egfUnstructuredField3d_getElement(byref(field), i, inds)
	assert ier == 0
	print 'element i = ', i, ' inds = ', inds[:]

# Check print
ier = lib.egfUnstructuredField3d_print(byref(field))
assert ier == 0

# Set the form



# Destroy field
ier = lib.egfUnstructuredField3d_del(byref(field))
assert ier == 0

# Destroy grid
ier = lib.egfUnstructuredGrid3d_del(byref(grid))
assert ier == 0

