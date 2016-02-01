import argparse
import os
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int

parser = argparse.ArgumentParser(description='Translate shape.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
parser.add_argument('--input', dest='input', default='',
                    help='Specify input file')
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
ier = lib.egfCellLocator_new(byref(locator))
assert ier == 0

# Set the grid
ier = lib.egfCellLocator_setGrid(byref(locator), grid)
assert ier == 0

# Set the number of cells per bucket
ier = lib.egfCellLocator_setNumberOfCellsPerBucket(byref(locator), 20)
assert ier == 0

# Find the cell within a box
bbox = (6*c_double)(-10, 10, -5, 5, -2, 2)
ier = lib.egfCellLocator_findCellsInBox(byref(locator), bbox)
assert ier == 0

# Get the number of ofund cells
numCells = c_int(0)
ier = lib.egfCellLocator_getNumberOfCells(byref(locator), byref(numCells))
assert ier == 0

# Fill in the cell IDs
cellIds = (numCells.value*c_int)()
ier = lib.egfCellLocator_fillInCellIds(byref(locator), cellIds)
assert ier == 0

for i in range(numCells.value):
	print '{} cell ID: {}'.format(i, cellIds[i])

# Destructor
ier = lib.egfCellLocator_del(byref(locator))
assert ier == 0

# Delete grid
ier = lib.egfGrid_del(byref(grid))
