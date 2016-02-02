import argparse
import os
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long, c_int

parser = argparse.ArgumentParser(description='Test cell locator.')

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
domMin = (3*c_double)(0, 0, 0)
domMax = (3*c_double)(0, 0, 0)
ier = lib.egfGrid_getDomainBounds(byref(grid), domMin, domMax)
posMin = (3*c_double)(0, 0, 0)
posMax = (3*c_double)(0, 0, 0)
for i in range(3):
    posMin[i] = domMin[i] + 0.55*(domMax[i] - domMin[i])
    posMax[i] = domMin[i] + 0.65*(domMax[i] - domMin[i])
    print '{} min/max domain: {}/{} box: {}/{}'.format(i, 
        domMin[i], domMax[i], posMin[i], posMax[i])
ier = lib.egfCellLocator_findCellsInBox(byref(locator), posMin, posMax)
assert ier == 0

# Get the number of found cells
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
