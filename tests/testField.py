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
args = parser.parse_args()

suffix = 'so'
if os.uname()[0] == 'Darwin':
    suffix = 'dylib'
elif os.uname()[0] == 'Windows':
    suffix = 'dll'

libName = 'libegface'
# Open the shared library
lib = cdll.LoadLibrary(args.build_dir + '/' + libName + '.' + suffix)

print >> sys.stderr, '1'
# Opaque handle
field = c_void_p(0)

print >> sys.stderr, '2'
# Constructor
ier = lib.egfField_new(byref(field))
assert ier == 0

print >> sys.stderr, '3'
# Set the staggering
ier = lib.egfField_setOrder(byref(field), args.order)

print >> sys.stderr, '4'
# Load the grid from file
grid = c_void_p(0)
ier = lib.egfGrid_new(byref(grid))
ier = lib.egfGrid_loadFromFile(byref(grid), args.input)
assert ier == 0

print >> sys.stderr, '5'
# Set the grid
ier = lib.egfField_setGrid(byref(field), grid)

print >> sys.stderr, '6'
ier = lib.egfGrid_del(byref(grid))
assert ier == 0

print >> sys.stderr, '7'
# Destructor
ier = lib.egfField_del(byref(field))
assert ier == 0

print >> sys.stderr, '8'