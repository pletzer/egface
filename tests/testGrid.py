import argparse
import os
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long

parser = argparse.ArgumentParser(description='Translate shape.')

parser.add_argument('--build_dir', dest='build_dir', default='./',
                    help='Specify build directory')
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

# Destructor
ier = lib.egfGrid_del(byref(handle))
assert ier == 0