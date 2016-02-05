import vtk
import numpy
import argparse

parser = argparse.ArgumentParser(description='Generate cube grid.')
parser.add_argument('--nx', dest='nx', type=int, default=1,
                    help='Specify the number of cells in the x direction')
parser.add_argument('--ny', dest='ny', type=int, default=2,
                    help='Specify the number of cells in the y direction')
parser.add_argument('--nz', dest='nz', type=int, default=3,
                    help='Specify the number of cells in the z direction')
parser.add_argument('--output', dest='output', type=str, default='cube.vtk',
                    help='Specify the output file')
args = parser.parse_args()

xmin, xmax = 0., 1.
ymin, ymax = 0., 1.
zmin, zmax = 0., 1.

assert args.nx > 0 and xmax > xmin
assert args.ny > 0 and ymax > ymin
assert args.nz > 0 and zmax > zmin
nx, ny, nz = args.nx, args.ny, args.nz
nx1, ny1 , nz1 = nx + 1, ny + 1, nz + 1

dx = (xmax - xmin)/float(nx)
dy = (ymax - ymin)/float(ny)
dz = (zmax - zmin)/float(nz)

# Create points
points = vtk.vtkPoints()
numPoints = nx1*ny1*nz1
points.SetNumberOfPoints(numPoints)
index = 0
for k in range(nz1):
    z = zmin + k*dz
    for j in range(ny1):
        y = ymin + j*dy
        for i in range(nx1):
            x = xmin + i*dx
            points.SetPoint(index, x, y, z)
            index += 1

# Create grid
ug = vtk.vtkUnstructuredGrid()
ug.SetPoints(points)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(8)
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            index = i + j*nx1 + k*nx1*ny1
            ptIds.SetId(0, index + 0 + 0*nx1 + 0*nx1*ny1)
            ptIds.SetId(1, index + 1 + 0*nx1 + 0*nx1*ny1)
            ptIds.SetId(2, index + 1 + 1*nx1 + 0*nx1*ny1)
            ptIds.SetId(3, index + 0 + 1*nx1 + 0*nx1*ny1)
            ptIds.SetId(4, index + 0 + 0*nx1 + 1*nx1*ny1)
            ptIds.SetId(5, index + 1 + 0*nx1 + 1*nx1*ny1)
            ptIds.SetId(6, index + 1 + 1*nx1 + 1*nx1*ny1)
            ptIds.SetId(7, index + 0 + 1*nx1 + 1*nx1*ny1)
            ug.InsertNextCell(vtk.VTK_HEXAHEDRON, ptIds)

# Save grid to file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(args.output)
writer.SetInput(ug)
writer.Update()


