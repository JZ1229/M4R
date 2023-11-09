from firedrake import *

def ImmersedRectangularMesh(nx, ny, x0, x1, r, quadrilateral=False):

    m = UnitSquareMesh(nx, ny, quadrilateral=quadrilateral)

    coord_family = "DQ" if quadrilateral else "DG"

    cell = "quadrilateral" if quadrilateral else "triangle"

    coord_fs = VectorFunctionSpace(m, FiniteElement(coord_family, cell, 1, variant="equispaced"), dim=3)

    x_old, y_old = SpatialCoordinate(m)

    longitude = x0[0] + x_old*(x1[0]- x0[0])
    latitude = x0[1] + y_old*(x1[1]- x0[1])

    x = r*cos(latitude)*cos(longitude)
    y = r*cos(latitude)*sin(longitude)    
    z = r*sin(latitude)

    new_coordinates = Function(coord_fs).interpolate(as_vector([x,y,z]))

    return Mesh(new_coordinates)

