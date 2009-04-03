#ifndef POINT_SET_DEMO_TYPES_H
#define POINT_SET_DEMO_TYPES_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// surface mesh
#include <CGAL/Polyhedron_3.h>

#include "Point_set_demo_types_fwd.h"

typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Plane_3 Plane_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

#endif // POINT_SET_DEMO_TYPES_H
