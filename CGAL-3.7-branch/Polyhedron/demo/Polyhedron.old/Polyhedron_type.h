#ifndef POLYHEDRON_TYPE_H
#define POLYHEDRON_TYPE_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#ifdef USE_FORWARD_DECL
struct Kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#endif

// surface mesh
#include <CGAL/Polyhedron_3.h>

#include "Polyhedron_type_fwd.h"

typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

#endif // POLYHEDRON_TYPE_H
