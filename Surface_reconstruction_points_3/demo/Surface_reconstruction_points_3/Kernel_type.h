#ifndef KERNEL_TYPE_H
#define KERNEL_TYPE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// (Main) CGAL kernel used by this demo
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid;
typedef Kernel::Plane_3 Plane_3;

#endif // KERNEL_TYPE_H
