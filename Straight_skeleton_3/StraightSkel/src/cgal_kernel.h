/**
 * @file   cgal_kernel.h
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#ifndef CGAL_KERNEL_H
#define CGAL_KERNEL_H

#include "config.h"

#include <CGAL/double.h>
#include <CGAL/float.h>
#include <CGAL/number_utils.h> // to_double, to_interval, approximate_sqrt

#ifdef USE_CGAL

// #include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

namespace CGAL {

// typedef Cartesian<double> K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt EPECK_w_sqrt;

typedef EPICK K;

typedef K::FT             FT;

typedef K::Point_2        Point2;
typedef K::Vector_2       Vector2;
typedef K::Line_2         Line2;
typedef K::Segment_2      Segment2;

typedef K::Point_3        Point3;
typedef K::Vector_3       Vector3;
typedef K::Line_3         Line3;
typedef K::Segment_3      Segment3;
typedef K::Plane_3        Plane3;
typedef K::Sphere_3       Sphere3;

} // namespace CGAL

#else /* USE_CGAL */

namespace CGAL {

typedef double            FT;

} // namespace CGAL

#endif  /* USE_CGAL */

#endif /* CGAL_KERNEL_H */
