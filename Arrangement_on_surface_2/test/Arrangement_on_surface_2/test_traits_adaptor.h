#ifndef CGAL_TEST_TRAITS_ADAPTOR_H
#define CGAL_TEST_TRAITS_ADAPTOR_H

#include "test_kernel.h"
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

// ============================================================================
// Traits includes:
// ============================================================================
#if TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS
#include <CGAL/Arr_linear_traits_2.h>

#elif TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#else
#error No geometry traits adaptor (GEOM_TRAITS) specified!
#endif

// ============================================================================
// Traits typedef:
// ============================================================================
#if TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>                 Traits_parameter;
#define GEOM_TRAITS_TYPE "Segments"

#elif TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS
typedef CGAL::Arr_linear_traits_2<Kernel>                  Traits_parameter;
#define GEOM_TRAITS_TYPE "Linear Lines"

#elif TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>  Traits_parameter;
#define GEOM_TRAITS_TYPE "Geodesic Arc on Sphere"

#else
#error No traits (TRAITS) specified!
#endif

typedef CGAL::Arr_traits_basic_adaptor_2<Traits_parameter> Geom_traits;

#endif
