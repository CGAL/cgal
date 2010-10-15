#ifndef CGAL_TEST_TRAITS_H
#define CGAL_TEST_TRAITS_H

#include "test_kernel.h"
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

// ============================================================================
// Traits includes:
// ============================================================================
#if TEST_TRAITS == SEGMENT_TRAITS
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_TRAITS == LINEAR_TRAITS
#include <CGAL/Arr_linear_traits_2.h>

#elif TEST_TRAITS == SPHERICAL_ARC_TRAITS
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#else
#error No traits (TRAITS) specified!
#endif

// ============================================================================
// Traits typedef:
// ============================================================================
#if TEST_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_parameter;
#define TRAITS_TYPE "Segments"

#elif TEST_TRAITS == LINEAR_TRAITS
typedef CGAL::Arr_linear_traits_2<Kernel>               Traits_parameter;
#define TRAITS_TYPE "Linear Lines"

#elif TEST_TRAITS == SPHERICAL_ARC_TRAITS
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>        Traits_parameter;
#define TRAITS_TYPE "Spherical Arc"

#else
#error No traits (TRAITS) specified!
#endif

typedef CGAL::Arr_traits_basic_adaptor_2 < Traits_parameter >    Traits;

#endif
