#ifndef CGAL_TEST_TRAITS_H
#define CGAL_TEST_TRAITS_H

#include "test_kernel.h"

// ============================================================================
// Traits includes:
// ============================================================================
#if TEST_TRAITS == SEGMENT_TRAITS
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_TRAITS == NON_CACHING_SEGMENT_TRAITS
#include <CGAL/Arr_non_caching_segment_traits_2.h>

#elif TEST_TRAITS == POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>

#elif TEST_TRAITS == CORE_CONIC_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>

#else
#error No traits (TRAITS) specified!
#endif

// ============================================================================
// Traits typedef:
// ============================================================================
#if TEST_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define TRAITS_TYPE "Segments"

#elif TEST_TRAITS == NON_CACHING_SEGMENT_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Traits;
#define TRAITS_TYPE "Non Caching Segments"

#elif TEST_TRAITS == POLYLINE_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Polylines"

#elif TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Non Caching Polylines"

#elif TEST_TRAITS == CORE_CONIC_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Traits;
#define TRAITS_TYPE "Conics"

#else
#error No traits (TRAITS) specified!
#endif

#endif
