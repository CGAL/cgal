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

#elif TEST_TRAITS == LINEAR_TRAITS
#include <CGAL/Arr_linear_traits_2.h>

#elif TEST_TRAITS == CORE_CONIC_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>

#elif TEST_TRAITS == LINE_ARC_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_line_arc_traits_2.h>

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS
#include <CGAL/Arr_circle_segment_traits_2.h>

#elif TEST_TRAITS == BEZIER_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>

#elif TEST_TRAITS == SPHERICAL_ARC_TRAITS
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#elif TEST_TRAITS == RATIONAL_ARC_TRAITS
//#include <CGAL/CORE_algebraic_number_traits.h>  //OS - new
//#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>


#elif TEST_TRAITS == ALGEBRAIC_TRAITS
#include <CGAL/Arr_algebraic_segment_traits_2.h>

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
typedef Traits::Point_2                                 Point_2;
typedef Traits::Rat_point_2                             Rat_point;
typedef Traits::Rat_circle_2                            Rat_circle;
typedef Traits::Rat_segment_2                           Rat_segment;
#define TRAITS_TYPE "Conics"

#elif TEST_TRAITS == LINEAR_TRAITS
typedef CGAL::Arr_linear_traits_2<Kernel>               Traits;
#define TRAITS_TYPE "Linear Lines"

#elif TEST_TRAITS == LINE_ARC_TRAITS
typedef Kernel                                                Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>   Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                              Circular_kernel;
typedef Circular_kernel::Line_arc_2                           Line_arc_2;
typedef CGAL::Arr_line_arc_traits_2<Circular_kernel>          Traits;

#define TRAITS_TYPE "Line Arc"

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS
typedef Kernel                                                Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>   Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                              Circular_kernel;
typedef CGAL::Arr_circular_arc_traits_2<Circular_kernel>      Traits;

#define TRAITS_TYPE "Circular Arc"

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

typedef Kernel                                                Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>   Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                              Circular_kernel;
typedef Circular_kernel::Circular_arc_2                       Circular_arc_2;
typedef Circular_kernel::Line_arc_2                           Line_arc_2;
typedef CGAL::Arr_circular_line_arc_traits_2<Circular_kernel> Traits;

#define TRAITS_TYPE "Circular Line Arc"

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS
typedef CGAL::Cartesian<Number_type>                          Rat_kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>             Traits;
typedef Rat_kernel::FT                                        Rat_nt;
typedef Rat_kernel::Circle_2                                  Circle_2;
typedef Rat_kernel::Line_2                                    Line_2;
typedef Rat_kernel::Segment_2                                 Segment_2;
typedef Rat_kernel::Point_2                                   Rat_point_2;
typedef Traits::Point_2                                       Point_2;
#define TRAITS_TYPE "Circle Segments"

#elif TEST_TRAITS == BEZIER_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits;
typedef Traits::Bezier_cache                            Bezier_cache;

#define TRAITS_TYPE "Bezier"

#elif TEST_TRAITS == SPHERICAL_ARC_TRAITS
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>        Traits;
#define TRAITS_TYPE "Spherical Arc"

#elif TEST_TRAITS == RATIONAL_ARC_TRAITS

// TODO: This is for old Ron's traits---remove.
//typedef CGAL::CORE_algebraic_number_traits                 Nt_traits; //OS - new
//typedef Nt_traits::Rational                                Rational;
//typedef Nt_traits::Algebraic                               Algebraic;
//typedef CGAL::Arr_rational_arc_traits_2<Kernel,Nt_traits>  Traits;
//typedef Traits::Rat_vector                                 Rat_vector;
//typedef Traits::Point_2                                    Point_2;
//#define TRAITS_TYPE "Rational Arc"

typedef CGAL::Arr_rational_function_traits_2<Kernel>	     Traits;
typedef Traits::Rational                                     Rational;
typedef Traits::Algebraic_real_1                             Algebraic_real_1;
typedef Traits::Point_2                                      Point_2;
typedef Traits::Rat_vector                                   Rat_vector;
#define TRAITS_TYPE "Rational Arc"


#elif TEST_TRAITS == ALGEBRAIC_TRAITS
typedef CGAL::Arr_algebraic_segment_traits_2<Number_type>  Traits;
typedef Traits::Point_2                                    Point_2;
typedef Traits::Curve_2                                    Curve_2;
typedef Traits::X_monotone_curve_2                         X_monotone_curve_2;

#define TRAITS_TYPE "Algebraic"


#else
#error No traits (TRAITS) specified!
#endif

#endif
