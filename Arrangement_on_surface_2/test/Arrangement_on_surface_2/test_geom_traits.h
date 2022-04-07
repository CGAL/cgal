#ifndef CGAL_TEST_GEOM_TRAITS_H
#define CGAL_TEST_GEOM_TRAITS_H

#include "test_kernel.h"

// ============================================================================
// Geometry Traits includes:
// ============================================================================
#if TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == NON_CACHING_SEGMENT_GEOM_TRAITS
#include <CGAL/Arr_non_caching_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS
#include <CGAL/Arr_linear_traits_2.h>

#elif TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>

#elif TEST_GEOM_TRAITS == LINE_ARC_GEOM_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_line_arc_traits_2.h>

#elif TEST_GEOM_TRAITS == CIRCULAR_ARC_GEOM_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>

#elif TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>

#elif TEST_GEOM_TRAITS == CIRCLE_SEGMENT_GEOM_TRAITS
#include <CGAL/Arr_circle_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>

#elif TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#elif TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS
//#include <CGAL/CORE_algebraic_number_traits.h>  //OS - new
//#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>


#elif TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == FLAT_TORUS_GEOM_TRAITS
#include <CGAL/Arr_flat_torus_traits_2.h>

//for the time being this will run for conic polycurves only
#elif TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>


#elif TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

#elif TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>

#else
#error No geometry traits (GEOM_TRAITS) specified!
#endif

// ============================================================================
// Geometry Traits typedef:
// ============================================================================
#if TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Base_geom_traits;
#define GEOM_TRAITS_TYPE "Segments"

#elif TEST_GEOM_TRAITS == NON_CACHING_SEGMENT_GEOM_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Base_geom_traits;
#define GEOM_TRAITS_TYPE "Non Caching Segments"

#elif TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Base_geom_traits;
// Poly curves needs some testing where Segments and X-monotone segments are
// required instead of polycurves/x-monotone polycurves.
typedef Base_geom_traits::Subcurve_2                    Subcurve_2;
typedef Base_geom_traits::X_monotone_subcurve_2         X_monotone_subcurve_2;
#define GEOM_TRAITS_TYPE "Polylines"

#elif TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Base_geom_traits;
// Poly curves needs some testing where Segments and X-monotone segments are
// required instead of polycurves/x-monotone polycurves.
typedef Base_geom_traits::Subcurve_2                    Subcurve_2;
typedef Base_geom_traits::X_monotone_subcurve_2         X_monotone_subcurve_2;
#define GEOM_TRAITS_TYPE "Non Caching Polylines"

#elif TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Base_geom_traits;
typedef Base_geom_traits::Rat_point_2                   Rat_point;
typedef Base_geom_traits::Rat_circle_2                  Rat_circle;
typedef Base_geom_traits::Rat_segment_2                 Rat_segment;
#define GEOM_TRAITS_TYPE "Conics"

#elif TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS
typedef CGAL::Arr_linear_traits_2<Kernel>               Base_geom_traits;
#define GEOM_TRAITS_TYPE "Linear Lines"

#elif TEST_GEOM_TRAITS == LINE_ARC_GEOM_TRAITS
typedef Kernel                                          Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>
                                                        Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                        Circular_kernel;
typedef Circular_kernel::Line_arc_2                     Line_arc_2;
typedef CGAL::Arr_line_arc_traits_2<Circular_kernel>    Base_geom_traits;

#define GEOM_TRAITS_TYPE "Line Arc"

#elif TEST_GEOM_TRAITS == CIRCULAR_ARC_GEOM_TRAITS
typedef Kernel                                          Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>
                                                        Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                        Circular_kernel;
typedef CGAL::Arr_circular_arc_traits_2<Circular_kernel>
                                                        Base_geom_traits;

#define GEOM_TRAITS_TYPE "Circular Arc"

#elif TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS

typedef Kernel                                          Linear_kernel;
typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>
                                                        Algebraic_kernel;
typedef CGAL::Circular_kernel_2<Linear_kernel,Algebraic_kernel>
                                                        Circular_kernel;
typedef Circular_kernel::Circular_arc_2                 Circular_arc_2;
typedef Circular_kernel::Line_arc_2                     Line_arc_2;
typedef CGAL::Arr_circular_line_arc_traits_2<Circular_kernel>
                                                        Base_geom_traits;

#define GEOM_TRAITS_TYPE "Circular Line Arc"

#elif TEST_GEOM_TRAITS == CIRCLE_SEGMENT_GEOM_TRAITS
typedef CGAL::Cartesian<Number_type>                    Rat_kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>       Base_geom_traits;
typedef Rat_kernel::FT                                  Rat_nt;
typedef Rat_kernel::Circle_2                            Circle_2;
typedef Rat_kernel::Line_2                              Line_2;
typedef Rat_kernel::Segment_2                           Segment_2;
typedef Rat_kernel::Point_2                             Rat_point_2;
#define GEOM_TRAITS_TYPE "Circle Segments"

#elif TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Base_geom_traits;
typedef Base_geom_traits::Bezier_cache                  Bezier_cache;

#define GEOM_TRAITS_TYPE "Bezier"

#elif TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>
                                                        Base_geom_traits;
#define GEOM_TRAITS_TYPE "Geodesic Arc on Sphere"

#elif TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS

// TODO: This is for old Ron's traits---remove.
//typedef CGAL::CORE_algebraic_number_traits               Nt_traits; //OS - new
//typedef Nt_traits::Rational                                Rational;
//typedef Nt_traits::Algebraic                               Algebraic;
//typedef CGAL::Arr_rational_arc_traits_2<Kernel,Nt_traits>  Base_geom_traits;
//typedef Base_geom_traits::Rat_vector                       Rat_vector;
//#define GEOM_TRAITS_TYPE "Rational Arc"

typedef CGAL::Arr_rational_function_traits_2<Kernel>    Base_geom_traits;
typedef Base_geom_traits::Rational                      Rational;
typedef Base_geom_traits::Algebraic_real_1              Algebraic_real_1;
typedef Base_geom_traits::Rat_vector                    Rat_vector;
#define GEOM_TRAITS_TYPE "Rational Arc"

#elif TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS
typedef CGAL::Arr_algebraic_segment_traits_2<Number_type>
                                                        Base_geom_traits;
#define GEOM_TRAITS_TYPE "Algebraic"

#elif TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Conic_traits_2;
typedef CGAL::Arr_polycurve_traits_2<Conic_traits_2>    Base_geom_traits;

//needed when we want to construct ellips, parabola and hyperbolas
typedef Conic_traits_2::Rat_point_2                     Rat_point;
typedef Conic_traits_2::Rat_circle_2                    Rat_circle;
typedef Conic_traits_2::Rat_segment_2                   Rat_segment;

// Poly curves needs some testing where Segments and X-monotone segments are
// required instead of polycurves/x-monotone polycurves.
typedef Base_geom_traits::Subcurve_2                    Subcurve_2;
typedef Base_geom_traits::X_monotone_subcurve_2         X_monotone_subcurve_2;
#define GEOM_TRAITS_TYPE "Polycurves_conics"

#elif TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS

typedef CGAL::Cartesian<Number_type>                    Rat_kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>       Circle_segment_traits_2;
typedef CGAL::Arr_polycurve_traits_2<Circle_segment_traits_2>
                                                        Base_geom_traits;

typedef Rat_kernel::FT                                  Rat_nt;
typedef Rat_kernel::Circle_2                            Circle_2;
typedef Rat_kernel::Line_2                              Line_2;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Circle_segment_traits_2::Point_2                Point_2;

// Poly curves needs some testing where Segments and X-monotone segments are
// required instead of polycurves/x-monotone polycurves.
typedef Base_geom_traits::Subcurve_2                    Subcurve_2;
typedef Base_geom_traits::X_monotone_subcurve_2         X_monotone_subcurve_2;
#define GEOM_TRAITS_TYPE "Polycurves_circular_arcs"

#elif TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Bezier_tratis;
typedef Bezier_tratis::Bezier_cache                     Bezier_cache;
typedef Rat_kernel::Point_2                             Control_point_2;
typedef Bezier_tratis::Point_2                          Point_2;

typedef CGAL::Arr_polycurve_traits_2<Bezier_tratis>     Base_geom_traits;
// Poly curves needs some testing where Segments and X-monotone segments are
// required instead of polycurves/x-monotone polycurves.
typedef Base_geom_traits::Subcurve_2                    Subcurve_2;
typedef Base_geom_traits::X_monotone_subcurve_2         X_monotone_subcurve_2;


#define GEOM_TRAITS_TYPE "polycurve_bezier"

#else
#error No geometry traits (GEOM_TRAITS) specified!
#endif

typedef Base_geom_traits::Point_2                       Base_point_2;
typedef Base_geom_traits::Curve_2                       Base_curve_2;
typedef Base_geom_traits::X_monotone_curve_2            Base_x_monotone_curve_2;

#endif
