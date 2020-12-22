#ifndef CGAL_TEST_CONFIGURATION_H
#define CGAL_TEST_CONFIGURATION_H

/*! This files contains define statements, include statement, and typedef
 * of the number types, kernel, and traits used.
 *
 */

//note that these values shloud match to the values in cgal_test script

#define DOUBLE_NT                          0
#define MP_FLOAT_NT                        1
#define GMPZ_NT                            2
#define LEDA_RAT_NT                        3
#define QUOTIENT_MP_FLOAT_NT               4
#define QUOTIENT_CGAL_GMPZ_NT              5
#define CGAL_GMPQ_NT                       6
#define LAZY_LEDA_RAT_NT                   7
#define LAZY_CGAL_GMPQ_NT                  8
#define LAZY_QUOTIENT_MP_FLOAT_NT          9
#define LEDA_REAL_NT                       10
#define CORE_EXPR_NT                       11
#define LAZY_GMPZ_NT                       12
#define LEDA_INT_NT                        13
#define CGAL_GMPZ_NT                       14
#define CORE_INT_NT                        15
#define CORE_RAT_NT                        16

#define CARTESIAN_KERNEL                   0
#define SIMPLE_CARTESIAN_KERNEL            1
#define UNIVARIATE_ALGEBRAIC_KERNEL        2

#define SEGMENT_GEOM_TRAITS                0
#define NON_CACHING_SEGMENT_GEOM_TRAITS    1
#define POLYLINE_GEOM_TRAITS               2
#define NON_CACHING_POLYLINE_GEOM_TRAITS   3
#define LINEAR_GEOM_TRAITS                 4
#define CORE_CONIC_GEOM_TRAITS             5
#define LINE_ARC_GEOM_TRAITS               6
#define CIRCULAR_ARC_GEOM_TRAITS           7
#define CIRCULAR_LINE_ARC_GEOM_TRAITS      8
#define CIRCLE_SEGMENT_GEOM_TRAITS         9
#define BEZIER_GEOM_TRAITS                 10
#define GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS 11
#define RATIONAL_ARC_GEOM_TRAITS           12
#define ALGEBRAIC_GEOM_TRAITS              13
#define POLYCURVE_CONIC_GEOM_TRAITS                      14
#define POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS 15
#define POLYCURVE_BEZIER_GEOM_TRAITS       16
#define FLAT_TORUS_GEOM_TRAITS             17

#define PLANAR_BOUNDED_TOPOL_TRAITS        0
#define PLANAR_UNBOUNDED_TOPOL_TRAITS      1
#define SPHERICAL_TOPOL_TRAITS             2

// Default value based on dependencies:
#if (TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) || \
         (TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS)
  #if !defined(TEST_NT)
  #define TEST_NT CORE_EXPR_NT
  #endif
#endif

// Default values:
#if !defined(TEST_KERNEL)
#define TEST_KERNEL CARTESIAN_KERNEL
#endif

#if !defined(TEST_NT)
#define TEST_NT QUOTIENT_MP_FLOAT_NT
#endif

#if !defined(TEST_GEOM_TRAITS)
#define TEST_GEOM_TRAITS SEGMENT_GEOM_TRAITS
#endif

#if !defined(TEST_TOPOL_TRAITS)
#if (TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS)
#define TEST_TOPOL_TRAITS SPHERICAL_TOPOL_TRAITS
#elif (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS)
#define TEST_TOPOL_TRAITS PLANAR_UNBOUNDED_TOPOL_TRAITS
#else
#define TEST_TOPOL_TRAITS PLANAR_BOUNDED_TOPOL_TRAITS
#endif
#endif

// Illegal combinations:
#if (TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS) || (TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS)
  #if(TEST_NT != CORE_EXPR_NT)
                #error "Core conic traits implies Core Expr number type!"
  #endif
#endif

#if (TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS) && (TEST_NT != CORE_EXPR_NT)
#error "Core bezier traits implies Core Expr number type!"
#endif

#if (TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS) && \
  (TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS)
#error "Geodesic geometry traits implies spherical topology traits"
#endif

#if ((TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS) ||      \
     (TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS)) && \
  (TEST_TOPOL_TRAITS != PLANAR_UNBOUNDED_TOPOL_TRAITS)
#error "Linear and Rational geometry traits implies planar unbounded topology traits"
#endif

#endif
