#ifndef ARR_CONICS_H
#define ARR_CONICS_H

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point;
typedef Rat_kernel::Segment_2                         Rat_segment;
typedef Rat_kernel::Circle_2                          Rat_circle;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                      Traits;
typedef Traits::Point_2                               Point;
typedef Traits::Curve_2                               Conic_arc;
typedef Traits::X_monotone_curve_2                    X_monotone_conic_arc;
typedef CGAL::Arrangement_2<Traits>                   Arrangement;

#endif
