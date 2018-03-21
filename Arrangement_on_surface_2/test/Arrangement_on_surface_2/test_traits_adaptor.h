#ifndef CGAL_TEST_TRAITS_ADAPTOR_H
#define CGAL_TEST_TRAITS_ADAPTOR_H

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include "test_geom_traits.h"

typedef CGAL::Arr_traits_adaptor_2<Base_geom_traits>       Geom_traits;
typedef Geom_traits::Point_2                               Point_2;
typedef Geom_traits::Curve_2                               Curve_2;
typedef Geom_traits::X_monotone_curve_2                    X_monotone_curve_2;

#endif
