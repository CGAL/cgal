#ifndef CGAL_TEST_TRAITS_H
#define CGAL_TEST_TRAITS_H

#include "test_geom_traits.h"

// Define Base_geom_traits to be the geometry traits, namely, Geom_traits.
typedef Base_geom_traits                Geom_traits;
typedef Geom_traits::Point_2            Point_2;
typedef Geom_traits::Curve_2            Curve_2;
typedef Geom_traits::X_monotone_curve_2 X_monotone_curve_2;

#include "test_topol_traits.h"

#endif
