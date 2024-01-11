#ifndef ARR_CIRCULAR_H
#define ARR_CIRCULAR_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;
using Traits = CGAL::Arr_circle_segment_traits_2<Kernel>;
using CoordNT = Traits::CoordNT;
using Point = Traits::Point_2;
using Curve = Traits::Curve_2;
using X_monotone_curve = Traits::X_monotone_curve_2;
using Rational_point = Traits::Rational_point_2;
using Segment = Traits::Rational_segment_2;
using Circle = Traits::Rational_circle_2;
using Arrangement = CGAL::Arrangement_2<Traits>;

#endif
