#ifndef APPROX_Y_AT_X_H
#define APPROX_Y_AT_X_H
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Segment_traits_2 = CGAL::Arr_segment_traits_2<Exact_kernel>;
using Arr = CGAL::Arrangement_2<Segment_traits_2>;
double approx_y_at_x(const Arr::Traits_2& traits, const Arr::X_monotone_curve_2& curve, double x) {}
#endif // APPROX_Y_AT_X_H