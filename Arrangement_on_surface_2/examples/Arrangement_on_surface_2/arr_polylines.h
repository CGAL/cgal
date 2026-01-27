#ifndef ARR_POLYLINES_H
#define ARR_POLYLINES_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;

using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
using Point = Traits::Point_2;
using Segment = Traits::Segment_2;
using My_polyline = Traits::Curve_2;
using Arrangement = CGAL::Arrangement_2<Traits>;

#endif
