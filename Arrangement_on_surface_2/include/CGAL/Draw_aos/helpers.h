#ifndef CGAL_DRAW_AOS_HELPERS_H
#define CGAL_DRAW_AOS_HELPERS_H
#include "CGAL/Arr_linear_traits_2.h"
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
namespace CGAL {

using Exact_kernel = CGAL::Exact_predicates_exact_constructions_kernel;
// using Geom_traits = Arr_linear_traits_2<Exact_kernel>;
using Geom_traits = CGAL::Arr_linear_traits_2<Exact_kernel>;
using Arrangement = Arrangement_2<Geom_traits>;

struct Inner_ccb_tag
{};
struct Outer_ccb_tag
{};

} // namespace CGAL

#endif // CGAL_DRAW_AOS_HELPERS_H