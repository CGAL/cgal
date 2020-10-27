#ifndef IS_IN_X_RANGE_H
#define IS_IN_X_RANGE_H

// Check whether the given point is in the x-range of the curve associated
// with the given halfedge.
template <typename Arrangement>
bool is_in_x_range(typename Arrangement::Halfedge_const_handle he,
                   const typename Arrangement::Point_2& p,
                   const typename Arrangement::Traits_2& traits)
{
  auto cmp = traits.compare_x_2_object();

  // Compare p with the source vertex (which may lie at x = +/- oo).
  CGAL::Arr_parameter_space src_px = he->source()->parameter_space_in_x();
  CGAL::Comparison_result res_s =
    (src_px == CGAL::ARR_LEFT_BOUNDARY) ? CGAL::SMALLER :
    ((src_px == CGAL::ARR_RIGHT_BOUNDARY) ? CGAL::LARGER :
     cmp(he->source()->point(), p));
  if (res_s == CGAL::EQUAL) return true;

  // Compare p with the target vertex (which may lie at x = +/- oo).
  CGAL::Arr_parameter_space trg_px = he->target()->parameter_space_in_x();
  CGAL::Comparison_result  res_t =
    (trg_px == CGAL::ARR_LEFT_BOUNDARY) ? CGAL::SMALLER :
    ((trg_px == CGAL::ARR_RIGHT_BOUNDARY) ? CGAL::LARGER :
     cmp(he->target()->point(), p));

  // p lies in the x-range of the halfedge iff its source and target lie
  // at opposite x-positions.
  return (res_s != res_t);
}

#endif
