#ifndef IS_IN_X_RANGE_H
#define IS_IN_X_RANGE_H

// Check whether the given point is in the x-range of the given curve that
// represents a great-circle arc.
template <typename GeometryTraits>
bool is_in_x_range(const typename GeometryTraits::X_monotone_curve_2& c,
                   const typename GeometryTraits::Point_2& p,
                   const GeometryTraits& traits) {
  CGAL_assertion(! traits.is_on_y_identification_2_object()(p));

  if (traits.is_on_y_identification_2_object()(c)) return false;

  auto cmp_x_boundary = traits.compare_x_on_boundary_2_object();
  auto psy = traits.parameter_space_in_y_2_object();
  auto psy_min = psy(c, CGAL::ARR_MIN_END);
  if (psy_min == CGAL::ARR_BOTTOM_BOUNDARY)
    return cmp_x_boundary(p, c, CGAL::ARR_MIN_END) == CGAL::EQUAL;
  auto psy_max = psy(c, CGAL::ARR_MAX_END);
  if (psy_max == CGAL::ARR_TOP_BOUNDARY)
    return cmp_x_boundary(p, c, CGAL::ARR_MAX_END) == CGAL::EQUAL;
  auto psx = traits.parameter_space_in_x_2_object();
  auto psx_min = psx(c, CGAL::ARR_MIN_END);
  auto psx_max = psx(c, CGAL::ARR_MAX_END);
  if ((psx_min == CGAL::ARR_LEFT_BOUNDARY) &&
      (psx_min == CGAL::ARR_RIGHT_BOUNDARY))
    return true;
  auto cmp_x = traits.compare_x_2_object();
  if (psx_min == CGAL::ARR_LEFT_BOUNDARY) {
    const auto& p_right = traits.construct_max_vertex_2_object()(c);
    auto res = cmp_x(p, p_right);
    return ((res == CGAL::SMALLER) || (res == CGAL::EQUAL));
  }
  if (psx_max == CGAL::ARR_RIGHT_BOUNDARY) {
    const auto& p_left = traits.construct_min_vertex_2_object()(c);
    auto res = cmp_x(p_left, p);
    return ((res == CGAL::SMALLER) || (res == CGAL::EQUAL));
  }
  const auto& p_left = traits.construct_min_vertex_2_object()(c);
  auto res = cmp_x(p_left, p);
  if (res == CGAL::LARGER) return false;
  const auto& p_right = traits.construct_max_vertex_2_object()(c);
  res = cmp_x(p, p_right);
  if (res == CGAL::LARGER) return false;
  return true;
}

#endif
