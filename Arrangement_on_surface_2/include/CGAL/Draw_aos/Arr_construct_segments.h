#ifndef CGAL_DRAW_AOS_ARR_CONSTRUCT_SEGMENTS_H
#define CGAL_DRAW_AOS_ARR_CONSTRUCT_SEGMENTS_H
#include "CGAL/Draw_aos/helpers.h"

namespace CGAL {
class Arr_construct_vertical_segment
{
  using Point_2 = Geom_traits::Point_2;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using FT = Geom_traits::FT;

public:
  Arr_construct_vertical_segment(const Geom_traits& traits)
      : m_traits(traits) {}

  X_monotone_curve_2 operator()(const FT& x, const FT& ymin, const FT& ymax) const {
    auto cst_x_curve = m_traits.construct_x_monotone_curve_2_object();
    return cst_x_curve(Point_2(x, ymin), Point_2(x, ymax));
  }

private:
  const Geom_traits& m_traits;
};

class Arr_construct_horizontal_segment
{
  using Point_2 = Geom_traits::Point_2;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using FT = Geom_traits::FT;

public:
  Arr_construct_horizontal_segment(const Geom_traits& traits)
      : m_traits(traits) {}

  X_monotone_curve_2 operator()(FT y, FT xmin, FT xmax) const {
    auto cst_x_curve = m_traits.construct_x_monotone_curve_2_object();
    return cst_x_curve(Point_2(xmin, y), Point_2(xmax, y));
  }

private:
  const Geom_traits& m_traits;
};
} // namespace CGAL

#endif