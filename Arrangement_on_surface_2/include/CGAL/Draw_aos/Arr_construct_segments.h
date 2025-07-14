#ifndef CGAL_DRAW_AOS_ARR_CONSTRUCT_SEGMENTS_H
#define CGAL_DRAW_AOS_ARR_CONSTRUCT_SEGMENTS_H

#include "CGAL/CORE/BigRat.h"
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

template <typename GeomTraits, bool HasConstructXMonotoneCurve2>
class Arr_construct_segment_impl;

// Default implementation for traits that models Construct_x_monotone_curve_2
template <typename GeomTraits>
class Arr_construct_segment_impl<GeomTraits, true>
{
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<GeomTraits>::X_monotone_curve_2;
  using Construct_x_monotone_curve_2 = typename GeomTraits::Construct_x_monotone_curve_2;
  using FT = typename Traits_adaptor<GeomTraits>::FT;

public:
  Arr_construct_segment_impl(const GeomTraits& traits)
      : m_cst_x_curve(traits.construct_x_monotone_curve_2_object()) {}

  X_monotone_curve_2 operator()(FT x1, FT y1, FT x2, FT y2) const {
    return m_cst_x_curve(Point_2(x1, y1), Point_2(x2, y2));
  }

private:
  const Construct_x_monotone_curve_2 m_cst_x_curve;
};

// Specialization for Arr_circle_segment_traits_2
template <typename Kernel>
class Arr_construct_segment_impl<Arr_circle_segment_traits_2<Kernel>, false>
{
  using Geom_traits = Arr_circle_segment_traits_2<Kernel>;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using Line_2 = typename X_monotone_curve_2::Line_2;
  using Approximate_2 = typename Geom_traits::Approximate_2;
  using NT = typename Traits_adaptor<Geom_traits>::FT;
  using FT = typename Kernel::FT;

public:
  Arr_construct_segment_impl(const Geom_traits& traits) {}

  X_monotone_curve_2 operator()(NT x1, NT y1, NT x2, NT y2) const {
    using Kernel_point_2 = typename Kernel::Point_2;
    return X_monotone_curve_2(Kernel_point_2(CGAL::to_double(x1), CGAL::to_double(y1)),
                              Kernel_point_2(CGAL::to_double(x2), CGAL::to_double(y2)));
  }
};

template <typename GeomTraits>
class Arr_construct_segment_impl<GeomTraits, false>
{
  using Geom_traits = GeomTraits;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using FT = typename Traits_adaptor<Geom_traits>::FT;

public:
  Arr_construct_segment_impl(const GeomTraits&) { CGAL_assertion_msg(false, "Not implemented yet"); }

  X_monotone_curve_2 operator()(FT, FT, FT, FT) const {}
};

template <typename GeomTraits>
using Arr_construct_segment =
    Arr_construct_segment_impl<GeomTraits, has_construct_x_monotone_curve_2<GeomTraits>::value>;

template <typename GeomTraits>
class Arr_construct_vertical_segment
{
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<GeomTraits>::X_monotone_curve_2;
  using FT = typename Traits_adaptor<GeomTraits>::FT;

public:
  Arr_construct_vertical_segment(const GeomTraits& traits)
      : m_cst_seg(traits) {}

  X_monotone_curve_2 operator()(FT x, FT ymin, FT ymax) const { return m_cst_seg(x, ymin, x, ymax); }

private:
  const Arr_construct_segment<GeomTraits> m_cst_seg;
};

template <typename GeomTraits>
class Arr_construct_horizontal_segment
{
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<GeomTraits>::X_monotone_curve_2;
  using FT = typename Traits_adaptor<GeomTraits>::FT;

public:
  Arr_construct_horizontal_segment(const GeomTraits& traits)
      : m_cst_seg(traits) {}

  X_monotone_curve_2 operator()(FT y, FT xmin, FT xmax) const { return m_cst_seg(xmin, y, xmax, y); }

private:
  const Arr_construct_segment<GeomTraits> m_cst_seg;
};

// Arr_construct_vertical_segment Specialization for Arr_rational_function_traits_2
template <typename Kernel>
class Arr_construct_vertical_segment<Arr_rational_function_traits_2<Kernel>>
{
  using Geom_traits = Arr_rational_function_traits_2<Kernel>;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using Construct_x_monotone_curve_2 = typename Geom_traits::Construct_x_monotone_curve_2;
  using FT = typename Traits_adaptor<Geom_traits>::FT;
  using Bound = typename Geom_traits::Bound;
  using Polynomial_1 = typename Geom_traits::Polynomial_1;

public:
  Arr_construct_vertical_segment(const Geom_traits&) {}

  X_monotone_curve_2 operator()(FT x0, FT ymin, FT ymax) const {
    Geom_traits traits;
    auto cst_x_curve = traits.construct_x_monotone_curve_2_object();

    // We could only construct a near vertical segment when dealing with rational functions.
    Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1);
    Polynomial_1 x0_num(CORE::numerator(x0.lower()));
    Polynomial_1 x0_denum(CORE::denominator(x0.lower()));
    double k = 100;
    Bound xmin = ymin.lower() / k + x0.lower();
    Bound xmax = ymax.lower() / k + x0.lower();
    Polynomial_1 p_num = x0_num * k * x - k * x0_denum;
    Polynomial_1 p_denum = x0_denum;
    auto cv = cst_x_curve(p_num, p_denum, FT(xmin), FT(xmax));
    return cv;
  }
};

// Arr_construct_horizontal_segment Specialization for Arr_rational_function_traits_2
template <typename Kernel>
class Arr_construct_horizontal_segment<Arr_rational_function_traits_2<Kernel>>
{
  using Geom_traits = Arr_rational_function_traits_2<Kernel>;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using FT = typename Traits_adaptor<Geom_traits>::FT;
  using Polynomial_1 = typename Geom_traits::Polynomial_1;
  using Construct_x_monotone_curve_2 = typename Geom_traits::Construct_x_monotone_curve_2;

public:
  Arr_construct_horizontal_segment(const Geom_traits&) {}

  X_monotone_curve_2 operator()(FT y, FT xmin, FT xmax) const {
    // The traits object is stateful. Currently there's a problem with internal cache that an exception will be
    // triggered after constructing a number of segments.
    // So we create a new traits object instead of reusing an old one.
    Geom_traits traits;
    auto cst_x_curve = traits.construct_x_monotone_curve_2_object();
    // TODO: only works for CORE now, support other number types
    Polynomial_1 num(CORE::numerator(y.lower()));
    Polynomial_1 denum(CORE::denominator(y.lower()));
    return cst_x_curve(num, denum, xmin, xmax);
  }
};

} // namespace draw_aos
} // namespace CGAL

#endif