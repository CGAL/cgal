// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#include "ConstructBoundingBox.h"
#include "ForwardDeclarations.h"
#include "ArrangementTypes.h"

#include <CGAL/number_utils.h>
#include <limits>

static const CGAL::Bbox_2 inf_bbox = {
  -std::numeric_limits<double>::infinity(),
  -std::numeric_limits<double>::infinity(),
  std::numeric_limits<double>::infinity(),
  std::numeric_limits<double>::infinity()};

static constexpr double inf_double = std::numeric_limits<double>::infinity();

template <typename Traits_>
struct ConstructBoundingBox_impl
{
  using Traits = Traits_;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
#ifdef CGAL_USE_CORE
    using Zero_resultant_exception = CGAL::internal::Zero_resultant_exception<
      typename demo_types::DemoTypes::Alg_seg_traits::Polynomial_2>;
#endif

    CGAL::Bbox_2 bbox;
    try
    {
      bbox = curve.bbox();
    }
    // algebraic traits sometimes crash when calling bbox
    // example: xy=0
    catch (std::exception& ex)
    {
      std::cerr << ex.what() << '\n';
      std::cerr << __FILE__ << ':' << __LINE__ << '\n';
      bbox = inf_bbox;
    }
#ifdef CGAL_USE_CORE
    catch (Zero_resultant_exception&)
    {
      std::cerr << "Exception thrown of type \"Zero_resultant_exception\"\n";
      std::cerr << __FILE__ << ':' << __LINE__ << '\n';
      bbox = inf_bbox;
    }
#endif
    catch (...)
    {
      std::cerr << "Exception thrown of unknown type!\n";
      std::cerr << __FILE__ << ':' << __LINE__ << '\n';
      bbox = inf_bbox;
    }
    return bbox;
  }

  CGAL::Bbox_2
  operator()(const Point_2& point)
  {
    double x = CGAL::to_double(point.x());
    double y = CGAL::to_double(point.y());
    return {x, y, x, y};
  }
};

// We currently avoid using bbox function in Arr_sgegment_2 because it creates
// its own kernel object
// TODO: remove this class and the polyline one once it's fixed
// and use bbox directly
template <typename Kernel_>
struct ConstructBoundingBox_impl<CGAL::Arr_segment_traits_2<Kernel_>>
{
  using Traits = CGAL::Arr_segment_traits_2<Kernel_>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
    return this->operator()(curve.source(), curve.target());
  }

  CGAL::Bbox_2 operator()(const Point_2& p1, const Point_2& p2)
  {
    CGAL::Bbox_2 bbox;
    double x1 = CGAL::to_double(p1.x());
    double y1 = CGAL::to_double(p1.y());
    double x2 = CGAL::to_double(p2.x());
    double y2 = CGAL::to_double(p2.y());

    double min_x, max_x, min_y, max_y;
    if (x1 < x2)
    {
      min_x = x1;
      max_x = x2;
    }
    else
    {
      min_x = x2;
      max_x = x1;
    }
    if (y1 < y2)
    {
      min_y = y1;
      max_y = y2;
    }
    else
    {
      min_y = y2;
      max_y = y1;
    }
    return {min_x, min_y, max_x, max_y};
  }

  CGAL::Bbox_2
  operator()(const Point_2& point)
  {
    double x = CGAL::to_double(point.x());
    double y = CGAL::to_double(point.y());
    return {x, y, x, y};
  }
};

template <typename SegmentTraits_2_>
struct ConstructBoundingBox_impl<CGAL::Arr_polyline_traits_2<SegmentTraits_2_>>
{
  using Traits = CGAL::Arr_polyline_traits_2<SegmentTraits_2_>;
  using SegmentTraits_2 = SegmentTraits_2_;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
    ConstructBoundingBox_impl<SegmentTraits_2> construct_bounding_box;

    auto n = curve.number_of_subcurves();
    CGAL::Bbox_2 bbox;
    for (std::size_t i = 0; i < n; ++i)
      bbox += construct_bounding_box(curve[i]);

    return bbox;
  }

  CGAL::Bbox_2 operator()(const Point_2& p1, const Point_2& p2)
  {
    CGAL::Bbox_2 bbox;
    double x1 = CGAL::to_double(p1.x());
    double y1 = CGAL::to_double(p1.y());
    double x2 = CGAL::to_double(p2.x());
    double y2 = CGAL::to_double(p2.y());

    double min_x, max_x, min_y, max_y;
    if (x1 < x2)
    {
      min_x = x1;
      max_x = x2;
    }
    else
    {
      min_x = x2;
      max_x = x1;
    }
    if (y1 < y2)
    {
      min_y = y1;
      max_y = y2;
    }
    else
    {
      min_y = y2;
      max_y = y1;
    }
    return {min_x, min_y, max_x, max_y};
  }

  CGAL::Bbox_2
  operator()(const Point_2& point)
  {
    double x = CGAL::to_double(point.x());
    double y = CGAL::to_double(point.y());
    return {x, y, x, y};
  }
};

template <typename Kernel_>
struct ConstructBoundingBox_impl<CGAL::Arr_linear_traits_2<Kernel_>>
{
  using Traits = CGAL::Arr_linear_traits_2<Kernel_>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
    if (curve.is_segment())
    {
      return ConstructBoundingBox_impl<CGAL::Arr_segment_traits_2<Kernel_>>{}(
        curve.source(), curve.target());
    }
    else if (curve.is_line())
    {
      return inf_bbox;
    }
    // ray
    else
    {
      auto&& ray = curve.ray();
      auto&& src = ray.source();
      double src_x = CGAL::to_double(src.x());
      double src_y = CGAL::to_double(src.y());
      auto&& dir = ray.direction();
      bool dx = CGAL::is_positive(dir.dx());
      bool dy = CGAL::is_positive(dir.dy());

      if (dx && dy)
        return {src_x, src_y, inf_double, inf_double};
      else if (!dx && dy)
        return {-inf_double, src_y, src_x, inf_double};
      else if (!dx && !dy)
        return {-inf_double, -inf_double, src_x, src_y};
      else // if (dx && !dy)
        return {src_x, -inf_double, inf_double, src_y};
    }
  }

  CGAL::Bbox_2 operator()(const Point_2& point)
  {
    double x = CGAL::to_double(point.x());
    double y = CGAL::to_double(point.y());
    return {x, y, x, y};
  }
};


template <
  typename RatKernel_, typename AlgKernel_, typename NtTraits_,
  typename BoundingTraits_>
struct ConstructBoundingBox_impl<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel_, AlgKernel_, NtTraits_, BoundingTraits_>>
{
  using Traits = typename CGAL::Arr_Bezier_curve_traits_2<
    RatKernel_, AlgKernel_, NtTraits_, BoundingTraits_>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
    // TODO: find a way to find the bounding box of a bezier X_monotone_curve
    return curve.supporting_curve().bbox();
  }

  CGAL::Bbox_2
  operator()(const Point_2& point)
  {
    std::pair<double, double> p = point.approximate();
    return {p.first, p.second, p.first, p.second};
  }
};

template <typename AlgebraicKernel_d_1>
struct ConstructBoundingBox_impl<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
  using Traits = CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Curve_2 = typename Traits::Curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2
  operator()(const X_monotone_curve_2& curve)
  {
    double min_x, max_x, min_y, max_y;
    if (curve.left_parameter_space_in_x() == CGAL::ARR_INTERIOR)
      min_x = CGAL::to_double(curve.left_x());
    else
      min_x = -std::numeric_limits<double>::infinity();

    if (curve.right_parameter_space_in_x() == CGAL::ARR_INTERIOR)
      max_x = CGAL::to_double(curve.right_x());
    else
      max_x = std::numeric_limits<double>::infinity();

    // TODO: calculate min/max points of rational functions to better
    // bound y coordinate
    // the problem is where to store them?
    // will also be useful for rendering

    min_y = -std::numeric_limits<double>::infinity();
    max_y = std::numeric_limits<double>::infinity();

    return {min_x, min_y, max_x, max_y};
  }

  CGAL::Bbox_2
  operator()(const Point_2& point)
  {
    double x = CGAL::to_double(point.x());
    double y = CGAL::to_double(point.y());
    return {x, y, x, y};
  }
};

template <typename Arr_>
CGAL::Bbox_2
ConstructBoundingBox<Arr_>::operator()(const X_monotone_curve_2& curve)
{
  return ConstructBoundingBox_impl<Traits>{}(curve);
}
template <typename Arr_>
CGAL::Bbox_2 ConstructBoundingBox<Arr_>::operator()(const Point_2& point)
{
  return ConstructBoundingBox_impl<Traits>{}(point);
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(ConstructBoundingBox)
