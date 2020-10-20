// Copyright (c) 2008, 2012, 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_TYPED_H
#define CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_TYPED_H

#include <CGAL/Object.h>
#include <boost/optional/optional_fwd.hpp>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Callback.h"
#include "CurveInputMethods.h"
#include "ForwardDeclarations.h"
#include "GraphicsSceneMixin.h"
#include "GraphicsViewCurveInput.h"
#include "PointSnapper.h"

class QEvent;

namespace demo_types
{
enum class TraitsType : int;
}

namespace CGAL
{
namespace Qt
{

template <typename ArrTraits_>
class CurveGeneratorBase
{
public:
  using ArrTraits = ArrTraits_;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

  void setTraits(const ArrTraits* traits_);

  boost::optional<Curve_2>
  generate(const std::vector<Point_2>& clickedPoints, CurveType type);

  virtual boost::optional<Curve_2>
  generateSegment(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateRay(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateLine(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generatePolyline(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateCircle(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateEllipse(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateThreePointCircularArc(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateFivePointConicArc(const std::vector<Point_2>&) { return {}; }
  virtual boost::optional<Curve_2>
  generateBezier(const std::vector<Point_2>&) { return {}; }

  const ArrTraits* traits;
};

template <typename ArrTraits_>
struct CurveGenerator : public CurveGeneratorBase<ArrTraits_>
{
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>> :
    public CurveGeneratorBase<CGAL::Arr_segment_traits_2<Kernel_>>
{
  using ArrTraits = CGAL::Arr_segment_traits_2<Kernel_>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Kernel = Kernel_;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;

  boost::optional<Curve_2>
  generateSegment(const std::vector<Point_2>&) override;
};

template <typename SegmentTraits>
struct CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>> :
    public CurveGeneratorBase<CGAL::Arr_polyline_traits_2<SegmentTraits>>
{
  using ArrTraits = CGAL::Arr_polyline_traits_2<SegmentTraits>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;

  boost::optional<Curve_2>
  generatePolyline(const std::vector<Point_2>&) override;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>> :
    public CurveGeneratorBase<
      Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>
{
  using ArrTraits = Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Kernel = AlgKernel;
  using Segment_2 = typename Kernel::Segment_2;
  using Rat_FT = typename RatKernel::FT;
  using Rat_point_2 = typename RatKernel::Point_2;
  using Rat_segment_2 = typename RatKernel::Segment_2;
  using Rat_circle_2 = typename RatKernel::Circle_2;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;

  boost::optional<Curve_2>
  generateSegment(const std::vector<Point_2>&) override;
  boost::optional<Curve_2> generateCircle(const std::vector<Point_2>&) override;
  boost::optional<Curve_2>
  generateEllipse(const std::vector<Point_2>&) override;
  boost::optional<Curve_2>
  generateThreePointCircularArc(const std::vector<Point_2>&) override;
  boost::optional<Curve_2>
  generateFivePointConicArc(const std::vector<Point_2>&) override;
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>> :
    public CurveGeneratorBase<CGAL::Arr_linear_traits_2<Kernel_>>
{
  using Kernel = Kernel_;
  using KernelPoint = typename Kernel::Point_2;
  using ArrTraits = CGAL::Arr_linear_traits_2<Kernel>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Ray_2 = typename Kernel::Ray_2;
  using Line_2 = typename Kernel::Line_2;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;

  boost::optional<Curve_2>
  generateSegment(const std::vector<Point_2>&) override;
  boost::optional<Curve_2> generateRay(const std::vector<Point_2>&) override;
  boost::optional<Curve_2> generateLine(const std::vector<Point_2>&) override;
};

template <typename Coefficient_>
struct CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>> :
    public CurveGeneratorBase<
      CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
  using Coefficient = Coefficient_;
  using ArrTraits = CGAL::Arr_algebraic_segment_traits_2<Coefficient>;
  using X_monotone_curve_2 = typename ArrTraits::X_monotone_curve_2;
  using Polynomial_2 = typename ArrTraits::Polynomial_2;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Algebraic_real_1 = typename ArrTraits::Algebraic_real_1;
  using Rational = typename Algebraic_real_1::Rational;
  using RationalTraits = Rational_traits<Rational>;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;

  boost::optional<Curve_2> generateLine(const std::vector<Point_2>&) override;
  boost::optional<Curve_2> generateCircle(const std::vector<Point_2>&) override;
  boost::optional<Curve_2>
  generateEllipse(const std::vector<Point_2>&) override;

private:
  boost::optional<Curve_2> generateEllipse_(const Point_2&, Rational, Rational);
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
struct CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>> :
    public CurveGeneratorBase<
      Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
  using ArrTraits =
    Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>;
  using Super = CurveGeneratorBase<ArrTraits>;
  using Point_2 = typename Super::Point_2;
  using Curve_2 = typename ArrTraits::Curve_2;

  boost::optional<Curve_2> generateBezier(const std::vector<Point_2>&) override;
};

template <typename ArrTraits>
struct GraphicsViewCurveInputTypeHelper
{
  using InputMethodTuple = std::tuple<>;
};

template <typename Kernel_>
struct GraphicsViewCurveInputTypeHelper<CGAL::Arr_segment_traits_2<Kernel_>>
{
  using InputMethodTuple = std::tuple<SegmentInputMethod>;
};

template <typename SegmentTraits>
struct GraphicsViewCurveInputTypeHelper<
  CGAL::Arr_polyline_traits_2<SegmentTraits>>
{
  using InputMethodTuple = std::tuple<PolylineInputMethod>;
};

template <typename Kernel_>
struct GraphicsViewCurveInputTypeHelper<CGAL::Arr_linear_traits_2<Kernel_>>
{
  using InputMethodTuple =
    std::tuple<SegmentInputMethod, RayInputMethod, LineInputMethod>;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct GraphicsViewCurveInputTypeHelper<
  CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>
{
  using InputMethodTuple = std::tuple<
    SegmentInputMethod, CircleInputMethod, EllipseInputMethod,
    ThreePointCircularInputMethod, FivePointConicInputMethod>;
};

template <typename Coefficient_>
struct GraphicsViewCurveInputTypeHelper<
  CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
  using InputMethodTuple =
    std::tuple<LineInputMethod, CircleInputMethod, EllipseInputMethod>;
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
struct GraphicsViewCurveInputTypeHelper<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
  using InputMethodTuple = std::tuple<BezierInputMethod>;
};

template <typename Arr_>
class GraphicsViewCurveInput :
    public GraphicsViewCurveInputBase,
    public CurveInputMethodCallback
{
  using Self = GraphicsViewCurveInput<Arr_>;
  using Arrangement = Arr_;
  using ArrTraits = typename Arrangement::Geometry_traits_2;
  using Curve_2 = typename Arrangement::Curve_2;
  using InputMethodTuple =
    typename GraphicsViewCurveInputTypeHelper<ArrTraits>::InputMethodTuple;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

public:
  GraphicsViewCurveInput(
    Arrangement* arrangement_, QObject* parent, QGraphicsScene* scene);
  void setCurveType(CurveType type) override;
  void setPointSnapper(PointSnapperBase*) override;
  void generate(CGAL::Object o) override;
  void curveInputDoneEvent(
    const std::vector<Point_2>& clickedPoints, CurveType type) override;

private:
  template <typename = ArrTraits>
  void setDefaultInputMethod(std::true_type);
  void setDefaultInputMethod(std::false_type);

  Arrangement* arrangement;
  InputMethodTuple inputMethods;
  CurveGenerator<ArrTraits> curveGenerator;
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
