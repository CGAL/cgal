// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>,
//            Saurabh Singh <ssingh@cs.iitr.ac.in>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H

#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Callback.h"
#include "CurveInputMethods.h"
#include "ForwardDeclarations.h"
#include "GraphicsSceneMixin.h"
#include "PointSnapper.h"

class QEvent;

namespace CGAL
{
namespace Qt
{

class CurveGeneratorBase : public QObject, public CurveInputMethodCallback
{
  Q_OBJECT

public:
  using Point_2 = PointSnapperBase::Point_2;

  void curveInputDoneEvent(
    const std::vector<Point_2>& clickedPoints, CurveType type) override;
  virtual CGAL::Object generateSegment(const std::vector<Point_2>&);
  virtual CGAL::Object generateRay(const std::vector<Point_2>&);
  virtual CGAL::Object generateLine(const std::vector<Point_2>&);
  virtual CGAL::Object generatePolyline(const std::vector<Point_2>&);
  virtual CGAL::Object generateCircle(const std::vector<Point_2>&);
  virtual CGAL::Object generateEllipse(const std::vector<Point_2>&);
  virtual CGAL::Object
  generateThreePointCircularArc(const std::vector<Point_2>&);
  virtual CGAL::Object generateFivePointConicArc(const std::vector<Point_2>&);
  virtual CGAL::Object generateBezier(const std::vector<Point_2>&);

Q_SIGNALS:
  void generate(CGAL::Object);
};

template <typename ArrTraits_>
struct CurveGeneratorTypedBase : public CurveGeneratorBase
{
  using ArrTraits = ArrTraits_;
  void setTraits(const ArrTraits* traits_);
  const ArrTraits* traits;
};

template <typename ArrTraits_>
struct CurveGenerator : public CurveGeneratorTypedBase<ArrTraits_>
{
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>> :
    public CurveGeneratorTypedBase<CGAL::Arr_segment_traits_2<Kernel_>>
{
  using ArrTraits = CGAL::Arr_segment_traits_2<Kernel_>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Kernel = Kernel_;
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generateSegment(const std::vector<Point_2>&) override;
};

template <typename SegmentTraits>
struct CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>> :
    public CurveGeneratorTypedBase<CGAL::Arr_polyline_traits_2<SegmentTraits>>
{
  using ArrTraits = CGAL::Arr_polyline_traits_2<SegmentTraits>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generatePolyline(const std::vector<Point_2>&) override;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>> :
    public CurveGeneratorTypedBase<
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
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generateSegment(const std::vector<Point_2>&) override;
  CGAL::Object generateCircle(const std::vector<Point_2>&) override;
  CGAL::Object generateEllipse(const std::vector<Point_2>&) override;
  CGAL::Object
  generateThreePointCircularArc(const std::vector<Point_2>&) override;
  CGAL::Object generateFivePointConicArc(const std::vector<Point_2>&) override;
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>> :
    public CurveGeneratorTypedBase<CGAL::Arr_linear_traits_2<Kernel_>>
{
  using Kernel = Kernel_;
  using KernelPoint = typename Kernel::Point_2;
  using ArrTraits = CGAL::Arr_linear_traits_2<Kernel>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Ray_2 = typename Kernel::Ray_2;
  using Line_2 = typename Kernel::Line_2;
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generateSegment(const std::vector<Point_2>&) override;
  CGAL::Object generateRay(const std::vector<Point_2>&) override;
  CGAL::Object generateLine(const std::vector<Point_2>&) override;
};

template <typename Coefficient_>
struct CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>> :
    public CurveGeneratorTypedBase<
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
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generateLine(const std::vector<Point_2>&) override;
  CGAL::Object generateCircle(const std::vector<Point_2>&) override;
  CGAL::Object generateEllipse(const std::vector<Point_2>&) override;

private:
  CGAL::Object generateEllipse_(const Point_2&, Rational, Rational);
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
struct CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>> :
    public CurveGeneratorTypedBase<
      Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
  using Point_2 = CurveGeneratorBase::Point_2;

  CGAL::Object generateBezier(const std::vector<Point_2>&) override;
};

class GraphicsViewCurveInputBase : public Callback
{
  Q_OBJECT

public:
  GraphicsViewCurveInputBase(QObject* parent, QGraphicsScene* scene);

  void setColor(QColor c);
  void reset() override;
  bool eventFilter(QObject* obj, QEvent* event) override;
  virtual void setCurveType(CurveType type) = 0;
  virtual void setPointSnapper(PointSnapperBase*) = 0;

Q_SIGNALS:
  void generate(CGAL::Object);

protected:
  void setInputMethod(CurveInputMethod*);

private:
  // active input method
  CurveInputMethod* inputMethod;
}; // class GraphicsViewCurveInputBase

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

template <typename ArrTraits_>
class GraphicsViewCurveInput : public GraphicsViewCurveInputBase
{
  using ArrTraits = ArrTraits_;
  using InputMethodTuple =
    typename GraphicsViewCurveInputTypeHelper<ArrTraits>::InputMethodTuple;

public:
  GraphicsViewCurveInput(
    const ArrTraits* traits, QObject* parent, QGraphicsScene* scene);
  void setCurveType(CurveType type) override;
  void setPointSnapper(PointSnapperBase*) override;

private:
  template <typename = ArrTraits>
  void setDefaultInputMethod(std::true_type);
  void setDefaultInputMethod(std::false_type);

  InputMethodTuple inputMethods;
  CurveGenerator<ArrTraits> curveGenerator;
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
