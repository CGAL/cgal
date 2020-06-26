// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>,
//                 Saurabh Singh <ssingh@cs.iitr.ac.in>

#ifndef CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Callback.h"
#include "GraphicsSceneMixin.h"
#include "ISnappable.h"
#include "PointsGraphicsItem.h"

class QEvent;

namespace CORE
{
class BigRat;
} // namespace CORE

namespace CGAL
{

template <typename T>
class Arr_segment_traits_2;
template <typename T>
class Arr_polyline_traits_2;
template <typename T, typename U, typename I>
class Arr_conic_traits_2;
template <typename T, typename U, typename I, typename S>
class Arr_Bezier_curve_traits_2;
template <typename T>
class Arr_linear_traits_2;
template <typename T>
class Arr_algebraic_segment_traits_2;
template <typename T>
class Rational_traits;

namespace Qt
{
enum class CurveType
{
  Segment,
  Polyline,
  Ray,
  Line,
  Circle,
  Ellipse,
  ThreePointCircularArc,
  FivePointConicArc,
  AlgebraicEquation,
  Bezier,
};

class CurveGeneratorBase : public QObject
{
  Q_OBJECT

public:
  void generate(const std::vector<QPointF>& clickedPoints, CurveType type);

  virtual boost::optional<CGAL::Object>
  generateSegment(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateRay(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateLine(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generatePolyline(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateCircle(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateEllipse(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateThreePointCircularArc(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateFivePointConicArc(const std::vector<QPointF>&);

  virtual boost::optional<CGAL::Object>
  generateBezier(const std::vector<QPointF>&);

Q_SIGNALS:
  void generate(CGAL::Object);
};

template <typename ArrTraits_>
struct CurveGenerator : public CurveGeneratorBase
{
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>> :
    public CurveGeneratorBase
{
  using ArrTraits = CGAL::Arr_segment_traits_2<Kernel_>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = typename ArrTraits::Point_2;
  using Kernel = Kernel_;

  boost::optional<CGAL::Object>
  generateSegment(const std::vector<QPointF>&) override;
};

template <typename SegmentTraits>
struct CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>> :
    public CurveGeneratorBase
{
  using ArrTraits = CGAL::Arr_polyline_traits_2<SegmentTraits>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = typename ArrTraits::Point_2;

  boost::optional<CGAL::Object>
  generatePolyline(const std::vector<QPointF>&) override;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
struct CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>> :
    public CurveGeneratorBase
{
  using ArrTraits = Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = typename ArrTraits::Point_2;
  using Kernel = AlgKernel;
  using Segment_2 = typename Kernel::Segment_2;
  using Rat_FT = typename RatKernel::FT;
  using Rat_point_2 = typename RatKernel::Point_2;
  using Rat_segment_2 = typename RatKernel::Segment_2;
  using Rat_circle_2 = typename RatKernel::Circle_2;

  boost::optional<CGAL::Object>
  generateSegment(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateCircle(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateEllipse(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateThreePointCircularArc(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateFivePointConicArc(const std::vector<QPointF>&) override;
};

template <typename Kernel_>
struct CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>> :
    public CurveGeneratorBase
{
  using Kernel = Kernel_;
  using ArrTraits = CGAL::Arr_linear_traits_2<Kernel>;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Point_2 = typename Kernel::Point_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Ray_2 = typename Kernel::Ray_2;
  using Line_2 = typename Kernel::Line_2;

  boost::optional<CGAL::Object>
  generateSegment(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateRay(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateLine(const std::vector<QPointF>&) override;
};

template <typename Coefficient_>
struct CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>> :
    public CurveGeneratorBase
{
  using Coefficient = Coefficient_;
  using ArrTraits = CGAL::Arr_algebraic_segment_traits_2<Coefficient>;
  using X_monotone_curve_2 = typename ArrTraits::X_monotone_curve_2;
  using Point_2 = typename ArrTraits::Point_2;
  using Polynomial_2 = typename ArrTraits::Polynomial_2;
  using Curve_2 = typename ArrTraits::Curve_2;
  using Rational = CORE::BigRat; // FIX: should probably query rational type
  using RationalTraits = Rational_traits<Rational>;

  boost::optional<CGAL::Object>
  generateLine(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateCircle(const std::vector<QPointF>&) override;

  boost::optional<CGAL::Object>
  generateEllipse(const std::vector<QPointF>&) override;

private:
  boost::optional<CGAL::Object> generateEllipse_(const QPointF&, float, float);
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
struct CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>> :
    public CurveGeneratorBase
{
  boost::optional<CGAL::Object>
  generateBezier(const std::vector<QPointF>&) override;
};

class CurveInputMethod : public QGraphicsSceneMixin
{
public:
  CurveInputMethod(CurveType, int numPoints_ = -1);
  virtual ~CurveInputMethod() { }

  void setCurveGenerator(CurveGeneratorBase*);
  CurveType curveType() const;

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  virtual void mousePressEvent(QGraphicsSceneMouseEvent* event);
  void beginInput_();
  void resetInput_();

  void setColor(QColor);

  // override this to snap to the points you like
  virtual QPointF snapPoint(QGraphicsSceneMouseEvent* event);

protected:
  virtual void beginInput();
  virtual void resetInput();
  virtual void updateVisualGuideNewPoint(const std::vector<QPointF>&);
  virtual void
  updateVisualGuideMouseMoved(const std::vector<QPointF>&, const QPointF&);
  void appendGraphicsItem(QGraphicsItem* item);

protected:
  QColor color;

private:
  const int numPoints;
  std::vector<QPointF> clickedPoints;
  PointsGraphicsItem pointsGraphicsItem;
  CurveGeneratorBase* curveGenerator;
  const CurveType type;
  std::vector<QGraphicsItem*> items;
  bool itemsAdded;
};

class SegmentInputMethod : public CurveInputMethod
{
public:
  SegmentInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem segmentGuide;
};

class RayInputMethod : public CurveInputMethod
{
public:
  RayInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem rayGuide;
};

class LineInputMethod : public CurveInputMethod
{
public:
  LineInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsLineItem lineGuide;
};

class PolylineInputMethod : public CurveInputMethod
{
public:
  PolylineInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

  void updateVisualGuideNewPoint(const std::vector<QPointF>&) override;

private:
  QGraphicsPathItem polylineGuide;
  QGraphicsLineItem lastLine;
  QPainterPath painterPath;
};

class CircleInputMethod : public CurveInputMethod
{
public:
  CircleInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsEllipseItem circleGuide;
};

class EllipseInputMethod : public CurveInputMethod
{
public:
  EllipseInputMethod();

  void beginInput() override;

  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsEllipseItem ellipseGuide;
};

class ThreePointCircularInputMethod : public CurveInputMethod
{
public:
  ThreePointCircularInputMethod();
};

class FivePointConicInputMethod : public CurveInputMethod
{
public:
  FivePointConicInputMethod();
};

class BezierInputMethod: public CurveInputMethod
{
public:
  BezierInputMethod();
  void beginInput() override;
  void updateVisualGuideNewPoint(const std::vector<QPointF>&) override;
  void updateVisualGuideMouseMoved(
    const std::vector<QPointF>& clickedPoints,
    const QPointF& movePoint) override;

private:
  QGraphicsPathItem bezierOldGuide;
  QGraphicsPathItem bezierGuide;
  QPainterPath painterOldPath;
  QPainterPath painterPath;
  std::vector<QPointF> controlPoints;
  std::vector<QPointF> cache;
};

class GraphicsViewCurveInputBase :
    public GraphicsViewInput,
    public ISnappable,
    public QGraphicsSceneMixin
{
public:
  GraphicsViewCurveInputBase(QObject* parent, QGraphicsScene* scene);

  void setSnappingEnabled(bool b);
  void setSnapToGridEnabled(bool b);
  void setColor(QColor c);
  void reset();
  virtual void setCurveType(CurveType type) = 0;

protected:
  GraphicsViewCurveInputBase(QObject* parent);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  virtual void mousePressEvent(QGraphicsSceneMouseEvent* event);
  virtual bool eventFilter(QObject* obj, QEvent* event);

  bool snappingEnabled;
  bool snapToGridEnabled;
  CurveInputMethod* inputMethod; // active input method
};                               // class GraphicsViewCurveInputBase

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

template <typename ArrTraits>
class GraphicsViewCurveInput : public GraphicsViewCurveInputBase
{
  using InputMethodTuple =
    typename GraphicsViewCurveInputTypeHelper<ArrTraits>::InputMethodTuple;

public:
  GraphicsViewCurveInput(QObject* parent, QGraphicsScene* scene);
  void setCurveType(CurveType type) override;

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
