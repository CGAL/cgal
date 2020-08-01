// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "GraphicsViewCurveInput.h"
#include "ArrangementTypes.h"
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/CORE/BigRat.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Qt/Converter.h>
#include <QEvent>
#include "QtMetaTypes.h"

#include <QGraphicsView>

// TODO(Ahmed Essam): move these somewhere else!
template <std::size_t I = 0, typename FuncT, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
for_each(std::tuple<Tp...>&, FuncT)
{
}

template <std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if <
  I<sizeof...(Tp), void>::type for_each(std::tuple<Tp...>& t, FuncT f)
{
  f(std::get<I>(t));
  for_each<I + 1, FuncT, Tp...>(t, f);
}

namespace CGAL
{
namespace Qt
{

GraphicsViewCurveInputBase::GraphicsViewCurveInputBase(
  QObject* parent, QGraphicsScene* scene) :
    GraphicsViewInput(parent),
    QGraphicsSceneMixin(scene), inputMethod(nullptr)
{
}

void GraphicsViewCurveInputBase::reset()
{
  if (this->inputMethod)
  {
    this->inputMethod->resetInput_();
    this->inputMethod = nullptr;
  }
}

void GraphicsViewCurveInputBase::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  if (this->inputMethod) this->inputMethod->mouseMoveEvent(event);
}

void GraphicsViewCurveInputBase::mousePressEvent(
  QGraphicsSceneMouseEvent* event)
{
  if (this->inputMethod) this->inputMethod->mousePressEvent(event);
}

bool GraphicsViewCurveInputBase::eventFilter(QObject* obj, QEvent* event)
{
  if (event->type() == QEvent::GraphicsSceneMouseMove)
  {
    QGraphicsSceneMouseEvent* mouseEvent =
      static_cast<QGraphicsSceneMouseEvent*>(event);
    this->mouseMoveEvent(mouseEvent);
    return true;
  }
  else if (event->type() == QEvent::GraphicsSceneMousePress)
  {
    QGraphicsSceneMouseEvent* mouseEvent =
      static_cast<QGraphicsSceneMouseEvent*>(event);
    this->mousePressEvent(mouseEvent);
    return true;
  }

  return QObject::eventFilter(obj, event);
}

void GraphicsViewCurveInputBase::setColor(QColor c)
{
  this->inputMethod->setColor(c);
}

CurveInputMethod::CurveInputMethod(CurveType type_, int numPoints_) :
    QGraphicsSceneMixin(), numPoints{numPoints_}, type{type_}, itemsAdded{false}
{
  if (numPoints > 0) clickedPoints.reserve(numPoints);
  this->pointsGraphicsItem.setZValue(100);
}

void CurveInputMethod::setColor(QColor c)
{
  this->color = c;
  this->pointsGraphicsItem.setColor(c);
}

void CurveInputMethod::setCurveGenerator(CurveGeneratorBase* generator)
{
  this->curveGenerator = generator;
}

void CurveInputMethod::setPointSnapper(PointSnapperBase* snapper_)
{
  this->snapper = snapper_;
}

CurveType CurveInputMethod::curveType() const { return type; }


void CurveInputMethod::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  if (this->clickedPoints.size() > 0)
    this->updateVisualGuideMouseMoved(
      this->clickedPoints, this->snapQPoint(event));
}

void CurveInputMethod::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  if (event->button() == ::Qt::LeftButton)
  {
    Point_2 big_point = this->snapPoint(event);
    QPointF point = {
      CGAL::to_double(big_point.x()), CGAL::to_double(big_point.y())};

    // only accept unique consecutive points
    if (!this->clickedPoints.empty() && this->clickedPoints.back() == point)
      return;
    this->clickedPoints.push_back(point);
    this->clickedBigPoints.push_back(big_point);

    if (this->clickedPoints.size() < static_cast<size_t>(this->numPoints))
    {
      if (this->clickedPoints.size() == 1) this->beginInput_();
      this->pointsGraphicsItem.insert(point);
      this->updateVisualGuideNewPoint(this->clickedPoints);
    }
    else
    {
      curveGenerator->generate(this->clickedBigPoints, this->type);
      this->resetInput_();
    }
  }
  else
  {
    if (this->numPoints == -1)
      curveGenerator->generate(this->clickedBigPoints, this->type);
    this->resetInput_();
  }
}

void CurveInputMethod::beginInput_()
{
  this->beginInput();
  this->itemsAdded = true;
  this->getScene()->addItem(&(this->pointsGraphicsItem));
  for (auto& item : items) this->getScene()->addItem(item);
}

void CurveInputMethod::resetInput_()
{
  this->resetInput();
  this->clickedPoints.clear();
  this->clickedBigPoints.clear();
  this->pointsGraphicsItem.clear();
  if (this->itemsAdded)
  {
    this->getScene()->removeItem(&(this->pointsGraphicsItem));
    for (auto& item : items) this->getScene()->removeItem(item);
    this->itemsAdded = false;
  }
}

void CurveInputMethod::resetInput() { }
void CurveInputMethod::beginInput() { }

void CurveInputMethod::updateVisualGuideNewPoint(const std::vector<QPointF>&) {
}

void CurveInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>&, const QPointF&)
{
}

auto CurveInputMethod::snapPoint(QGraphicsSceneMouseEvent* event) -> Point_2
{
  return this->snapper->snapPoint(event->scenePos());
}

QPointF CurveInputMethod::snapQPoint(QGraphicsSceneMouseEvent* event)
{
  auto&& pt = this->snapPoint(event);
  return QPointF{CGAL::to_double(pt.x()), CGAL::to_double(pt.y())};
}

void CurveInputMethod::appendGraphicsItem(QGraphicsItem* item)
{
  items.push_back(item);
  item->setZValue(100);
}

// SegmentInputMethod
SegmentInputMethod::SegmentInputMethod() :
    CurveInputMethod(CurveType::Segment, 2)
{
  this->appendGraphicsItem(&(this->segmentGuide));
}

void SegmentInputMethod::beginInput()
{
  this->segmentGuide.setLine(0, 0, 0, 0);
  QPen pen = this->segmentGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->segmentGuide.setPen(pen);
}

void SegmentInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  QLineF segment{clickedPoints[0], movePoint};
  this->segmentGuide.setLine(segment);
}

// PolylineInputMethod
PolylineInputMethod::PolylineInputMethod() :
    CurveInputMethod(CurveType::Polyline, -1)
{
  this->appendGraphicsItem(&(this->polylineGuide));
  this->appendGraphicsItem(&(this->lastLine));
}

void PolylineInputMethod::beginInput()
{
  this->painterPath.clear();
  this->polylineGuide.setPath(this->painterPath);
  this->lastLine.setLine(0, 0, 0, 0);
  QPen pen = this->polylineGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->polylineGuide.setPen(pen);
  this->lastLine.setPen(pen);
}

void PolylineInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  this->lastLine.setLine(QLineF{clickedPoints.back(), movePoint});
}

void PolylineInputMethod::updateVisualGuideNewPoint(
  const std::vector<QPointF>& points)
{
  if (points.size() == 1)
    this->painterPath.moveTo(points.back());
  else
    this->painterPath.lineTo(points.back());

  this->polylineGuide.setPath(this->painterPath);
}

// RayInputMethod
RayInputMethod::RayInputMethod() : CurveInputMethod(CurveType::Ray, 2)
{
  this->appendGraphicsItem(&(this->rayGuide));
}

void RayInputMethod::beginInput()
{
  this->rayGuide.setLine(0, 0, 0, 0);
  QPen pen = this->rayGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->rayGuide.setPen(pen);
}

void RayInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  float length = QLineF{clickedPoints[0], movePoint}.length();
  float dx = (movePoint.x() - clickedPoints[0].x()) / length;
  float dy = (movePoint.y() - clickedPoints[0].y()) / length;
  QRectF viewport = this->viewportRect();
  float r = std::sqrt(
    viewport.width() * viewport.width() +
    viewport.height() * viewport.height());
  QPointF endPoint = {
    clickedPoints[0].x() + dx * r, clickedPoints[0].y() + dy * r};
  QLineF segment{clickedPoints[0], endPoint};
  this->rayGuide.setLine(segment);
}

// LineInputMethod
LineInputMethod::LineInputMethod() : CurveInputMethod(CurveType::Line, 2)
{
  this->appendGraphicsItem(&(this->lineGuide));
}

void LineInputMethod::beginInput()
{
  this->lineGuide.setLine(0, 0, 0, 0);
  QPen pen = this->lineGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->lineGuide.setPen(pen);
}

void LineInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  float length = QLineF{clickedPoints[0], movePoint}.length();
  float dx = (clickedPoints[0].x() - movePoint.x()) / length;
  float dy = (clickedPoints[0].y() - movePoint.y()) / length;
  QRectF viewport = this->viewportRect();
  float r = std::sqrt(
    viewport.width() * viewport.width() +
    viewport.height() * viewport.height());
  QPointF endPoint = {
    clickedPoints[0].x() + dx * r, clickedPoints[0].y() + dy * r};
  QPointF firstPoint = {
    clickedPoints[0].x() - dx * r, clickedPoints[0].y() - dy * r};
  QLineF segment{firstPoint, endPoint};
  this->lineGuide.setLine(segment);
}

// CircleInputMethod
CircleInputMethod::CircleInputMethod() : CurveInputMethod(CurveType::Circle, 2)
{
  this->appendGraphicsItem(&(this->circleGuide));
}

void CircleInputMethod::beginInput()
{
  this->circleGuide.setRect(0, 0, 0, 0);
  QPen pen = this->circleGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->circleGuide.setPen(pen);
}

void CircleInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  auto& center = clickedPoints.front();
  float radius = QLineF{center, movePoint}.length();
  this->circleGuide.setRect(
    center.x() - radius, center.y() - radius, 2 * radius, 2 * radius);
}

// EllipseInputMethod
EllipseInputMethod::EllipseInputMethod() :
    CurveInputMethod(CurveType::Ellipse, 2)
{
  this->appendGraphicsItem(&(this->ellipseGuide));
}

void EllipseInputMethod::beginInput()
{
  this->ellipseGuide.setRect(0, 0, 0, 0);
  QPen pen = this->ellipseGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->ellipseGuide.setPen(pen);
}

void EllipseInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  this->ellipseGuide.setRect(QRectF{clickedPoints[0], movePoint});
}

// ThreePointCircularInputMethod
ThreePointCircularInputMethod::ThreePointCircularInputMethod() :
    CurveInputMethod(CurveType::ThreePointCircularArc, 3)
{
}

// FivePointConicInputMethod
FivePointConicInputMethod::FivePointConicInputMethod() :
    CurveInputMethod(CurveType::FivePointConicArc, 5)
{
}

// BezierInputMethod
BezierInputMethod::BezierInputMethod() : CurveInputMethod(CurveType::Bezier, -1)
{
  this->appendGraphicsItem(&(this->bezierGuide));
  this->appendGraphicsItem(&(this->bezierOldGuide));
}

void BezierInputMethod::beginInput()
{
  this->painterOldPath.clear();
  this->painterPath.clear();
  this->bezierGuide.setPath(this->painterPath);
  this->bezierOldGuide.setPath(this->painterOldPath);

  QPen pen = this->bezierOldGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->bezierOldGuide.setPen(pen);

  // TODO: set bezier guide color dynamically
  pen = this->bezierGuide.pen();
  pen.setColor(QColorConstants::DarkGray);
  pen.setCosmetic(true);
  this->bezierGuide.setPen(pen);
}

static QPointF evalBezier(
  const std::vector<QPointF>& control_points, std::vector<QPointF>& cache,
  float t)
{
  // iterative de Casteljau Algorithm
  cache.clear();

  for (size_t j = 0; j + 1 < control_points.size(); j++)
    cache.push_back(control_points[j] * (1.0f - t) + control_points[j + 1] * t);

  for (size_t i = 1; i + 1 < control_points.size(); i++)
    for (size_t j = 0; j + i + 1 < control_points.size(); j++)
      cache[j] = cache[j] * (1.0f - t) + cache[j + 1] * t;

  return cache[0];
}

static float approx_pixel_length(
  const std::vector<QPointF>& control_points, const QTransform& worldTransform)
{
  float l = 0;
  for (size_t i = 0; i + 1 < control_points.size(); i++)
  {
    QPointF p1 = worldTransform.map(control_points[i]);
    QPointF p2 = worldTransform.map(control_points[i + 1]);
    l += QLineF{p1, p2}.length();
  }
  return l;
}

static void updateBezierPainterPath(
  const std::vector<QPointF>& controlPoints, std::vector<QPointF>& cache,
  const QTransform& worldTransform, QPainterPath& painterPath)
{
  painterPath.clear();
  if (controlPoints.size() < 2) return;

  float pixel_len = approx_pixel_length(controlPoints, worldTransform);
  static constexpr int PIXEL_DIV = 4;
  int num_segs = std::max(1, static_cast<int>(pixel_len / PIXEL_DIV));

  painterPath.moveTo(controlPoints[0]);
  for (int i = 0; i < num_segs; i++)
    painterPath.lineTo(
      evalBezier(controlPoints, cache, static_cast<float>(i + 1) / num_segs));
}

void BezierInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  this->controlPoints.clear();
  std::copy(
    clickedPoints.begin(), clickedPoints.end(),
    std::back_inserter(this->controlPoints));
  this->controlPoints.push_back(movePoint);

  updateBezierPainterPath(
    this->controlPoints, this->cache, this->getView()->transform(),
    this->painterPath);

  this->bezierGuide.setPath(this->painterPath);
}

void BezierInputMethod::updateVisualGuideNewPoint(
  const std::vector<QPointF>& clickedPoints)
{
  updateBezierPainterPath(
    clickedPoints, this->cache, this->getView()->transform(),
    this->painterOldPath);
  this->bezierOldGuide.setPath(this->painterPath);
}

template <typename ArrTraits>
GraphicsViewCurveInput<ArrTraits>::GraphicsViewCurveInput(
  QObject* parent, QGraphicsScene* scene) :
    GraphicsViewCurveInputBase(parent, scene)
{
  this->setDefaultInputMethod(
    std::integral_constant<
      bool, std::tuple_size<InputMethodTuple>::value != 0>{});
  QObject::connect(
    &curveGenerator, SIGNAL(generate(CGAL::Object)), this,
    SIGNAL(generate(CGAL::Object)));
  for_each(inputMethods, [&](auto&& it) {
    it.setScene(scene);
    it.setCurveGenerator(&curveGenerator);
  });
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setCurveType(CurveType type)
{
  this->reset();
  for_each(inputMethods, [&](auto&& it) {
    if (it.curveType() == type)
      this->inputMethod = static_cast<CurveInputMethod*>(&it);
  });
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setPointSnapper(
  PointSnapperBase* snapper_)
{
  for_each(inputMethods, [&](auto&& it) { it.setPointSnapper(snapper_); });
}


template <typename ArrTraits>
template <typename>
void GraphicsViewCurveInput<ArrTraits>::setDefaultInputMethod(std::true_type)
{
  this->inputMethod =
    static_cast<CurveInputMethod*>(&std::get<0>(inputMethods));
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setDefaultInputMethod(std::false_type)
{
}

// CurveGeneratorBase
void CurveGeneratorBase::generate(
  const std::vector<Point_2>& clickedPoints, CurveType type)
{
  CGAL::Object obj;
  switch (type)
  {
  case CurveType::Segment:
    obj = generateSegment(clickedPoints);
    break;
  case CurveType::Ray:
    obj = generateRay(clickedPoints);
    break;
  case CurveType::Line:
    obj = generateLine(clickedPoints);
    break;
  case CurveType::Polyline:
    obj = generatePolyline(clickedPoints);
    break;
  case CurveType::Circle:
    obj = generateCircle(clickedPoints);
    break;
  case CurveType::Ellipse:
    obj = generateEllipse(clickedPoints);
    break;
  case CurveType::ThreePointCircularArc:
    obj = generateThreePointCircularArc(clickedPoints);
    break;
  case CurveType::FivePointConicArc:
    obj = generateFivePointConicArc(clickedPoints);
    break;
  case CurveType::Bezier:
    obj = generateBezier(clickedPoints);
  default:
    break;
  }
  if (!obj.empty()) Q_EMIT generate(obj);
}

CGAL::Object CurveGeneratorBase::generateSegment(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generatePolyline(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateRay(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateLine(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateCircle(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateEllipse(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object
CurveGeneratorBase::generateThreePointCircularArc(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object
CurveGeneratorBase::generateFivePointConicArc(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateBezier(const std::vector<Point_2>&)
{
  return {};
}

// Curve Generator Segment Traits
template <typename Kernel_>
CGAL::Object
CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& clickedPoints)
{
  Curve_2 res{clickedPoints[0], clickedPoints[1]};
  return CGAL::make_object(res);
}

// Curve Generator Polyline Traits
template <typename SegmentTraits>
CGAL::Object
CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>>::generatePolyline(
  const std::vector<Point_2>& clickedPoints)
{
  if (clickedPoints.size() < 2) return {};

  ArrTraits poly_tr;
  auto construct_poly = poly_tr.construct_curve_2_object();
  Curve_2 res = construct_poly(clickedPoints.begin(), clickedPoints.end());
  return CGAL::make_object(res);
}

// Curve Generator Linear Traits
template <typename Kernel_>
CGAL::Object
CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Segment_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename Kernel_>
CGAL::Object CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateRay(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Ray_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename Kernel_>
CGAL::Object CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateLine(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Line_2(points[0], points[1]));
  return CGAL::make_object(res);
}

// CurveGenerator Conic Traits
template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateSegment(const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Rat_segment_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateCircle(const std::vector<Point_2>& points)
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  Curve_2 res =
    Curve_2(Rat_circle_2(points[0], sq_rad));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateEllipse(const std::vector<Point_2>& points)
{
  auto x1 = CGAL::min(points[0].x(), points[1].x());
  auto y1 = CGAL::min(points[0].y(), points[1].y());
  auto x2 = CGAL::max(points[0].x(), points[1].x());
  auto y2 = CGAL::max(points[0].y(), points[1].y());

  Rat_FT a = CGAL::abs(Rat_FT(x1) - Rat_FT(x2)) / 2;
  Rat_FT b = CGAL::abs(Rat_FT(y1) - Rat_FT(y2)) / 2;
  Rat_FT a_sq = a * a;
  Rat_FT b_sq = b * b;
  Rat_FT x0 = (x2 + x1) / 2;
  Rat_FT y0 = (y2 + y1) / 2;

  Rat_FT r = b_sq;
  Rat_FT s = a_sq;
  Rat_FT t = 0;
  Rat_FT u = -2 * x0 * b_sq;
  Rat_FT v = -2 * y0 * a_sq;
  Rat_FT ww = x0 * x0 * b_sq + y0 * y0 * a_sq - a_sq * b_sq;

  Curve_2 res = Curve_2(r, s, t, u, v, ww);
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateThreePointCircularArc(const std::vector<Point_2>& points)
{
  auto& qp1 = points[0];
  auto& qp2 = points[1];
  auto& qp3 = points[2];
  Rat_point_2 p1 = Rat_point_2(qp1.x(), qp1.y());
  Rat_point_2 p2 = Rat_point_2(qp2.x(), qp2.y());
  Rat_point_2 p3 = Rat_point_2(qp3.x(), qp3.y());
  RatKernel ker;
  if (!ker.collinear_2_object()(p1, p2, p3))
  {
    Curve_2 res(p1, p2, p3);
    return CGAL::make_object(res);
  }
  else
  {
    std::cout << "Points don't specify a valid conic." << std::endl;
    return {};
  }
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateFivePointConicArc(const std::vector<Point_2>& points)
{
  auto& qp0 = points[0];
  auto& qp1 = points[1];
  auto& qp2 = points[2];
  auto& qp3 = points[3];
  auto& qp4 = points[4];
  Rat_point_2 p0 = Rat_point_2(qp0.x(), qp0.y());
  Rat_point_2 p1 = Rat_point_2(qp1.x(), qp1.y());
  Rat_point_2 p2 = Rat_point_2(qp2.x(), qp2.y());
  Rat_point_2 p3 = Rat_point_2(qp3.x(), qp3.y());
  Rat_point_2 p4 = Rat_point_2(qp4.x(), qp4.y());
  try
  {
    Curve_2 res(p0, p1, p2, p3, p4);
    if (res.is_valid())
      return CGAL::make_object(res);
    else
      std::cout << "Points don't specify a valid conic. Try again!"
                << std::endl;
  }
  catch (...)
  {
    std::cout << "Points don't specify a valid conic. Try again!" << std::endl;
  }
  return {};
}

// CurveGenerator Algebraic Traits
template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateLine(const std::vector<Point_2>& points)
{
  RationalTraits ratTraits;
  ArrTraits arrTraits;

  Rational dx = points[1].x() - points[0].x();
  Rational dy = points[1].y() - points[0].y();
  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);

  Polynomial_2 poly;
  if (dx != 0)
  {
    Rational mRat = dy / dx;
    Rational cRat = points[0].y() - mRat * points[0].x();
    // y = (a/b) x + (e/f)
    auto a = ratTraits.numerator(mRat);
    auto b = ratTraits.denominator(mRat);
    auto e = ratTraits.numerator(cRat);
    auto f = ratTraits.denominator(cRat);

    poly = b * f * y - f * a * x - b * e;
  }
  // vertical line
  else
  {
    Rational xP = points[0].x();
    auto a = ratTraits.numerator(xP);
    auto b = ratTraits.denominator(xP);
    poly = b * x - a;
  }

  auto construct_curve = arrTraits.construct_curve_2_object();
  auto res = construct_curve(poly);
  return CGAL::make_object(res);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateCircle(const std::vector<Point_2>& points)
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  return this->generateEllipse_(points[0], sq_rad, sq_rad);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse(const std::vector<Point_2>& points)
{
  auto rx =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) / 4.;
  auto ry =
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y()) / 4.;
  Point_2 center = {
    (points[1].x() + points[0].x()) / 2, (points[1].y() + points[0].y()) / 2};
  return this->generateEllipse_(center, rx, ry);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse_(const Point_2& center, Rational rxRat, Rational ryRat)
{
  RationalTraits ratTraits;
  ArrTraits arrTraits;

  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
  Rational xCenter = center.x();
  Rational yCenter = center.y();

  // (g/h) (x - a/b)^2 + (e/f) (y - c/d)^2 = (e/f) * (g/h)
  auto a = ratTraits.numerator(xCenter);
  auto b = ratTraits.denominator(xCenter);
  auto c = ratTraits.numerator(yCenter);
  auto d = ratTraits.denominator(yCenter);
  auto e = ratTraits.numerator(rxRat);
  auto f = ratTraits.denominator(rxRat);
  auto g = ratTraits.numerator(ryRat);
  auto h = ratTraits.denominator(ryRat);

  Polynomial_2 poly = g * f * CGAL::ipower(b * d * x - d * a, 2) +
                      e * h * CGAL::ipower(b * d * y - b * c, 2) -
                      e * g * b * b * d * d;

  auto construct_curve = arrTraits.construct_curve_2_object();
  auto res = construct_curve(poly);
  return CGAL::make_object(res);
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
CGAL::Object CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
  generateBezier(const std::vector<Point_2>& clickedPoints)
{
  using Traits =
    Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>;

  if (clickedPoints.size() < 2) return {};

  return CGAL::make_object(
    typename Traits::Curve_2{clickedPoints.begin(), clickedPoints.end()});
}

template class GraphicsViewCurveInput<Seg_traits>;
template class GraphicsViewCurveInput<Pol_traits>;
template class GraphicsViewCurveInput<Conic_traits>;
template class GraphicsViewCurveInput<Lin_traits>;
template class GraphicsViewCurveInput<Alg_seg_traits>;
template class GraphicsViewCurveInput<Bezier_traits>;

} // namespace Qt
} // namespace CGAL
