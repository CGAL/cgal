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

#include <QGraphicsView>

namespace CGAL
{
namespace Qt
{

GraphicsViewCurveInputBase::GraphicsViewCurveInputBase(
  QObject* parent, QGraphicsScene* scene) :
    GraphicsViewInput(parent),
    QGraphicsSceneMixin(scene)
{
}

void GraphicsViewCurveInputBase::setSnappingEnabled(bool b)
{
  this->snappingEnabled = b;
}

void GraphicsViewCurveInputBase::setSnapToGridEnabled(bool b)
{
  this->snapToGridEnabled = b;
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

CurveInputMethod::CurveInputMethod(
  QGraphicsScene* scene, CurveType type_, int num_points_) :
    QGraphicsSceneMixin(scene),
    type{type_}, num_points{num_points_}
{
  if (num_points > 0) clickedPoints.reserve(num_points);
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

CurveType CurveInputMethod::curveType() const { return type; }

void CurveInputMethod::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  if (this->clickedPoints.size() > 0)
    this->updateVisualGuideMouseMoved(
      this->clickedPoints, this->snapPoint(event));
}

void CurveInputMethod::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  if (event->button() == ::Qt::LeftButton)
  {
    QPointF clickedPoint = this->snapPoint(event);
    this->clickedPoints.push_back(clickedPoint);
    if (this->clickedPoints.size() < this->num_points)
    {
      if (this->clickedPoints.size() == 1) this->beginInput_();
      this->pointsGraphicsItem.insert(clickedPoint);
      this->updateVisualGuideNewPoint(this->clickedPoints);
    }
    else
    {
      curveGenerator->generate(this->clickedPoints, this->type);
      this->resetInput_();
    }
  }
  else
  {
    if (this->num_points == -1)
      curveGenerator->generate(this->clickedPoints, this->type);
    this->resetInput_();
  }
}

void CurveInputMethod::beginInput_()
{
  this->beginInput();
  this->getScene()->addItem(&(this->pointsGraphicsItem));
  for (auto& item : items)
    this->getScene()->addItem(item);
}

void CurveInputMethod::resetInput_()
{
  this->resetInput();
  this->clickedPoints.clear();
  this->pointsGraphicsItem.clear();
  this->getScene()->removeItem(&(this->pointsGraphicsItem));
  for (auto& item : items)
    this->getScene()->removeItem(item);
}

void CurveInputMethod::resetInput() { }
void CurveInputMethod::beginInput() { }

void CurveInputMethod::updateVisualGuideNewPoint(const std::vector<QPointF>&) {
}

void CurveInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>&, const QPointF&)
{
}

QPointF CurveInputMethod::snapPoint(QGraphicsSceneMouseEvent* event)
{
  return event->scenePos();
}

void CurveInputMethod::appendGraphicsItem(QGraphicsItem* item)
{
  items.push_back(item);
  item->setZValue(100);
}

// SegmentInputMethod
SegmentInputMethod::SegmentInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Segment, 2)
{
  this->appendGraphicsItem(&(this->segmentGuide));
}

void SegmentInputMethod::beginInput()
{
  this->segmentGuide.setLine(0, 0, 0, 0);
  QPen pen = this->segmentGuide.pen();
  pen.setColor(this->color);
  this->segmentGuide.setPen(pen);
}

void SegmentInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  QLineF segment{clickedPoints[0], movePoint};
  this->segmentGuide.setLine(segment);
}

// PolylineInputMethod
PolylineInputMethod::PolylineInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Polyline, -1)
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
  this->polylineGuide.setPen(pen);
  this->lastLine.setPen(pen);
}

void PolylineInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  this->lastLine.setLine(QLineF{clickedPoints.back(), movePoint});
}

void PolylineInputMethod::updateVisualGuideNewPoint(const std::vector<QPointF>& points)
{
  if (points.size() == 1)
    this->painterPath.moveTo(points.back());
  else
    this->painterPath.lineTo(points.back());

  this->polylineGuide.setPath(this->painterPath);
}

// RayInputMethod
RayInputMethod::RayInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Ray, 2)
{
  this->appendGraphicsItem(&(this->rayGuide));
}

void RayInputMethod::beginInput()
{
  this->rayGuide.setLine(0, 0, 0, 0);
  QPen pen = this->rayGuide.pen();
  pen.setColor(this->color);
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
LineInputMethod::LineInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Line, 2)
{
  this->appendGraphicsItem(&(this->lineGuide));
}

void LineInputMethod::beginInput()
{
  this->lineGuide.setLine(0, 0, 0, 0);
  QPen pen = this->lineGuide.pen();
  pen.setColor(this->color);
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
CircleInputMethod::CircleInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Circle, 2)
{
  this->appendGraphicsItem(&(this->circleGuide));
}

void CircleInputMethod::beginInput()
{
  this->circleGuide.setRect(0, 0, 0, 0);
  QPen pen = this->circleGuide.pen();
  pen.setColor(this->color);
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
EllipseInputMethod::EllipseInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::Ellipse, 2)
{
  this->appendGraphicsItem(&(this->ellipseGuide));
}

void EllipseInputMethod::beginInput()
{
  this->ellipseGuide.setRect(0, 0, 0, 0);
  QPen pen = this->ellipseGuide.pen();
  pen.setColor(this->color);
  this->ellipseGuide.setPen(pen);
}

void EllipseInputMethod::updateVisualGuideMouseMoved(
  const std::vector<QPointF>& clickedPoints, const QPointF& movePoint)
{
  this->ellipseGuide.setRect(QRectF{clickedPoints[0], movePoint});
}

// ThreePointCircularInputMethod
ThreePointCircularInputMethod::ThreePointCircularInputMethod(
  QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::ThreePointCircularArc, 3)
{
}

// FivePointConicInputMethod
FivePointConicInputMethod::FivePointConicInputMethod(QGraphicsScene* scene) :
    CurveInputMethod(scene, CurveType::FivePointConicArc, 5)
{
}

template <typename ArrTraits>
template <size_t... Is>
GraphicsViewCurveInput<ArrTraits>::GraphicsViewCurveInput(
  QObject* parent, QGraphicsScene* scene, std::index_sequence<Is...>) :
    GraphicsViewCurveInputBase(parent, scene),
    inputMethods{std::make_tuple((Is, scene)...)}
{
  this->setDefaultInputMethod(
    std::integral_constant<
      bool, std::tuple_size<InputMethodTuple>::value != 0>{});
  QObject::connect(
    &curveGenerator, SIGNAL(generate(CGAL::Object)), this,
    SIGNAL(generate(CGAL::Object)));
}

template <typename ArrTraits>
GraphicsViewCurveInput<ArrTraits>::GraphicsViewCurveInput(
  QObject* parent, QGraphicsScene* scene) :
    GraphicsViewCurveInput(
      parent, scene,
      std::make_index_sequence<std::tuple_size<InputMethodTuple>::value>())
{
  for_each(
    inputMethods, [&](auto&& it) { it.setCurveGenerator(&curveGenerator); });
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setCurveType(CurveType type)
{
  this->inputMethod = nullptr;
  for_each(inputMethods, [&](auto&& it) {
    if (it.curveType() == type)
      this->inputMethod = static_cast<CurveInputMethod*>(&it);
  });
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
  const std::vector<QPointF>& clickedPoints, CurveType type)
{
  boost::optional<CGAL::Object> obj;
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
  }
  if (obj) Q_EMIT generate(*obj);
}

boost::optional<CGAL::Object>
CurveGeneratorBase::generateSegment(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generatePolyline(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateRay(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateLine(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateCircle(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateEllipse(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateThreePointCircularArc(const std::vector<QPointF>&)
{
  return {};
}
boost::optional<CGAL::Object>
CurveGeneratorBase::generateFivePointConicArc(const std::vector<QPointF>&)
{
  return {};
}

// Curve Generator Segment Traits
template <typename Kernel_>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>>::generateSegment(
  const std::vector<QPointF>& clickedPoints)
{
  Converter<Kernel> convert;
  auto p0 = convert(clickedPoints[0]);
  auto p1 = convert(clickedPoints[1]);
  Curve_2 res(p0, p1);
  return CGAL::make_object(res);
}

// Curve Generator Polyline Traits
template <typename SegmentTraits>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>>::generatePolyline(
  const std::vector<QPointF>& clickedPoints)
{
  ArrTraits poly_tr;
  auto construct_poly = poly_tr.construct_curve_2_object();
  std::vector<Point_2> points;
  for(auto& p : clickedPoints)
    points.emplace_back(p.x(), p.y());
  Curve_2 res = construct_poly(points.begin(), points.end());
  return CGAL::make_object(res);
}

// Curve Generator Linear Traits
template <typename Kernel_>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateSegment(
  const std::vector<QPointF>& points)
{
  Point_2 p0{points[0].x(), points[0].y()};
  Point_2 p1{points[1].x(), points[1].y()};
  Curve_2 res = Curve_2(Segment_2(p0, p1));
  return CGAL::make_object(res);
}

template <typename Kernel_>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateRay(
  const std::vector<QPointF>& points)
{
  Point_2 p0{points[0].x(), points[0].y()};
  Point_2 p1{points[1].x(), points[1].y()};
  Curve_2 res = Curve_2(Ray_2(p0, p1));
  return CGAL::make_object(res);
}

template <typename Kernel_>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateLine(
  const std::vector<QPointF>& points)
{
  Point_2 p0{points[0].x(), points[0].y()};
  Point_2 p1{points[1].x(), points[1].y()};
  Curve_2 res = Curve_2(Line_2(p0, p1));
  return CGAL::make_object(res);
}

// CurveGenerator Conic Traits
template <typename RatKernel, typename AlgKernel, typename NtTraits>
boost::optional<CGAL::Object>
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateSegment(const std::vector<QPointF>& points)
{
  Curve_2 res = Curve_2(Rat_segment_2(
    Rat_point_2(points[0].x(), points[0].y()),
    Rat_point_2(points[1].x(), points[1].y())));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
boost::optional<CGAL::Object>
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateCircle(const std::vector<QPointF>& points)
{
  float sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  Curve_2 res =
    Curve_2(Rat_circle_2(Rat_point_2(points[0].x(), points[0].y()), sq_rad));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
boost::optional<CGAL::Object>
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateEllipse(const std::vector<QPointF>& points)
{
  double x1 = (std::min)(points[0].x(), points[1].x());
  double y1 = (std::min)(points[0].y(), points[1].y());
  double x2 = (std::max)(points[0].x(), points[1].x());
  double y2 = (std::max)(points[0].y(), points[1].y());

  Rat_FT a = CORE::abs(Rat_FT(x1) - Rat_FT(x2)) / 2;
  Rat_FT b = CORE::abs(Rat_FT(y1) - Rat_FT(y2)) / 2;
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
boost::optional<CGAL::Object>
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateThreePointCircularArc(const std::vector<QPointF>& points)
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
boost::optional<CGAL::Object>
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateFivePointConicArc(const std::vector<QPointF>& points)
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

// CurveGenerator Circular Traits
template <typename CircularKernel>
boost::optional<CGAL::Object>
CurveGenerator<CGAL::Arr_circular_arc_traits_2<CircularKernel>>::
  generateThreePointCircularArc(const std::vector<QPointF>& points)
{
  auto& qp0 = points[0];
  auto& qp1 = points[1];
  auto& qp2 = points[2];
  Non_arc_point_2 pp0(qp0.x(), qp0.y());
  Non_arc_point_2 pp1(qp1.x(), qp1.y());
  Non_arc_point_2 pp2(qp2.x(), qp2.y());
  CircularKernel ker;
  if (!ker.collinear_2_object()(pp0, pp1, pp2))
  {
    Circle_2 circle(pp0, pp1, pp2);
    Circular_arc_2 arc(circle, pp0, pp2);
    std::vector<CGAL::Object> subarcs;
    CGAL::make_x_monotone(arc, std::back_inserter(subarcs));
    typename CircularKernel::Has_on_2 has_on;
    bool isOn = false;
    for (unsigned int i = 0; i < subarcs.size(); ++i)
    {
      Circular_arc_2 subarc;
      if (CGAL::assign(subarc, subarcs[i]))
      {
        if (has_on(subarc, pp1))
        {
          isOn = true;
          break;
        }
      }
    }

    if (isOn)
    {
      Curve_2 res(circle, pp0, pp2);
      return CGAL::make_object(res);
    }
    else
    {
      Curve_2 res(circle, pp2, pp0);
      return CGAL::make_object(res);
    }
  }
  else
  {
    std::cout << "Points don't specify a valid circular arc. Try again!\n";
  }
  return {};
}

template class GraphicsViewCurveInput<Seg_traits>;
template class GraphicsViewCurveInput<Pol_traits>;
template class GraphicsViewCurveInput<Conic_traits>;
template class GraphicsViewCurveInput<Lin_traits>;
template class GraphicsViewCurveInput<Arc_traits>;
template class GraphicsViewCurveInput<Alg_seg_traits>;

} // namespace Qt
} // namespace CGAL
