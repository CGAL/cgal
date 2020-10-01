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

#include "CurveInputMethods.h"
#include <QEvent>
#include <QKeyEvent>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QtGlobal>

namespace CGAL
{
namespace Qt
{

CurveInputMethod::CurveInputMethod(CurveType type_, int numPoints_) :
    Callback(nullptr), numPoints{numPoints_}, callback{nullptr}, type{type_},
    itemsAdded{false}
{
  if (numPoints > 0) clickedPoints.reserve(numPoints);
  this->pointsGraphicsItem.setZValue(100);
}

void CurveInputMethod::setColor(QColor c)
{
  this->color = c;
  this->pointsGraphicsItem.setColor(c);
}

QColor CurveInputMethod::getColor() const
{
  return this->color;
}

void CurveInputMethod::setCallback(CurveInputMethodCallback* callback_)
{
  this->callback = callback_;
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
      if (this->callback)
        callback->curveInputDoneEvent(this->clickedBigPoints, this->type);
      this->reset();
    }
  }
  else if (event->button() == ::Qt::RightButton)
  {
    if (this->numPoints == -1 && this->callback)
      callback->curveInputDoneEvent(this->clickedBigPoints, this->type);
    this->reset();
  }
}

void CurveInputMethod::keyPressEvent(QKeyEvent* event)
{
  if (event->key() == ::Qt::Key_Escape)
    this->reset();
}

void CurveInputMethod::beginInput_()
{
  this->beginInput();
  this->itemsAdded = true;
  this->getScene()->addItem(&(this->pointsGraphicsItem));
  for (auto& item : items) this->getScene()->addItem(item);
}

static inline void clearPainterPath(QPainterPath& ppath)
{
#if (QT_VERSION >= QT_VERSION_CHECK(5, 13, 0))
  ppath.clear();
#else
  ppath = {};
#endif
}

void CurveInputMethod::reset()
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
  clearPainterPath(this->painterPath);
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
  if (clickedPoints[0] == movePoint)
  {
     this->lineGuide.setLine({});
     return;
  }
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
  clearPainterPath(this->painterOldPath);
  clearPainterPath(this->painterPath);
  this->bezierGuide.setPath(this->painterPath);
  this->bezierOldGuide.setPath(this->painterOldPath);

  QPen pen = this->bezierOldGuide.pen();
  pen.setColor(this->color);
  pen.setCosmetic(true);
  this->bezierOldGuide.setPen(pen);

  // TODO: set bezier guide color dynamically
  pen = this->bezierGuide.pen();
  pen.setColor(::Qt::darkGray);
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
  clearPainterPath(painterPath);
  if (controlPoints.size() < 2) return;

  float pixel_len = approx_pixel_length(controlPoints, worldTransform);
  static constexpr int PIXEL_DIV = 4;
  int num_segs = (std::max)(1, static_cast<int>(pixel_len / PIXEL_DIV));

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

} // namespace Qt
} // namespace CGAL
