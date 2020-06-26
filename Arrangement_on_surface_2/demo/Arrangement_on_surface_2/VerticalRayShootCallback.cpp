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

#include "VerticalRayShootCallback.h"

VerticalRayShootCallbackBase::VerticalRayShootCallbackBase(QObject* parent_) :
    CGAL::Qt::Callback(parent_), shootingUp(true)
{
}

//! displays the direction of the arrow relative to the point being selected
/*!
  \param isShootingUp boolean value to determine the direction
*/
void VerticalRayShootCallbackBase::setShootingUp(bool isShootingUp)
{
  this->shootingUp = isShootingUp;
}

template <typename Arr_>
VerticalRayShootCallback<Arr_>::VerticalRayShootCallback(
  Arrangement* arr_, QObject* parent_) :
    VerticalRayShootCallbackBase(parent_),
    arr(arr_), intersectCurves(this->traits.intersect_2_object()),
    pointLocationStrategy(CGAL::make_object(new Walk_pl_strategy(*arr_))),
    highlightedCurves(new CGAL::Qt::CurveGraphicsItem<Traits>())
{
  this->rayGraphicsItem.setZValue(100);

  this->highlightedCurves->setEdgeColor(this->rayGraphicsItem.color());
  this->highlightedCurves->setVertexColor(this->rayGraphicsItem.color());
  this->highlightedCurves->setZValue(100);

  QObject::connect(
    this, SIGNAL(modelChanged()), this->highlightedCurves,
    SLOT(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this, SLOT(slotModelChanged()));
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::setScene(QGraphicsScene* scene_)
{
  CGAL::Qt::Callback::setScene(scene_);
  this->highlightedCurves->setScene(scene_);
  if (scene_)
  {
    this->scene->addItem(this->highlightedCurves);
    this->scene->addItem(&this->rayGraphicsItem);
  }
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::slotModelChanged()
{
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::reset()
{
  this->rayGraphicsItem.reset();
  this->highlightedCurves->clear();
  Q_EMIT modelChanged();
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::mousePressEvent(
  QGraphicsSceneMouseEvent* event)
{
  this->highlightPointLocation(event);
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::mouseMoveEvent(
  QGraphicsSceneMouseEvent* /* event */)
{
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::highlightPointLocation(
  QGraphicsSceneMouseEvent* event)
{
  this->highlightedCurves->clear();
  this->queryPt = this->convert(event->scenePos());
  CGAL::Object pointLocationResult;
  if (this->shootingUp)
  { pointLocationResult = this->rayShootUp(this->queryPt); }
  else
  {
    pointLocationResult = this->rayShootDown(this->queryPt);
  }
  if (pointLocationResult.is_empty()) { return; }

  QRectF viewportRect = this->viewportRect();
  FT y2;
  if (this->shootingUp)
  { // +y in Qt is towards the bottom
    y2 = FT(viewportRect.bottom());
  }
  else
  {
    y2 = FT(viewportRect.top());
  }
  Face_const_handle unboundedFace;
  Halfedge_const_handle halfedge;
  Vertex_const_handle vertex;
  if (CGAL::assign(unboundedFace, pointLocationResult))
  {
    Kernel_point_2 p2(FT(this->queryPt.x()), y2);
    Segment_2 lineSegment(this->queryPt, p2);
    this->rayGraphicsItem.setSource(event->scenePos());
    this->rayGraphicsItem.setTargetY(CGAL::to_double(y2));
    this->rayGraphicsItem.setIsInfinite(true);
  }
  else if (CGAL::assign(halfedge, pointLocationResult))
  {
    this->highlightedCurves->insert(halfedge->curve());

    // draw a ray from the clicked point to the hit curve
    Arr_compute_y_at_x_2<Traits> compute_y_at_x_2;
    compute_y_at_x_2.setScene(this->getScene());
    CoordinateType x(this->queryPt.x());
    double yApprox =
      CGAL::to_double(compute_y_at_x_2.approx(halfedge->curve(), x));
    FT yInt(yApprox);
    Kernel_point_2 p2(this->queryPt.x(), yInt);
    Segment_2 seg(this->queryPt, p2);
    this->rayGraphicsItem.setSource(event->scenePos());
    this->rayGraphicsItem.setTargetY(CGAL::to_double(yInt));
    this->rayGraphicsItem.setIsInfinite(false);
  }
  else if (CGAL::assign(vertex, pointLocationResult))
  {
    if (!vertex->is_at_open_boundary())
    {
      Point_2 pt = vertex->point();
      this->highlightedCurves->insert(pt);
    }
  }

  Q_EMIT modelChanged();
}

template <typename Arr_>
typename VerticalRayShootCallback<Arr_>::Face_const_handle
VerticalRayShootCallback<Arr_>::getFace(const CGAL::Object& obj)
{
  Face_const_handle f;
  if (CGAL::assign(f, obj)) return f;

  Halfedge_const_handle he;
  if (CGAL::assign(he, obj)) return (he->face());

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign(v, obj));
  CGAL::assign(v, obj);
  if (v->is_isolated()) return v->face();
  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges();
  return (eit->face());
}

template <typename Arr_>
CGAL::Object
VerticalRayShootCallback<Arr_>::rayShootUp(const Kernel_point_2& point)
{
  typename Supports_landmarks<Arrangement>::Tag supportsLandmarks;
  return this->rayShootUp(point, supportsLandmarks);
}

template <typename Arr_>
CGAL::Object VerticalRayShootCallback<Arr_>::rayShootUp(
  const Kernel_point_2& pt, CGAL::Tag_true)
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy* landmarksStrategy;

  Point_2 point = this->toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
  { pointLocationResult = walkStrategy->ray_shoot_up(point); }
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_up(point);
  }
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = simpleStrategy->ray_shoot_up(point);
  }
  else if (CGAL::assign(landmarksStrategy, this->pointLocationStrategy))
  {
    // pointLocationResult = landmarksStrategy->locate( point );
    std::cerr << "Warning: landmarks point location strategy doesn't support "
                 "ray shooting"
              << std::endl;
    return CGAL::Object();
  }
  else
  {
    std::cout << "Didn't find the right strategy\n";
  }

  return pointLocationResult;
}

template <typename Arr_>
CGAL::Object VerticalRayShootCallback<Arr_>::rayShootUp(
  const Kernel_point_2& pt, CGAL::Tag_false)
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;

  Point_2 point = this->toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
  { pointLocationResult = walkStrategy->ray_shoot_up(point); }
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_up(point);
  }
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = simpleStrategy->ray_shoot_up(point);
  }
  else
  {
    std::cout << "Didn't find the right strategy\n";
  }

  return pointLocationResult;
}

template <typename Arr_>
CGAL::Object
VerticalRayShootCallback<Arr_>::rayShootDown(const Kernel_point_2& point)
{
  typename Supports_landmarks<Arrangement>::Tag supportsLandmarks;
  return this->rayShootDown(point, supportsLandmarks);
}

template <typename Arr_>
CGAL::Object VerticalRayShootCallback<Arr_>::rayShootDown(
  const Kernel_point_2& pt, CGAL::Tag_true)
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy* landmarksStrategy;

  Point_2 point = this->toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
  { pointLocationResult = walkStrategy->ray_shoot_down(point); }
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_down(point);
  }
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = simpleStrategy->ray_shoot_down(point);
  }
  else if (CGAL::assign(landmarksStrategy, this->pointLocationStrategy))
  {
    // pointLocationResult = landmarksStrategy->locate( point );
    std::cerr << "Warning: landmarks point location strategy doesn't support "
                 "ray shooting"
              << std::endl;
    return CGAL::Object();
  }
  return pointLocationResult;
}

template <typename Arr_>
CGAL::Object VerticalRayShootCallback<Arr_>::rayShootDown(
  const Kernel_point_2& pt, CGAL::Tag_false)
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;

  Point_2 point = this->toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
  { pointLocationResult = walkStrategy->ray_shoot_down(point); }
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_down(point);
  }
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
  {
    pointLocationResult = simpleStrategy->ray_shoot_down(point);
  }
  return pointLocationResult;
}

template class VerticalRayShootCallback<Seg_arr>;
template class VerticalRayShootCallback<Pol_arr>;
template class VerticalRayShootCallback<Conic_arr>;
template class VerticalRayShootCallback<Lin_arr>;
template class VerticalRayShootCallback<Alg_seg_arr>;
template class VerticalRayShootCallback<Bezier_arr>;
