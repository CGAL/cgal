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
#include "CurveGraphicsItem.h"
#include "Utils.h"
#include "PointLocationFunctions.h"

#include <CGAL/Qt/Converter.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

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
    arr(arr_),
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
void VerticalRayShootCallback<Arr_>::setEdgeWidth( int width )
{
  this->highlightedCurves->setEdgeWidth( width );
  this->rayGraphicsItem.setWidth( width );
}

template <typename Arr_>
void VerticalRayShootCallback<Arr_>::setEdgeColor( const QColor& color )
{
  this->highlightedCurves->setEdgeColor( color );
  this->rayGraphicsItem.setColor( color );
}

template <typename Arr_>
const QColor& VerticalRayShootCallback<Arr_>::edgeColor( ) const
{
  return this->highlightedCurves->edgeColor( );
}

template <typename Arr_>
int VerticalRayShootCallback<Arr_>::edgeWidth( ) const
{
  return this->highlightedCurves->edgeWidth( );
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
  // minimizing #includes in in the header file
  // Utils.h and ArrangementTypes.h are disasters
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::FT                           FT;

  this->highlightedCurves->clear();
  QPointF queryQPt = event->scenePos();
  Kernel_point_2 queryPt = CGAL::Qt::Converter<Kernel>{}(queryQPt);

  CGAL::Object pointLocationResult;
  if (this->shootingUp)
    pointLocationResult =
      PointLocationFunctions<Arrangement>{}.rayShootUp(arr, queryQPt);
  else
    pointLocationResult =
      PointLocationFunctions<Arrangement>{}.rayShootDown(arr, queryQPt);

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
    Kernel_point_2 p2(FT(queryPt.x()), y2);
    Segment_2 lineSegment(queryPt, p2);
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
    CoordinateType x(queryPt.x());
    double yApprox =
      CGAL::to_double(compute_y_at_x_2.approx(halfedge->curve(), x));
    FT yInt(yApprox);
    Kernel_point_2 p2(queryPt.x(), yInt);
    Segment_2 seg(queryPt, p2);
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

template class VerticalRayShootCallback<Seg_arr>;
template class VerticalRayShootCallback<Pol_arr>;
template class VerticalRayShootCallback<Conic_arr>;
template class VerticalRayShootCallback<Lin_arr>;
template class VerticalRayShootCallback<Alg_seg_arr>;
template class VerticalRayShootCallback<Bezier_arr>;
