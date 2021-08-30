// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "VerticalRayShootCallback.h"
#include "VerticalRayGraphicsItem.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "CurveGraphicsItem.h"
#include "Utils/PointLocationFunctions.h"
#include "Utils/Utils.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

template <typename Arr_>
class VerticalRayShootCallback : public VerticalRayShootCallbackBase
{
public:
  typedef VerticalRayShootCallbackBase Superclass;
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle Vertex_const_handle;
  typedef typename Arrangement::Geometry_traits_2 Traits;

  VerticalRayShootCallback(Arrangement* arr_, QObject* parent_);
  void reset() override;
  void setScene(QGraphicsScene* scene_) override;
  void slotModelChanged() override;
  void setEdgeWidth(int width) override;
  void setEdgeColor(const QColor& color) override;
  const QColor& edgeColor() const override;
  int edgeWidth() const override;

protected:
  void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
  void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;
  void highlightPointLocation(QGraphicsSceneMouseEvent* event);

  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem<Traits>* highlightedCurves;
  VerticalRayGraphicsItem rayGraphicsItem;
}; // class VerticalRayShootCallback

VerticalRayShootCallbackBase::VerticalRayShootCallbackBase(QObject* parent_) :
    CGAL::Qt::Callback(parent_), shootingUp(true)
{
}

// msvc2015 doesn't play well with polymorphic lambdas
namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    res = new VerticalRayShootCallback<Arrangement>(arr, parent);
  }

  VerticalRayShootCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

VerticalRayShootCallbackBase* VerticalRayShootCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  VerticalRayShootCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

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
  this->highlightedCurves->clear();
  QPointF queryQPt = event->scenePos();

  CGAL::Object pointLocationResult;
  if (this->shootingUp)
    pointLocationResult =
      PointLocationFunctions<Arrangement>{}.rayShootUp(arr, queryQPt);
  else
    pointLocationResult =
      PointLocationFunctions<Arrangement>{}.rayShootDown(arr, queryQPt);

  if (pointLocationResult.is_empty()) { return; }

  QRectF viewportRect = this->viewportRect();
  qreal y2;
  if (this->shootingUp)
  { // +y in Qt is towards the bottom
    y2 = viewportRect.bottom();
  }
  else
  {
    y2 = viewportRect.top();
  }
  Face_const_handle unboundedFace;
  Halfedge_const_handle halfedge;
  Vertex_const_handle vertex;
  if (CGAL::assign(unboundedFace, pointLocationResult))
  {
    this->rayGraphicsItem.setSource(queryQPt);
    this->rayGraphicsItem.setTargetY(y2);
    this->rayGraphicsItem.setIsInfinite(true);
  }
  else if (CGAL::assign(halfedge, pointLocationResult))
  {
    this->highlightedCurves->insert(halfedge->curve());

    // draw a ray from the clicked point to the hit curve
    Arr_compute_y_at_x_2<Traits> compute_y_at_x_2{arr->traits()};
    compute_y_at_x_2.setScene(this->getScene());
    double yApprox = compute_y_at_x_2.approx(halfedge->curve(), queryQPt.x());
    this->rayGraphicsItem.setSource(queryQPt);
    this->rayGraphicsItem.setTargetY(yApprox);
    this->rayGraphicsItem.setIsInfinite(false);
  }
  else if (CGAL::assign(vertex, pointLocationResult))
  {
    if (!vertex->is_at_open_boundary())
    {
      auto pt = vertex->point();
      this->highlightedCurves->insert(pt);
    }
  }

  Q_EMIT modelChanged();
}
