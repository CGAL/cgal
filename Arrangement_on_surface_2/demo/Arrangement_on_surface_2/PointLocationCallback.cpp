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

#include "PointLocationCallback.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "CurveGraphicsItem.h"
#include "Utils/PointLocationFunctions.h"

#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

/**
   Supports visualization of point location on arrangements.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template <typename Arr_>
class PointLocationCallback : public PointLocationCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  typedef typename Arrangement::Face_handle Face_handle;
  typedef typename Arrangement::Face_const_handle Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle Vertex_const_handle;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Curve_handle Curve_handle;
  typedef
    typename Arrangement::Originating_curve_iterator Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Hole_const_iterator Hole_const_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  PointLocationCallback(Arrangement* arr_, QObject* parent_);
  void reset();
  void setScene(QGraphicsScene* scene_);

protected:
  void mousePressEvent(QGraphicsSceneMouseEvent* event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  void highlightPointLocation(QGraphicsSceneMouseEvent* event);

  void highlightPointLocation(
    QGraphicsSceneMouseEvent* event, CGAL::Arr_oblivious_side_tag);
  void highlightPointLocation(
    QGraphicsSceneMouseEvent* event, CGAL::Arr_open_side_tag);

  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem<Traits>* highlightedCurves;
}; // class PointLocationCallback

namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    res = new PointLocationCallback<Arrangement>(arr, parent);
  }

  PointLocationCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

PointLocationCallbackBase* PointLocationCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  PointLocationCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

/*! Constructor */
template <typename Arr_>
PointLocationCallback<Arr_>::PointLocationCallback(
  Arrangement* arr_, QObject* parent_) :
    PointLocationCallbackBase(parent_),
    arr(arr_), highlightedCurves(new CGAL::Qt::CurveGraphicsItem<Traits>())
{
  QObject::connect(
    this, SIGNAL(modelChanged()), this->highlightedCurves,
    SLOT(modelChanged()));
}

template <typename Arr_>
void PointLocationCallback<Arr_>::setScene(QGraphicsScene* scene_)
{
  this->scene = scene_;
  this->highlightedCurves->setScene(scene_);
  if (this->scene) { this->scene->addItem(this->highlightedCurves); }
}

template <typename Arr_>
void PointLocationCallback<Arr_>::reset()
{
  this->highlightedCurves->clear();
  Q_EMIT modelChanged();
}

template <typename Arr_>
void PointLocationCallback<Arr_>::mousePressEvent(
  QGraphicsSceneMouseEvent* event)
{
  this->highlightPointLocation(event);
}

template <typename Arr_>
void PointLocationCallback<Arr_>::mouseMoveEvent(
  QGraphicsSceneMouseEvent* /* event */)
{
}

template <typename Arr_>
void PointLocationCallback<Arr_>::highlightPointLocation(
  QGraphicsSceneMouseEvent* event)
{
  typename Traits::Left_side_category category;
  this->highlightPointLocation(event, category);

  Q_EMIT modelChanged();
}

template <typename Arr_>
void PointLocationCallback<Arr_>::highlightPointLocation(
  QGraphicsSceneMouseEvent* event, CGAL::Arr_oblivious_side_tag)
{
  Face_const_handle face =
    PointLocationFunctions<Arrangement>{}.getFace(this->arr, event->scenePos());

  this->highlightedCurves->clear();
  if (!face->is_unbounded())
  { // it is an interior face; highlight its border
    Ccb_halfedge_const_circulator cc = face->outer_ccb();
    do
    {
      X_monotone_curve_2 curve = cc->curve();
      this->highlightedCurves->insert(curve);
    } while (++cc != face->outer_ccb());
  }
  Hole_const_iterator hit;
  Hole_const_iterator eit = face->holes_end();
  for (hit = face->holes_begin(); hit != eit; ++hit)
  { // highlight any holes inside this face
    Ccb_halfedge_const_circulator cc = *hit;
    do
    {
      X_monotone_curve_2 curve = cc->curve();
      this->highlightedCurves->insert(curve);
      cc++;
    } while (cc != *hit);
  }
}

template <typename Arr_>
void PointLocationCallback<Arr_>::highlightPointLocation(
  QGraphicsSceneMouseEvent* event, CGAL::Arr_open_side_tag)
{
  Face_const_handle face =
    PointLocationFunctions<Arrangement>{}.getFace(this->arr, event->scenePos());

  this->highlightedCurves->clear();
  Ccb_halfedge_const_circulator cc = face->outer_ccb();
  do
  {
    if (!cc->is_fictitious())
    {
      X_monotone_curve_2 curve = cc->curve();
      this->highlightedCurves->insert(curve);
    }
  } while (++cc != face->outer_ccb());
  Hole_const_iterator hit;
  Hole_const_iterator eit = face->holes_end();
  for (hit = face->holes_begin(); hit != eit; ++hit)
  { // highlight any holes inside this face
    Ccb_halfedge_const_circulator cc = *hit;
    do
    {
      X_monotone_curve_2 curve = cc->curve();
      this->highlightedCurves->insert(curve);
      cc++;
    } while (cc != *hit);
  }
}
