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

#include "MergeEdgeCallback.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "CurveGraphicsItem.h"
#include "Utils/SplitAndMerge.h"
#include "Utils/Utils.h"

#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>

template <typename Arr_>
class MergeEdgeCallback : public MergeEdgeCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  typedef typename Arrangement::Vertex_iterator Vertex_iterator;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  MergeEdgeCallback(Arrangement* arr_, QObject* parent_);
  void setScene(QGraphicsScene* scene_) override;
  void reset() override;

protected:
  void mousePressEvent(QGraphicsSceneMouseEvent* event) override;
  void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;
  Halfedge_handle getNearestMergeableCurve(QGraphicsSceneMouseEvent* event);
  Halfedge_handle
  getNearestMergeableCurve(Halfedge_handle h, QGraphicsSceneMouseEvent* event);

  CGAL::Qt::CurveGraphicsItem<Traits>* highlightedCurve;
  CGAL::Qt::CurveGraphicsItem<Traits>* highlightedCurve2;
  Arrangement* arr;
  Merge_edge<Arrangement> mergeEdge;
  Halfedge_handle mergeableHalfedge;
  bool isFirst;
}; // class MergeEdgeCallback

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
    res = new MergeEdgeCallback<Arrangement>(arr, parent);
  }

  MergeEdgeCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

MergeEdgeCallbackBase* MergeEdgeCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  MergeEdgeCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

/*! Constructor */
template <typename Arr_>
MergeEdgeCallback<Arr_>::MergeEdgeCallback(
  Arrangement* arr_, QObject* parent_) :
    MergeEdgeCallbackBase(parent_),
    highlightedCurve(new CGAL::Qt::CurveGraphicsItem<Traits>()),
    highlightedCurve2(new CGAL::Qt::CurveGraphicsItem<Traits>()), arr(arr_),
    isFirst(true)
{
  QObject::connect(
    this, SIGNAL(modelChanged()), this->highlightedCurve, SLOT(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this->highlightedCurve2,
    SLOT(modelChanged()));
}

template <typename Arr_>
void MergeEdgeCallback<Arr_>::setScene(QGraphicsScene* scene_)
{
  Callback::setScene(scene_);
  this->highlightedCurve->setScene(scene_);
  this->highlightedCurve2->setScene(scene_);

  if (scene_)
  {
    this->scene->addItem(this->highlightedCurve);
    this->scene->addItem(this->highlightedCurve2);
  }
}

template <typename Arr_>
void MergeEdgeCallback<Arr_>::reset()
{
  this->isFirst = true;
  this->highlightedCurve->clear();
  this->highlightedCurve2->clear();
  this->mergeableHalfedge = Halfedge_handle();
}

template <typename Arr_>
void MergeEdgeCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  if (this->isFirst)
  { // save the first edge if mergeable
    Halfedge_handle halfedge = this->getNearestMergeableCurve(event);
    if (halfedge == Halfedge_handle()) { return; }
    this->isFirst = false;
    this->mergeableHalfedge = halfedge;
  }
  else
  {
    Halfedge_handle nextHalfedge =
      this->getNearestMergeableCurve(this->mergeableHalfedge, event);
    this->mergeEdge.mergeEdge(this->arr, this->mergeableHalfedge, nextHalfedge);
    this->reset();
  }

  Q_EMIT modelChanged();
}

template <typename Arr_>
void MergeEdgeCallback<Arr_>::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  if (this->isFirst)
  {
    Halfedge_handle halfedge = this->getNearestMergeableCurve(event);
    if (halfedge == Halfedge_handle()) { return; }
    this->highlightedCurve->clear();
    this->highlightedCurve->insert(halfedge->curve());
    Q_EMIT modelChanged();
  }
  else
  {
    Halfedge_handle nextHalfedge =
      this->getNearestMergeableCurve(this->mergeableHalfedge, event);

    if (nextHalfedge != Halfedge_handle())
    {
      this->highlightedCurve2->clear();
      this->highlightedCurve2->insert(nextHalfedge->curve());
      Q_EMIT modelChanged();
    }
  }
}

template <typename Arr_>
typename MergeEdgeCallback<Arr_>::Halfedge_handle
MergeEdgeCallback<Arr_>::getNearestMergeableCurve(
  QGraphicsSceneMouseEvent* event)
{
  // find the nearest curve to the cursor that is adjacent to a curve that
  // can be merged with it
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;

  Kernel_point_2 p = CGAL::Qt::Converter<Kernel>{}(event->scenePos());
  double minDist = (std::numeric_limits<double>::max)();
  Halfedge_iterator nearestHei;
  bool found = false;

  for (Halfedge_iterator hei = this->arr->halfedges_begin();
       hei != this->arr->halfedges_end(); ++hei)
  {
    Vertex_iterator source = hei->source();
    Vertex_iterator target = hei->target();
    if (source->degree() != 2 && target->degree() != 2)
    { // then this halfedge has no mergeable neighbors
      continue;
    }
    Halfedge_handle h1 = hei->prev();
    Halfedge_handle h2 = hei->next();
    if (
      (!this->mergeEdge.areMergeable(this->arr, hei, h1)) &&
      (!this->mergeEdge.areMergeable(this->arr, hei, h2)))
    { continue; }

    X_monotone_curve_2 curve = hei->curve();
    Compute_squared_distance_2<Traits> squaredDistance;
    squaredDistance.setScene(this->getScene());
    double dist = CGAL::to_double(squaredDistance(p, curve));
    if (!found || dist < minDist)
    {
      found = true;
      minDist = dist;
      nearestHei = hei;
    }
  }

  if (!found)
  { // then we did not find a mergeable halfedge
    return Halfedge_handle();
  }
  return nearestHei;
}

template <typename Arr_>
typename MergeEdgeCallback<Arr_>::Halfedge_handle
MergeEdgeCallback<Arr_>::getNearestMergeableCurve(
  Halfedge_handle h, QGraphicsSceneMouseEvent* event)
{
  // find the nearest curve to the cursor that is adjacent to a curve that
  // can be merged with it
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;

  Kernel_point_2 p = CGAL::Qt::Converter<Kernel>{}(event->scenePos());
  Halfedge_handle h1 = h->prev();
  Halfedge_handle h2 = h->next();
  Vertex_iterator source = h->source();
  Vertex_iterator target = h->target();

  if (source->degree() != 2 && target->degree() != 2)
    return Halfedge_handle();
  else if (source->degree() != 2)
    return h2;
  else if (target->degree() != 2)
    return h1;
  else if (
    this->mergeEdge.areMergeable(arr, h, h1) &&
    this->mergeEdge.areMergeable(arr, h, h2))
  {
    X_monotone_curve_2 c1 = h1->curve();
    X_monotone_curve_2 c2 = h2->curve();
    Compute_squared_distance_2<Traits> squaredDistance;
    squaredDistance.setScene(this->getScene());
    double d1 = CGAL::to_double(squaredDistance(p, c1));
    double d2 = CGAL::to_double(squaredDistance(p, c2));

    return (d1 < d2) ? h1 : h2;
  }
  else if (this->mergeEdge.areMergeable(arr, h, h2))
    return h2;
  else if (this->mergeEdge.areMergeable(arr, h, h1))
    return h1;
  else
    return Halfedge_handle();
}
