#include "MergeEdgeCallback.h"
#include "ArrangementTypes.h"
#include "CurveGraphicsItem.h"
#include <CGAL/Arrangement_with_history_2.h>
#include <QEvent>
#include <QGraphicsSceneMouseEvent>

#include "Utils.h"

/*! Constructor */
template <typename Arr_>
MergeEdgeCallback<Arr_>::MergeEdgeCallback(
  Arrangement* arr_, QObject* parent_) :
    CGAL::Qt::Callback(parent_),
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
    this->arr->merge_edge(this->mergeableHalfedge, nextHalfedge);
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
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  Kernel_point_2 p = CGAL::Qt::Converter<Kernel>{}(event->scenePos());
  double minDist = std::numeric_limits<double>::max();
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
      (!this->arr->are_mergeable(hei, h1)) &&
      (!this->arr->are_mergeable(hei, h2)))
    { continue; }

    X_monotone_curve_2 curve = hei->curve();
    Compute_squared_distance_2< Traits > squaredDistance;
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
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Kernel_point_2;

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
  else if (this->arr->are_mergeable(h, h1) && this->arr->are_mergeable(h, h2))
  {
    X_monotone_curve_2 c1 = h1->curve();
    X_monotone_curve_2 c2 = h2->curve();
    Compute_squared_distance_2< Traits > squaredDistance;
    squaredDistance.setScene(this->getScene());
    double d1 = CGAL::to_double(squaredDistance(p, c1));
    double d2 = CGAL::to_double(squaredDistance(p, c2));

    return (d1 < d2) ? h1 : h2;
  }
  else if (this->arr->are_mergeable(h, h2))
    return h2;
  else if (this->arr->are_mergeable(h, h1))
    return h1;
  else
    return Halfedge_handle();
}

template class MergeEdgeCallback<Seg_arr>;
template class MergeEdgeCallback<Pol_arr>;
template class MergeEdgeCallback<Conic_arr>;
template class MergeEdgeCallback<Lin_arr>;
template class MergeEdgeCallback<Alg_seg_arr>;
template class MergeEdgeCallback<Bezier_arr>;
