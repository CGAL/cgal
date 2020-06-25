#include "PointLocationCallback.h"
#include "CurveGraphicsItem.h"

#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

/*! Constructor */
template <typename Arr_>
PointLocationCallback<Arr_>::PointLocationCallback(
  Arrangement* arr_, QObject* parent_) :
    CGAL::Qt::Callback(parent_),
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
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;

  CGAL::Qt::Converter<Kernel> convert;
  Kernel_point_2 point = convert(event->scenePos());

  CGAL::Object pointLocationResult = this->locate(point);
  Face_const_handle face = this->getFace(pointLocationResult);
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
  typedef typename ArrTraitsAdaptor<Traits>::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;

  CGAL::Qt::Converter<Kernel> convert;
  Kernel_point_2 point = convert(event->scenePos());

  CGAL::Object pointLocationResult = this->locate(point);
  Face_const_handle face = this->getFace(pointLocationResult);
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

template <typename Arr_>
typename PointLocationCallback<Arr_>::Face_const_handle
PointLocationCallback<Arr_>::getFace(const CGAL::Object& obj)
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
CGAL::Object PointLocationCallback<Arr_>::locate(const Kernel_point_2& point)
{
  typename Supports_landmarks<Arrangement>::Tag supportsLandmarks;
  return this->locate(point, supportsLandmarks);
}

template <typename Arr_>
template <typename>
CGAL::Object
PointLocationCallback<Arr_>::locate(const Kernel_point_2& pt, CGAL::Tag_true)
{
  typedef typename Supports_landmarks<Arrangement>::LandmarksType
    LandmarksPointLocationStrategy;

  Arr_construct_point_2<Traits> toArrPoint;
  Point_2 point = toArrPoint(pt);

  LandmarksPointLocationStrategy pointLocationStrategy{*arr};
  CGAL::Object pointLocationResult = pointLocationStrategy.locate(point);

  return pointLocationResult;
}

template <typename Arr_>
CGAL::Object
PointLocationCallback<Arr_>::locate(const Kernel_point_2& pt, CGAL::Tag_false)
{
  typedef typename CGAL::Arr_walk_along_line_point_location<Arrangement>
    Walk_pl_strategy;

  Arr_construct_point_2<Traits> toArrPoint;
  Point_2 point = toArrPoint(pt);

  Walk_pl_strategy pointLocationStrategy{*arr};
  CGAL::Object pointLocationResult = pointLocationStrategy.locate(point);

  return pointLocationResult;
}

template class PointLocationCallback<Seg_arr>;
template class PointLocationCallback<Pol_arr>;
template class PointLocationCallback<Conic_arr>;
template class PointLocationCallback<Lin_arr>;
template class PointLocationCallback<Alg_seg_arr>;
// template class PointLocationCallback<Bezier_arr>;
