#include "PointLocationCallback.h"
#include "PointLocationFunctions.h"
#include "CurveGraphicsItem.h"
#include "ArrangementTypes.h"

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
  QGraphicsSceneMouseEvent* event, const CGAL::Arr_oblivious_side_tag&)
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
  QGraphicsSceneMouseEvent* event, const CGAL::Arr_open_side_tag&)
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

template class PointLocationCallback<Seg_arr>;
template class PointLocationCallback<Pol_arr>;
template class PointLocationCallback<Conic_arr>;
template class PointLocationCallback<Lin_arr>;
template class PointLocationCallback<Alg_seg_arr>;
template class PointLocationCallback<Bezier_arr>;
