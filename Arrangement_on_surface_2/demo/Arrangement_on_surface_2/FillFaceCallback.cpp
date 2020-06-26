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

#include "FillFaceCallback.h"
#include "ArrangementTypes.h"
#include "Utils.h"

#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_with_history_2.h>

FillFaceCallbackBase::FillFaceCallbackBase(QObject* parent) :
    CGAL::Qt::Callback(parent), fillColor(::Qt::black)
{
}

//! sets the color of the selected viewport
/*!
  \param c A QColor object
*/
void FillFaceCallbackBase::setColor(QColor c)
{
  this->fillColor = c;
  Q_EMIT modelChanged();
}

//! gets the color of the fill color option.
/*!
  \return the color selected by the user
*/
QColor FillFaceCallbackBase::getColor() const { return this->fillColor; }

/*! Constructor */
template <class Arr_>
FillFaceCallback<Arr_>::FillFaceCallback(Arrangement* arr_, QObject* parent_) :
    FillFaceCallbackBase(parent_),
    pointLocationStrategy(CGAL::make_object(new Walk_pl_strategy(*arr_))),
    arr(arr_)
{
}

template <class Arr_>
void FillFaceCallback<Arr_>::reset()
{
  Q_EMIT modelChanged();
}

template <class Arr_>
void FillFaceCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  this->fillFace(event);
  Q_EMIT modelChanged();
}

template <class Arr_>
void FillFaceCallback<Arr_>::mouseMoveEvent(
  QGraphicsSceneMouseEvent* /* event */)
{
}

//! A Template type
//! Coloring the closed selected space
/*!
    \param event A QGrpahicsSceneMouseEvent pointer to the class
*/
template <class Arr_>
void FillFaceCallback<Arr_>::fillFace(QGraphicsSceneMouseEvent* event)
{
  Kernel_point_2 point = this->convert(event->scenePos());
  CGAL::Object pointLocationResult = this->locate(point);
  Face_const_handle face = this->getFace(pointLocationResult);
  Face_handle f = this->arr->non_const_handle(face);

  if (f->color() == ::Qt::white && this->fillColor.isValid())
    f->set_color(this->fillColor);
  else
    f->set_color(::Qt::white);
}

//! A Template type
//! get the selected face
/*!
    \param obj A CGAL::Object reference of the face
*/
template <class Arr_>
typename FillFaceCallback<Arr_>::Face_const_handle
FillFaceCallback<Arr_>::getFace(const CGAL::Object& obj)
{
  Face_const_handle f;
  if (CGAL::assign(f, obj)) { return f; }

  Halfedge_const_handle he;
  if (CGAL::assign(he, obj)) { return (he->face()); }

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign(v, obj));
  CGAL::assign(v, obj);
  if (v->is_isolated()) { return v->face(); }

  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges();
  return (eit->face());
}

//! A Template type
//! locating the mouse position
/*!
    \param point A Kernel_point_2 object
*/
template <class Arr_>
CGAL::Object FillFaceCallback<Arr_>::locate(const Kernel_point_2& point)
{
  typename Supports_landmarks<Arrangement>::Tag supportsLandmarks;
  return this->locate(point, supportsLandmarks);
}

//! A Template type
//! locating points given a query
/*!
    \param pt A kernel_point_2 object
    \param Tag_true A CGAL object available in class
*/
template <class Arr_>
template <typename>
CGAL::Object
FillFaceCallback<Arr_>::locate(const Kernel_point_2& pt, CGAL::Tag_true)
{
  typedef typename Supports_landmarks<Arrangement>::LandmarksType
    LandmarksPointLocationStrategy;

  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy*
    trapezoidStrategy; /*!< searching of a trapezoid from the given faces */
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy*
    landmarksStrategy; /*!< finds the neares neighbor point */

  Arr_construct_point_2<Traits> toArrPoint;
  Point_2 point = toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
    pointLocationResult = walkStrategy->locate(point);
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
    pointLocationResult = trapezoidStrategy->locate(point);
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
    pointLocationResult = simpleStrategy->locate(point);
  else if (CGAL::assign(landmarksStrategy, this->pointLocationStrategy))
    pointLocationResult = landmarksStrategy->locate(point);

  return pointLocationResult;
}

//! A Template type
//! locating points given a query
/*!
    \param pt A kernel_point_2 object
    \param Tag_false A CGAL object not available
*/
template <class Arr_>
CGAL::Object
FillFaceCallback<Arr_>::locate(const Kernel_point_2& pt, CGAL::Tag_false)
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;

  Arr_construct_point_2<Traits> toArrPoint;
  Point_2 point = toArrPoint(pt);

  if (CGAL::assign(walkStrategy, this->pointLocationStrategy))
    pointLocationResult = walkStrategy->locate(point);
  else if (CGAL::assign(trapezoidStrategy, this->pointLocationStrategy))
    pointLocationResult = trapezoidStrategy->locate(point);
  else if (CGAL::assign(simpleStrategy, this->pointLocationStrategy))
    pointLocationResult = simpleStrategy->locate(point);
  return pointLocationResult;
}

template class FillFaceCallback<Seg_arr>;
template class FillFaceCallback<Pol_arr>;
template class FillFaceCallback<Conic_arr>;
template class FillFaceCallback<Lin_arr>;
template class FillFaceCallback<Alg_seg_arr>;
template class FillFaceCallback<Bezier_arr>;
