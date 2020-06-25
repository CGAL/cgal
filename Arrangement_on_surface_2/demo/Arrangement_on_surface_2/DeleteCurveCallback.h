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

#ifndef DELETE_CURVE_CALLBACK_H
#define DELETE_CURVE_CALLBACK_H

#include "Callback.h"

class QGraphicsScene;
class QGraphicsMouseEvent;

namespace CGAL
{
namespace Qt
{
template <typename T>
class CurveGraphicsItem;
}
} // namespace CGAL

enum class DeleteMode
{
  DeleteOriginatingCuve,
  DeleteEdge,
};

/**
   Handles deletion of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
class DeleteCurveCallbackBase : public CGAL::Qt::Callback
{
public:
  using CGAL::Qt::Callback::Callback;

  void setDeleteMode(DeleteMode deleteMode_) { this->deleteMode = deleteMode_; }
  DeleteMode getDeleteMode() { return this->deleteMode; }

protected:
  DeleteMode deleteMode;
};

template < typename Arr_ >
class DeleteCurveCallback : public DeleteCurveCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Curve_handle            Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
  Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator   Induced_edge_iterator;

  DeleteCurveCallback( Arrangement* arr_, QObject* parent_ );
  void setScene(QGraphicsScene* scene_) override;
  void reset() override;

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightNearestCurve( QGraphicsSceneMouseEvent *event );

  CGAL::Qt::CurveGraphicsItem<Traits>* highlightedCurve;
  Arrangement* arr;
  Halfedge_handle removableHalfedge;
}; // class DeleteCurveCallback

#endif // DELETE_CURVE_CALLBACK_H
