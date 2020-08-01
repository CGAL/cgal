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

#ifndef POINT_LOCATION_CALLBACK_H
#define POINT_LOCATION_CALLBACK_H

#include "Callback.h"

class QGraphicsScene;
class QGraphicsSceneMouseEvent;

namespace CGAL
{
class Arr_oblivious_side_tag;
class Arr_open_side_tag;

namespace Qt
{
template <typename T>
class CurveGraphicsItem;
}
} // namespace CGAL

/**
   Supports visualization of point location on arrangements.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < typename Arr_ >
class PointLocationCallback : public CGAL::Qt::Callback
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
  typedef typename Arrangement::Originating_curve_iterator
    Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Hole_const_iterator Hole_const_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  PointLocationCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );
  void setScene( QGraphicsScene* scene_ );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event );

  void highlightPointLocation( QGraphicsSceneMouseEvent *event,
                               const CGAL::Arr_oblivious_side_tag& );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event,
                               const CGAL::Arr_open_side_tag& );

  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
}; // class PointLocationCallback

#endif // POINT_LOCATION_CALLBACK_H
