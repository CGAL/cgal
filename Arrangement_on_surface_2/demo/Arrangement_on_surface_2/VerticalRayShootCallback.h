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

#ifndef VERTICAL_RAY_SHOOT_CALLBACK_H
#define VERTICAL_RAY_SHOOT_CALLBACK_H

#include "Callback.h"
#include "VerticalRayGraphicsItem.h"

namespace CGAL
{
namespace Qt
{
template <typename T>
class CurveGraphicsItem;
}
}
class QGraphicsSceneMouseEvent;
class QGraphicsScene;

class VerticalRayShootCallbackBase : public CGAL::Qt::Callback
{
public:
  void setShootingUp( bool isShootingUp );

  virtual void setEdgeWidth( int width ) = 0;
  virtual void setEdgeColor( const QColor& color ) = 0;
  virtual const QColor& edgeColor( ) const = 0;
  virtual int edgeWidth( ) const = 0;

protected:
  VerticalRayShootCallbackBase( QObject* parent_ );
  bool shootingUp;
}; // class VerticalRayShootCallbackBase

/*
 * Supports visualization of vertical ray shooting on arrangements.
 *
 * The template parameter is a CGAL::Arrangement_with_history_2 of some type.
 */
template < typename Arr_ >
class VerticalRayShootCallback : public VerticalRayShootCallbackBase
{
public:
  typedef VerticalRayShootCallbackBase Superclass;
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Curve_handle            Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
    Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator   Induced_edge_iterator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;

  VerticalRayShootCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );
  void setScene( QGraphicsScene* scene_ );
  void slotModelChanged( );
  void setEdgeWidth( int width ) override;
  void setEdgeColor( const QColor& color ) override;
  const QColor& edgeColor( ) const override;
  int edgeWidth( ) const override;

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event );

  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
  VerticalRayGraphicsItem rayGraphicsItem;
}; // class VerticalRayShootCallback

#endif // VERTICAL_RAY_SHOOT_CALLBACK_H
