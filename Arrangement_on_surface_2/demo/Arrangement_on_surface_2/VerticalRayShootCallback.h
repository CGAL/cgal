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

#include <QEvent>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include "CurveGraphicsItem.h"
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <vector>

#include "Callback.h"
#include "Utils.h"
#include "VerticalRayGraphicsItem.h"

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
  using Callback::scene;
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
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::FT                           FT;
  typedef typename CGAL::Arr_trapezoid_ric_point_location< Arrangement >
    TrapezoidPointLocationStrategy;
  typedef typename CGAL::Arr_simple_point_location< Arrangement >
    SimplePointLocationStrategy;
  typedef typename CGAL::Arr_walk_along_line_point_location< Arrangement >
    Walk_pl_strategy;
  typedef typename CGAL::Arr_landmarks_point_location< Arrangement >
    LandmarksPointLocationStrategy;

  VerticalRayShootCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );
  void setScene( QGraphicsScene* scene_ );

  void slotModelChanged( );

  virtual void setEdgeWidth( int width )
  {
    this->highlightedCurves->setEdgeWidth( width );
    this->rayGraphicsItem.setWidth( width );
  }

  virtual void setEdgeColor( const QColor& color )
  {
    this->highlightedCurves->setEdgeColor( color );
    this->rayGraphicsItem.setColor( color );
  }

  virtual const QColor& edgeColor( ) const
  {
    return this->highlightedCurves->edgeColor( );
  }

  virtual int edgeWidth( ) const
  {
    return this->highlightedCurves->edgeWidth( );
  }

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event );
  Face_const_handle getFace( const CGAL::Object& o );
  CGAL::Object rayShootUp( const Kernel_point_2& point );
  CGAL::Object rayShootUp( const Kernel_point_2& point, CGAL::Tag_true );
  CGAL::Object rayShootUp( const Kernel_point_2& point, CGAL::Tag_false );
  CGAL::Object rayShootDown( const Kernel_point_2& point );
  CGAL::Object rayShootDown( const Kernel_point_2& point, CGAL::Tag_true );
  CGAL::Object rayShootDown( const Kernel_point_2& point, CGAL::Tag_false );

  using Superclass::scene;
  using Superclass::shootingUp;
  Traits traits;
  Arrangement* arr;
  Intersect_2 intersectCurves;
  CGAL::Qt::Converter< Kernel > convert;
  CGAL::Object pointLocationStrategy;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
  Kernel_point_2 queryPt;
  Arr_construct_point_2< Traits > toArrPoint;
  VerticalRayGraphicsItem rayGraphicsItem;
}; // class VerticalRayShootCallback

#endif // VERTICAL_RAY_SHOOT_CALLBACK_H
