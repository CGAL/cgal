// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
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
  CGAL::Object rayShootDown( const Kernel_point_2& point );

  using Superclass::scene;
  using Superclass::shootingUp;
  Traits traits;
  Arrangement* arr;
  Intersect_2 intersectCurves;
  CGAL::Qt::Converter< Kernel > convert;
  CGAL::Object pointLocationStrategy;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
  QGraphicsLineItem* activeRay;
  Kernel_point_2 queryPt;
  Arr_construct_point_2< Traits > toArrPoint;
  VerticalRayGraphicsItem rayGraphicsItem;
}; // class VerticalRayShootCallback

template < typename Arr_ >
VerticalRayShootCallback< Arr_ >::
VerticalRayShootCallback( Arrangement* arr_, QObject* parent_ ):
  VerticalRayShootCallbackBase( parent_ ),
  arr( arr_ ),
  intersectCurves( this->traits.intersect_2_object( ) ),
  pointLocationStrategy( CGAL::make_object( new Walk_pl_strategy( *arr_ ) ) ),
  highlightedCurves( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
  activeRay( new QGraphicsLineItem )
{
  this->rayGraphicsItem.setZValue( 100 );

  this->highlightedCurves->setEdgeColor( this->rayGraphicsItem.color( ) );
  this->highlightedCurves->setVertexColor( this->rayGraphicsItem.color( ) );
  this->highlightedCurves->setZValue( 100 );

  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this->highlightedCurves, SLOT( modelChanged( ) ) );
  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this, SLOT( slotModelChanged( ) ) );
}

template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
  this->highlightedCurves->setScene( scene_ );
  if ( this->scene )
  {
    this->scene->addItem( this->highlightedCurves );
    this->scene->addItem( this->activeRay );
    this->scene->addItem( &this->rayGraphicsItem );
  }
}


template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::slotModelChanged( )
{
  this->activeRay->update( );
}

template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::reset( )
{
  this->activeRay->setLine( 0, 0, 0, 0 );
  this->rayGraphicsItem.reset( );
  this->highlightedCurves->clear( );
  emit modelChanged( );
}

template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
  this->highlightPointLocation( event );
}

template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::
mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */)
{ }

template < typename Arr_ >
void VerticalRayShootCallback< Arr_ >::
highlightPointLocation( QGraphicsSceneMouseEvent* event )
{
  this->highlightedCurves->clear( );
  this->queryPt = this->convert( event->scenePos( ) );
  CGAL::Object pointLocationResult;
  if ( this->shootingUp )
  {
    pointLocationResult = this->rayShootUp( this->queryPt );
  }
  else
  {
    pointLocationResult = this->rayShootDown( this->queryPt );
  }
  if ( pointLocationResult.is_empty( ) )
  {
    return;
  }
    
  QRectF viewportRect = this->viewportRect( );
  FT y2;
  if ( this->shootingUp )
  { // +y in Qt is towards the bottom
    y2 = FT( viewportRect.bottom( ) );
  }
  else
  {
    y2 = FT( viewportRect.top( ) );
  }
  Face_const_handle unboundedFace;
  Halfedge_const_handle halfedge;
  Vertex_const_handle vertex;
  if ( CGAL::assign( unboundedFace, pointLocationResult ) )
  {
    Kernel_point_2 p2( FT( this->queryPt.x( ) ), y2 );
    Segment_2 lineSegment( this->queryPt, p2 );
    // QLineF qLineSegment = this->convert( lineSegment );
    // this->activeRay->setLine( qLineSegment );
    this->rayGraphicsItem.setSource( event->scenePos( ) );
    this->rayGraphicsItem.setTargetY( CGAL::to_double( y2 ) );
    this->rayGraphicsItem.setIsInfinite( true );
  }
  else if ( CGAL::assign( halfedge, pointLocationResult ) )
  {
    this->highlightedCurves->insert( halfedge->curve( ) );

    // draw a ray from the clicked point to the hit curve
    Arr_compute_y_at_x_2< Traits > compute_y_at_x_2;
    compute_y_at_x_2.setScene( this->getScene( ) );
    CoordinateType x( this->queryPt.x( ) );
    double yApprox =
      CGAL::to_double( compute_y_at_x_2.approx( halfedge->curve( ), x ) );
    FT yInt( yApprox );
    Kernel_point_2 p2( this->queryPt.x( ), yInt );
    Segment_2 seg( this->queryPt, p2 );
    // QLineF qseg = this->convert( seg );
    // this->activeRay->setLine( qseg );
    this->rayGraphicsItem.setSource( event->scenePos( ) );
    this->rayGraphicsItem.setTargetY( CGAL::to_double( yInt ) );
    this->rayGraphicsItem.setIsInfinite( false );
  }
  else if ( CGAL::assign( vertex, pointLocationResult ) )
  {
    if ( ! vertex->is_at_open_boundary( ) )
    {
      Point_2 pt = vertex->point( );
      this->highlightedCurves->insert( pt );
    }
  }

  emit modelChanged( );
}

template < typename Arr_ >
typename VerticalRayShootCallback< Arr_ >::Face_const_handle
VerticalRayShootCallback< Arr_ >::getFace( const CGAL::Object& obj )
{
  Face_const_handle f;
  if ( CGAL::assign( f, obj ) )
    return f;

  Halfedge_const_handle he;
  if (CGAL::assign( he, obj ))
    return (he->face( ));

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign( v, obj ));
  CGAL::assign( v, obj );
  if ( v->is_isolated( ) )
    return v->face( );
  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges( );
  return  (eit->face( ));
}

template < typename Arr_ >
CGAL::Object
VerticalRayShootCallback< Arr_ >::rayShootUp( const Kernel_point_2& pt )
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy* landmarksStrategy;

  Point_2 point = this->toArrPoint( pt );

  if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = walkStrategy->ray_shoot_up( point );
  }
  else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_up( point );
  }
  else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = simpleStrategy->ray_shoot_up( point );
  }
  else if ( CGAL::assign( landmarksStrategy, this->pointLocationStrategy ) )
  {
    // pointLocationResult = landmarksStrategy->locate( point );
    std::cerr << "Warning: landmarks point location strategy doesn't support ray shooting" << std::endl;
    return CGAL::Object( );
  }
  return pointLocationResult;
}

template < typename Arr_ >
CGAL::Object
VerticalRayShootCallback< Arr_ >::rayShootDown( const Kernel_point_2& pt )
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy* landmarksStrategy;

  Point_2 point = this->toArrPoint( pt );

  if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = walkStrategy->ray_shoot_down( point );
  }
  else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = trapezoidStrategy->ray_shoot_down( point );
  }
  else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = simpleStrategy->ray_shoot_down( point );
  }
  else if ( CGAL::assign( landmarksStrategy, this->pointLocationStrategy ) )
  {
    // pointLocationResult = landmarksStrategy->locate( point );
    std::cerr << "Warning: landmarks point location strategy doesn't support ray shooting" << std::endl;
    return CGAL::Object( );
  }
  return pointLocationResult;
}

#endif // VERTICAL_RAY_SHOOT_CALLBACK_H
