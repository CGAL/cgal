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

#ifndef POINT_LOCATION_CALLBACK_H
#define POINT_LOCATION_CALLBACK_H

#include "Callback.h"
#include <QEvent>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include "CurveGraphicsItem.h"
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_tags.h>

#include "Utils.h"

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
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename CGAL::Arr_trapezoid_ric_point_location< Arrangement >
    TrapezoidPointLocationStrategy;
  typedef typename CGAL::Arr_simple_point_location< Arrangement >
    SimplePointLocationStrategy;
  typedef typename CGAL::Arr_walk_along_line_point_location< Arrangement >
    Walk_pl_strategy;
  typedef typename Supports_landmarks< Arrangement >::LandmarksType
    LandmarksPointLocationStrategy;

  PointLocationCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );
  void setScene( QGraphicsScene* scene_ );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event,
                               CGAL::Arr_oblivious_side_tag );
  void highlightPointLocation( QGraphicsSceneMouseEvent *event,
                               CGAL::Arr_open_side_tag );
  Face_const_handle getFace( const CGAL::Object& o );
  CGAL::Object locate( const Kernel_point_2& point );
  CGAL::Object locate( const Kernel_point_2& point,
                       CGAL::Tag_false /*supportsLandmarks*/ );
  CGAL::Object locate( const Kernel_point_2& point,
                       CGAL::Tag_true /*doesNotSupportLandmarks*/ );

  using Callback::scene;
  CGAL::Qt::Converter< Kernel > convert;
  CGAL::Object pointLocationStrategy;
  Arrangement* arr;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurves;
  Arr_construct_point_2< Traits > toArrPoint;
}; // class PointLocationCallback

/*! Constructor */
template < typename Arr_ >
PointLocationCallback< Arr_ >::
PointLocationCallback( Arrangement* arr_, QObject* parent_ ) :
  CGAL::Qt::Callback( parent_ ),
  pointLocationStrategy( CGAL::make_object( new Walk_pl_strategy( *arr_ ) ) ),
  arr( arr_ ),
  highlightedCurves( new CGAL::Qt::CurveGraphicsItem< Traits >( ) )
{ 
  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this->highlightedCurves, SLOT( modelChanged( ) ) );
}

template < typename Arr_ >
void
PointLocationCallback< Arr_ >::
setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
  this->highlightedCurves->setScene( scene_ );
  if ( this->scene )
  {
    this->scene->addItem( this->highlightedCurves );
  }
}

template < typename Arr_ >
void
PointLocationCallback< Arr_ >::
reset( )
{
  this->highlightedCurves->clear( );
  emit modelChanged( );
}

template < typename Arr_ >
void 
PointLocationCallback< Arr_ >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
  this->highlightPointLocation( event );
}

template < typename Arr_ >
void PointLocationCallback< Arr_ >::
mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */)
{ }

template < typename Arr_ >
void PointLocationCallback< Arr_ >::
highlightPointLocation( QGraphicsSceneMouseEvent* event )
{
  typename Traits::Left_side_category category;
  this->highlightPointLocation( event, category );

  emit modelChanged( );
}

template < typename Arr_ >
void PointLocationCallback< Arr_ >::
highlightPointLocation( QGraphicsSceneMouseEvent *event,
                        CGAL::Arr_oblivious_side_tag )
{
  Kernel_point_2 point = this->convert( event->scenePos( ) );

  CGAL::Object pointLocationResult = this->locate( point );
  Face_const_handle face = this->getFace( pointLocationResult );
  this->highlightedCurves->clear( );
  if ( ! face->is_unbounded( ) )
  { // it is an interior face; highlight its border
    Ccb_halfedge_const_circulator cc = face->outer_ccb( );
    do
    {
      X_monotone_curve_2 curve = cc->curve( );
      this->highlightedCurves->insert( curve );
    } while ( ++cc != face->outer_ccb( ) );
  }
  Hole_const_iterator hit; 
  Hole_const_iterator eit = face->holes_end( );
  for ( hit = face->holes_begin( ); hit != eit; ++hit )
  { // highlight any holes inside this face
    Ccb_halfedge_const_circulator cc = *hit;
    do
    {
      X_monotone_curve_2 curve = cc->curve( );
      this->highlightedCurves->insert( curve );
      cc++;
    }
    while ( cc != *hit );
  }
}

template < typename Arr_ >
void PointLocationCallback< Arr_ >::
highlightPointLocation( QGraphicsSceneMouseEvent *event,
                        CGAL::Arr_open_side_tag )
{
  Kernel_point_2 point = this->convert( event->scenePos( ) );
  CGAL::Object pointLocationResult = this->locate( point );
  Face_const_handle face = this->getFace( pointLocationResult );
  this->highlightedCurves->clear( );
  Ccb_halfedge_const_circulator cc = face->outer_ccb( );
  do
  {
    if ( ! cc->is_fictitious( ) )
    {
      X_monotone_curve_2 curve = cc->curve( );
      this->highlightedCurves->insert( curve );
    }
  } while ( ++cc != face->outer_ccb( ) );
  Hole_const_iterator hit; 
  Hole_const_iterator eit = face->holes_end( );
  for ( hit = face->holes_begin( ); hit != eit; ++hit )
  { // highlight any holes inside this face
    Ccb_halfedge_const_circulator cc = *hit;
    do
    {
      X_monotone_curve_2 curve = cc->curve( );
      this->highlightedCurves->insert( curve );
      cc++;
    }
    while ( cc != *hit );
  }
}

template < typename Arr_ >
typename PointLocationCallback< Arr_ >::Face_const_handle
PointLocationCallback< Arr_ >::getFace( const CGAL::Object& obj )
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
CGAL::Object PointLocationCallback<Arr_>::locate(const Kernel_point_2& point)
{
  typename Supports_landmarks< Arrangement >::Tag supportsLandmarks;
  return this->locate( point, supportsLandmarks );
}

template < typename Arr_ >
CGAL::Object PointLocationCallback< Arr_ >::locate( const Kernel_point_2& pt,
                                                    CGAL::Tag_true )
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;
  LandmarksPointLocationStrategy* landmarksStrategy;

  Point_2 point = this->toArrPoint( pt );

  if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = walkStrategy->locate( point );
  }
  else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = trapezoidStrategy->locate( point );
  }
  else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = simpleStrategy->locate( point );
  }
  else if ( CGAL::assign( landmarksStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = landmarksStrategy->locate( point );
  }
  return pointLocationResult;
}

template < typename Arr_ >
CGAL::Object PointLocationCallback< Arr_ >::locate( const Kernel_point_2& pt,
                                                    CGAL::Tag_false )
{
  CGAL::Object pointLocationResult;
  Walk_pl_strategy* walkStrategy;
  TrapezoidPointLocationStrategy* trapezoidStrategy;
  SimplePointLocationStrategy* simpleStrategy;

  Point_2 point = this->toArrPoint( pt );

  if ( CGAL::assign( walkStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = walkStrategy->locate( point );
  }
  else if ( CGAL::assign( trapezoidStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = trapezoidStrategy->locate( point );
  }
  else if ( CGAL::assign( simpleStrategy, this->pointLocationStrategy ) )
  {
    pointLocationResult = simpleStrategy->locate( point );
  }
  return pointLocationResult;
}

#endif // POINT_LOCATION_CALLBACK_H
