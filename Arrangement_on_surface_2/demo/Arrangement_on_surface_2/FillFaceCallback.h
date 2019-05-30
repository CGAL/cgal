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

#ifndef FILL_FACE_CALLBACK_H
#define FILL_FACE_CALLBACK_H

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

class FillFaceCallbackBase : public CGAL::Qt::Callback
{
public:
  FillFaceCallbackBase( QObject* parent );

  void setColor( QColor c );
  QColor getColor( ) const;

protected:
  QColor fillColor;                       				/*!< Qcolor object to fill a selected space */
};

/**
   Supports visualization of point location on arrangements.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class Arr_ >
class FillFaceCallback : public FillFaceCallbackBase
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

  FillFaceCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );


protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );

  void fillFace( QGraphicsSceneMouseEvent* event );

  Face_const_handle getFace( const CGAL::Object& o );
  CGAL::Object locate( const Kernel_point_2& point );
  CGAL::Object locate( const Kernel_point_2& point,
					   CGAL::Tag_false/*supportsLandmarks*/ );
  CGAL::Object locate( const Kernel_point_2& point,
					   CGAL::Tag_true /*doesNotSupportLandmarks*/ );

  CGAL::Qt::Converter< Kernel > convert;
  CGAL::Object pointLocationStrategy;
  Arrangement* arr;
  Arr_construct_point_2< Traits > toArrPoint;
}; // class FillFaceCallback

/*! Constructor */
template < class Arr_ >
FillFaceCallback<Arr_>::FillFaceCallback(Arrangement* arr_, QObject* parent_):
  FillFaceCallbackBase( parent_ ),
  pointLocationStrategy( CGAL::make_object( new Walk_pl_strategy( *arr_ ) ) ),
  arr( arr_ )
{ }

template < class Arr_ >
void FillFaceCallback< Arr_ >::reset( )
{
  Q_EMIT modelChanged( );
}

template < class Arr_ >
void FillFaceCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  this->fillFace( event );
  Q_EMIT modelChanged( );

  QGraphicsView* view = this->scene->views( ).first( );
  view->scale(1.01, 1.01);
  view->scale(1/1.01, 1/1.01);
}

template < class Arr_ >
void
FillFaceCallback< Arr_ >::mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */)
{ }

//! A Template type
//! Coloring the closed selected space
/*!
 	\param event A QGrpahicsSceneMouseEvent pointer to the class
*/
template < class Arr_ >
void
FillFaceCallback< Arr_ >::
fillFace( QGraphicsSceneMouseEvent* event )
{
  Kernel_point_2 point = this->convert( event->scenePos( ) );
  CGAL::Object pointLocationResult = this->locate( point );
  Face_const_handle face = this->getFace( pointLocationResult );
  Face_handle f = this->arr->non_const_handle( face );

  if ( f->color() == ::Qt::white && this->fillColor.isValid() )
  {
	f->set_color( this->fillColor );
  }
  else
  {
	f->set_color( ::Qt::white );
  }
}

//! A Template type
//! get the selected face
/*!
 	\param obj A CGAL::Object reference of the face
*/
template < class Arr_ >
typename FillFaceCallback< Arr_ >::Face_const_handle
FillFaceCallback< Arr_ >::getFace( const CGAL::Object& obj )
{
  Face_const_handle f;
  if ( CGAL::assign( f, obj ) )
  {
	return f;
  }

  Halfedge_const_handle he;
  if (CGAL::assign( he, obj ))
  {
	return (he->face( ));
  }

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign( v, obj ));
  CGAL::assign( v, obj );
  if ( v->is_isolated( ) )
  {
	return v->face( );
  }

  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges( );
  return  (eit->face( ));
}

//! A Template type
//! locating the mouse position
/*!
 	\param point A Kernel_point_2 object
*/
template < class Arr_ >
CGAL::Object FillFaceCallback< Arr_ >::locate( const Kernel_point_2& point )
{
  typename Supports_landmarks< Arrangement >::Tag supportsLandmarks;
  return this->locate( point, supportsLandmarks );
}

//! A Template type
//! locating points given a query
/*!
 	\param pt A kernel_point_2 object
 	\param Tag_true A CGAL object available in class
*/
template < class Arr_ >
CGAL::Object
FillFaceCallback< Arr_ >::locate( const Kernel_point_2& pt, CGAL::Tag_true )
{
  CGAL::Object pointLocationResult;									
  Walk_pl_strategy* walkStrategy;									
  TrapezoidPointLocationStrategy* trapezoidStrategy;				/*!< searching of a trapezoid from the given faces */
  SimplePointLocationStrategy* simpleStrategy;						
  LandmarksPointLocationStrategy* landmarksStrategy;				/*!< finds the neares neighbor point */

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

//! A Template type
//! locating points given a query
/*!
 	\param pt A kernel_point_2 object
 	\param Tag_false A CGAL object not available
*/
template < class Arr_ >
CGAL::Object
FillFaceCallback< Arr_ >::locate( const Kernel_point_2& pt, CGAL::Tag_false )
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

#endif // FILL_FACE_CALLBACK_H
