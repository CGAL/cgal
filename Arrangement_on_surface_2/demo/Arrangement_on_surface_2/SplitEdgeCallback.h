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

#ifndef SPLIT_EDGE_CALLBACK_H
#define SPLIT_EDGE_CALLBACK_H

#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arrangement_with_history_2.h>

#include "CurveGraphicsItem.h"
#include "Utils.h"
#include "ISnappable.h"

class SplitEdgeCallbackBase : public CGAL::Qt::Callback, public ISnappable
{
public:
  void setSnappingEnabled( bool b );
  void setSnapToGridEnabled( bool b );

  virtual void setColor( QColor c );
  QColor getColor( ) const;

protected:
  SplitEdgeCallbackBase( QObject* parent );

  bool snappingEnabled;
  bool snapToGridEnabled;
  QColor color;

}; // SplitEdgeCallbackBase

/**
   Handles splitting of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class Arr_ >
class SplitEdgeCallback : public SplitEdgeCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  typedef typename Arrangement::Vertex_iterator Vertex_iterator;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Curve_handle Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
    Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef typename Traits::Intersect_2 Intersect_2;
  typedef typename Traits::Equal_2 Equal_2;
  typedef typename Traits::Multiplicity Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Point_2 Point_2;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::FT FT;

  SplitEdgeCallback( Arrangement* arr_, QObject* parent );
  void setScene( QGraphicsScene* scene_ );
  virtual void setColor( QColor c );
  void reset( );
    
  void slotModelChanged( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );

  virtual Point_2 snapPoint( QGraphicsSceneMouseEvent *event );
  template < class TTraits >
  Point_2 snapPoint( QGraphicsSceneMouseEvent *event, TTraits traits );
  template < class CircularKernel >
  Point_2 snapPoint( QGraphicsSceneMouseEvent* event,
                     CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

  template < class TTraits >
  void splitEdges( const Point_2& pt, TTraits traits );
  template < class CircularKernel >
  void splitEdges( const Point_2& pt,
                   CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

  template < class Coefficient_ >
  void splitEdges( const Point_2& pt,
                   CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits);

  template < class TTraits >
  void updateGuide( const Point_2& pt, TTraits traits );

  template < class CircularKernel >
  void updateGuide( const Point_2& pt,
                    CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

  Traits traits;
  CGAL::Qt::Converter< Kernel > convert;
  Arrangement* arr;
  bool hasFirstPoint;
  Point_2 p1;
  Point_2 p2;
  QGraphicsLineItem segmentGuide;

  Intersect_2 intersectCurves;
  Equal_2 areEqual;
  SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
  SnapToGridStrategy< typename Arrangement::Geometry_traits_2 >
    snapToGridStrategy;
}; // class SplitEdgeCallback

template < typename Arr_ >
SplitEdgeCallback< Arr_ >::SplitEdgeCallback( Arrangement* arr_,
                                              QObject* parent ):
  SplitEdgeCallbackBase( parent ),
  arr( arr_ ),
  hasFirstPoint( false ),
  intersectCurves( this->traits.intersect_2_object( ) ),
  areEqual( this->traits.equal_2_object( ) )
{
  this->segmentGuide.setZValue( 100 );
  QPen pen = this->segmentGuide.pen( );
  pen.setColor( this->color );
  this->segmentGuide.setPen( pen );

  this->snapToVertexStrategy.setArrangement( arr_ );

  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this, SLOT( slotModelChanged( ) ) );
}

template < typename Arr_ >
void SplitEdgeCallback< Arr_ >::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
  this->snapToVertexStrategy.setScene( scene_ );
  this->snapToGridStrategy.setScene( scene_ );
  if ( this->scene )
  {
    this->scene->addItem( &( this->segmentGuide ) );
  }
}

template < typename Arr_ >
void SplitEdgeCallback< Arr_ >::setColor( QColor c )
{
  this->color = c;

  QPen pen = this->segmentGuide.pen( );
  pen.setColor( c );
  this->segmentGuide.setPen( pen );
}

template < typename Arr_ >
void SplitEdgeCallback< Arr_ >::reset( )
{
  this->hasFirstPoint = false;
  this->segmentGuide.setLine( 0, 0, 0, 0 );
  emit modelChanged( );
}

template < typename Arr_ >
void SplitEdgeCallback< Arr_ >::slotModelChanged( )
{
  this->segmentGuide.update( );
}

template < typename Arr_ >
void
SplitEdgeCallback< Arr_ >::mousePressEvent( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->snapPoint( event );
  this->splitEdges( clickedPoint, Traits( ) );
}

template < typename Arr_ >
template < typename TTraits >
void SplitEdgeCallback< Arr_ >::
splitEdges( const Point_2& clickedPoint, TTraits traits )
{
  typename TTraits::Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
    traits.construct_x_monotone_curve_2_object( );
  if ( ! this->hasFirstPoint )
  {
    this->p1 = clickedPoint;
    this->hasFirstPoint = true;
  }
  else
  {
    this->p2 = clickedPoint;
    X_monotone_curve_2 splitCurve =
      construct_x_monotone_curve_2( this->p1, this->p2 );
    for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
          hei != this->arr->halfedges_end( ); ++hei )
    {
      X_monotone_curve_2 curve = hei->curve( );
      CGAL::Object res;
      CGAL::Oneset_iterator< CGAL::Object > oi( res );
      this->intersectCurves( splitCurve, curve, oi );
      std::pair< Point_2, Multiplicity > pair;
      if ( hei == this->arr->halfedges_end( ) )
        continue;
      if ( CGAL::assign( pair, res ) )
      {
        Point_2 splitPoint = pair.first;
        if ( ( ! hei->source( )->is_at_open_boundary( ) &&
               this->areEqual( hei->source( )->point( ), splitPoint ) ) ||
             ( ! hei->target( )->is_at_open_boundary( ) &&
               this->areEqual( hei->target( )->point( ), splitPoint ) ) )
        {
          continue;
        }
        this->arr->split_edge( hei, splitPoint );
      }
    }

    this->reset( );
  }

  emit modelChanged( );
}

template < typename Arr_ >
template < typename CircularKernel >
void SplitEdgeCallback< Arr_ >::
splitEdges(const Point_2& /* clickedPoint */,
           CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
  // std::cout << "Circular arc split edges stub" << std::endl;
}

template < typename Arr_ >
template < typename Coefficient_ >
void SplitEdgeCallback< Arr_ >::
splitEdges(const Point_2& /* clickedPoint */,
           CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
{
  // std::cout << "Algebraic segment split edges stub" << std::endl;
}

template < typename Arr_ >
void 
SplitEdgeCallback< Arr_ >::mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->snapPoint( event );
  this->updateGuide( clickedPoint, this->traits );
}

template < typename Arr_ >
template < typename TTraits >
void SplitEdgeCallback<Arr_>::updateGuide(const Point_2& clickedPoint,
                                          TTraits /* traits */)
{
  if ( this->hasFirstPoint )
  { // provide visual feedback for where the split line is
    Point_2 currentPoint = clickedPoint;
    typename Kernel::Point_2 pt1( CGAL::to_double( this->p1.x( ) ),
                                  CGAL::to_double( this->p1.y( ) ) );
    typename Kernel::Point_2 pt2( CGAL::to_double( currentPoint.x( ) ),
                                  CGAL::to_double( currentPoint.y( ) ) );
    Segment_2 currentSegment( pt1, pt2 );
    QLineF qSegment = this->convert( currentSegment );
    this->segmentGuide.setLine( qSegment );
    emit modelChanged( );
  }
}

template < typename Arr_ >
template < typename CircularKernel >
void SplitEdgeCallback< Arr_ >::
updateGuide(const Point_2& clickedPoint,
            CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
  if ( this->hasFirstPoint )
  { // provide visual feedback for where the split line is
    Point_2 currentPoint = clickedPoint;
    typename CircularKernel::Point_2 pt1( CGAL::to_double( this->p1.x( ) ),
                                          CGAL::to_double( this->p1.y( ) ) );
    typename CircularKernel::Point_2 pt2( CGAL::to_double( currentPoint.x( ) ),
                                          CGAL::to_double( currentPoint.y( ) ) );
    Segment_2 currentSegment( pt1, pt2 );
    QLineF qSegment = this->convert( currentSegment );
    this->segmentGuide.setLine( qSegment );
    emit modelChanged( );
  }
}

template < typename Arr_ >
typename SplitEdgeCallback< Arr_ >::Point_2
SplitEdgeCallback< Arr_ >::snapPoint( QGraphicsSceneMouseEvent* event )
{
  return this->snapPoint( event, Traits( ) );
}

template < typename Arr_ >
template < typename TTraits >
typename SplitEdgeCallback< Arr_ >::Point_2
SplitEdgeCallback<Arr_>::snapPoint(QGraphicsSceneMouseEvent* event,
                                   TTraits /* traits */)
{
  if ( this->snapToGridEnabled )
  {
    return this->snapToGridStrategy.snapPoint( event );
  }
  if ( this->snappingEnabled )
  {
    return this->snapToVertexStrategy.snapPoint( event );
  }
  else
  { // fallback "analog" selection
    typename Kernel::Point_2 pt = this->convert( event->scenePos( ) );
    CoordinateType x( pt.x( ) );
    CoordinateType y( pt.y( ) );
    return Point_2( x, y );
  }
}

template < typename Arr_ >
template < typename CircularKernel >
typename SplitEdgeCallback< Arr_ >::Point_2 SplitEdgeCallback< Arr_ >::
snapPoint(QGraphicsSceneMouseEvent* event,
          CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
{
  if ( this->snapToGridEnabled )
  {
    return this->snapToGridStrategy.snapPoint( event );
  }
  if ( this->snappingEnabled )
  {
    return this->snapToVertexStrategy.snapPoint( event );
  }
  else
  { // fallback "analog" selection
    typename Kernel::Point_2 pt = this->convert( event->scenePos( ) );
    return Point_2( pt );
  }
}

#endif // SPLIT_EDGE_CALLBACK_H
