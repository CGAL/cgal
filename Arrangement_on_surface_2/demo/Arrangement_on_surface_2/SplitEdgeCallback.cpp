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

#include "SplitEdgeCallback.h"

//! enable/disable snapping
/*!
  \param b boolean value to toggle the its currents state
*/
void SplitEdgeCallbackBase::setSnappingEnabled( bool b )
{
  this->snappingEnabled = b;
}

//! enable/disable of the snapping grid
/*!
  \param b boolean value to toggle its current state
*/
void SplitEdgeCallbackBase::setSnapToGridEnabled( bool b )
{
  this->snapToGridEnabled = b;
}

//! displays the color of the nodes and edges after the splitting has been performed
/*!
  \param c A QColor object for the selected scene
*/
void SplitEdgeCallbackBase::setColor( QColor c )
{
  this->color = c;
}

//! returns the color of nodes and the edges which is splitted
/*!
  \return a QColor object 
*/
QColor SplitEdgeCallbackBase::getColor( ) const
{
  return this->color;
}

SplitEdgeCallbackBase::SplitEdgeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent ),
  snappingEnabled( false ),
  snapToGridEnabled( false ),
  color( ::Qt::blue )
{ }


template <typename Arr_>
SplitEdgeCallback<Arr_>::SplitEdgeCallback(Arrangement* arr_, QObject* parent):
  SplitEdgeCallbackBase( parent ),
  arr( arr_ ),
  hasFirstPoint( false ),
  intersectCurves( this->traits.intersect_2_object( ) ),
  areEqual( this->traits.equal_2_object( ) ),
  segmentGuide(new QGraphicsLineItem())
{
  this->segmentGuide->setZValue( 100 );
  QPen pen = this->segmentGuide->pen( );
  pen.setColor( this->color );
  pen.setCosmetic(true);
  this->segmentGuide->setPen( pen );
  this->snapToVertexStrategy.setArrangement( arr_ );

  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this, SLOT( slotModelChanged( ) ) );
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
  this->snapToVertexStrategy.setScene( scene_ );
  this->snapToGridStrategy.setScene( scene_ );
  if ( this->scene )
  {
    this->scene->addItem( this->segmentGuide );
  }
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setColor( QColor c )
{
  this->color = c;

  QPen pen = this->segmentGuide->pen( );
  pen.setColor( c );
  this->segmentGuide->setPen( pen );
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::reset( )
{
  this->hasFirstPoint = false;
  this->segmentGuide->setLine(0,0,0,0);
  Q_EMIT modelChanged( );
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::slotModelChanged( )
{
  this->segmentGuide->update( );
}

template <typename Arr_>
void
SplitEdgeCallback<Arr_>::mousePressEvent( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->snapPoint( event );
  this->splitEdges( clickedPoint, Traits( ) );
}

template <typename Arr_>
template < typename TTraits >
void SplitEdgeCallback<Arr_>::
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

  Q_EMIT modelChanged( );
}

template <typename Arr_>
template < typename CircularKernel >
void SplitEdgeCallback<Arr_>::
splitEdges(const Point_2& clickedPoint,
           CGAL::Arr_circular_arc_traits_2<CircularKernel> /* traits */)
{
  typedef CircularKernel Kernel;
  typedef typename Kernel::Point_2 Ker_Point_2;
  typedef typename Kernel::Line_arc_2 Line_arc_2;

  if ( ! this->hasFirstPoint )
  {
    this->p1 = clickedPoint;
    this->hasFirstPoint = true;
  }
  else
  {
    this->p2 = clickedPoint;
    Ker_Point_2 ker_p1((CGAL::to_double(this->p1.x())), CGAL::to_double((this->p1.y())));
    Ker_Point_2 ker_p2((CGAL::to_double(this->p2.x())), CGAL::to_double((this->p2.y())));

    Line_arc_2 splitCurve( Segment_2( ker_p1, ker_p2 ) );

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

  Q_EMIT modelChanged( );
}

template <typename Arr_>
template <typename Coefficient_>
void SplitEdgeCallback<Arr_>::
splitEdges(const Point_2& clickedPoint,
           CGAL::Arr_algebraic_segment_traits_2<Coefficient_> /* traits */)
{
  typename Traits::Construct_x_monotone_segment_2 constructSegment =
      traits.construct_x_monotone_segment_2_object( );

  std::vector< X_monotone_curve_2 > curves;

  if ( ! this->hasFirstPoint )
  {
    this->p1 = clickedPoint;
    this->hasFirstPoint = true;
  }
  else
  {
    this->p2 = clickedPoint;
    constructSegment( this->p1, this->p2, std::back_inserter( curves ) );

    X_monotone_curve_2 splitCurve = curves[0];

    for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
          hei != this->arr->halfedges_end( ); ++hei )
    {
      X_monotone_curve_2 curve = hei->curve( );
      CGAL::Object res;
      CGAL::Oneset_iterator< CGAL::Object > oi( res );
      this->intersectCurves( splitCurve, curve, oi );
      std::pair< Point_2, Multiplicity > pair;
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

  Q_EMIT modelChanged( );
}

template <typename Arr_>
void
SplitEdgeCallback<Arr_>::mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->snapPoint( event );
  this->updateGuide( clickedPoint, this->traits );
}

template <typename Arr_>
template < typename TTraits >
void SplitEdgeCallback<Arr_>::
updateGuide(const Point_2& clickedPoint,
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

    this->segmentGuide->setLine( qSegment );
    Q_EMIT modelChanged( );
  }
}

template <typename Arr_>
template <typename CircularKernel>
void SplitEdgeCallback<Arr_>::
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

    this->segmentGuide->setLine( qSegment );
    Q_EMIT modelChanged( );
  }
}

template <typename Arr_>
typename SplitEdgeCallback<Arr_>::Point_2
SplitEdgeCallback<Arr_>::snapPoint( QGraphicsSceneMouseEvent* event )
{
  return this->snapPoint( event, Traits( ) );
}

template <typename Arr_>
template <typename TTraits>
typename SplitEdgeCallback<Arr_>::Point_2
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

template <typename Arr_>
template <typename CircularKernel>
typename SplitEdgeCallback<Arr_>::Point_2 SplitEdgeCallback<Arr_>::
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

template class SplitEdgeCallback<Seg_arr>;
template class SplitEdgeCallback<Pol_arr>;
template class SplitEdgeCallback<Conic_arr>;
template class SplitEdgeCallback<Lin_arr>;
template class SplitEdgeCallback<Arc_arr>;
template class SplitEdgeCallback<Alg_seg_arr>;
