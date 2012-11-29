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

#ifndef DELETE_CURVE_CALLBACK_H
#define DELETE_CURVE_CALLBACK_H

#include "Callback.h"

//#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include <CGAL/Qt/Converter.h>
//#include <CGAL/Arrangement_with_history_2.h>

#include "CurveGraphicsItem.h"
#include "Utils.h"

/**
   Handles deletion of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < typename Arr_ >
class DeleteCurveCallback : public CGAL::Qt::Callback
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Curve_handle            Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
  Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator   Induced_edge_iterator;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Point;
  typedef typename Kernel::Segment_2                    Segment;

  DeleteCurveCallback( Arrangement* arr_, QObject* parent_ );
  void setScene( QGraphicsScene* scene_ );
  QGraphicsScene* getScene( ) const;
  void reset( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightNearestCurve( QGraphicsSceneMouseEvent *event );

  Compute_squared_distance_2< Traits > squaredDistance;
  CGAL::Qt::Converter< Kernel > convert;
  QGraphicsScene* scene;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve;
  Arrangement* arr;
  Halfedge_handle removableHalfedge;
}; // class DeleteCurveCallback

/*! Constructor */
template < typename Arr_ >
DeleteCurveCallback< Arr_ >::
DeleteCurveCallback( Arrangement* arr_, QObject* parent_ ) :
  CGAL::Qt::Callback( parent_ ),
  scene( NULL ),
  highlightedCurve( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
  arr( arr_ )
{
  QObject::connect( this, SIGNAL( modelChanged( ) ),
                    this->highlightedCurve, SLOT( modelChanged( ) ) );
}

template < typename Arr_ >
void DeleteCurveCallback< Arr_ >::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
  this->highlightedCurve->setScene( scene_ );
  if ( this->scene )
  {
    this->scene->addItem( this->highlightedCurve );
  }
}

template < typename Arr_ >
QGraphicsScene* DeleteCurveCallback< Arr_ >::getScene( ) const
{
  return this->scene;
}

template < typename Arr_ >
void DeleteCurveCallback< Arr_ >::reset( )
{
  this->highlightedCurve->clear( );
  this->removableHalfedge = Halfedge_handle( );
  emit modelChanged( );
}

template < typename Arr_ >
void 
DeleteCurveCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* /* event */)
{
  if ( this->removableHalfedge == Halfedge_handle( ) )
  {
    return;
  }

  bool deleteOriginatingCurve = 1;
  if ( deleteOriginatingCurve )
  {
    Originating_curve_iterator it =
      this->arr->originating_curves_begin( this->removableHalfedge );
    Originating_curve_iterator it_end =
      this->arr->originating_curves_end( this->removableHalfedge );
    while ( it != it_end )
    {
      Originating_curve_iterator temp = it;
      ++temp;
      CGAL::remove_curve( *(this->arr), it );
      it = temp;
    }
  }
  else
  {
    //CGAL::remove_edge( *(this->arr), this->removableHalfedge->curve( ) );
    this->arr->remove_edge( this->removableHalfedge );
  }

  this->reset( );
}

template < typename Arr_ >
void DeleteCurveCallback< Arr_ >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
  this->highlightNearestCurve( event );
}

template < typename Arr_ >
void 
DeleteCurveCallback< Arr_ >::
highlightNearestCurve( QGraphicsSceneMouseEvent* event )
{
  // find the nearest curve to the cursor to be the new highlighted curve
  Point p = this->convert( event->scenePos( ) );
  //bool isFirst = true;
  //double minDist = 0.0;
  //Halfedge_iterator nearestHei;

#if 0
  for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
        hei != this->arr->halfedges_end( );
        ++hei )
  {
    X_monotone_curve_2 curve = hei->curve( );
    double dist = CGAL::to_double( this->squaredDistance( p, curve ) );
    // std::cout << dist << std::endl;
    if ( isFirst || dist < minDist )
    {
      isFirst = false;
      minDist = dist;
      nearestHei = hei;
    }
  }
#endif
  Find_nearest_edge< Arr_ > findNearestEdge( this->arr );
  findNearestEdge.setScene( this->scene );
  Halfedge_const_handle nearestEdge = findNearestEdge( p );
  this->removableHalfedge = this->arr->non_const_handle( nearestEdge );

  // now 'removableHalfedge' holds the closest halfedge to the point of the mouse
  //this->removableHalfedge = nearestHei;
  //if ( isFirst )
  if ( this->removableHalfedge == Halfedge_handle( ) )
  {
    // std::cout << "no curve found" << std::endl;
    return;
  }

  // create a curve graphics item and add it to the scene
  bool deleteOriginatingCurve = 1;
  this->highlightedCurve->clear( );
  if ( deleteOriginatingCurve )
  { // highlight the originating curve
    Originating_curve_iterator ocit, temp;
    ocit = this->arr->originating_curves_begin( this->removableHalfedge );
    while (ocit != this->arr->originating_curves_end(this->removableHalfedge))
    {
      temp = ocit;
      ++temp;

      Curve_handle ch = ocit;
      Induced_edge_iterator itr;
      for ( itr = this->arr->induced_edges_begin( ch );
            itr != this->arr->induced_edges_end( ch );
            ++itr )
      {
        X_monotone_curve_2 curve = (*itr)->curve( );
        this->highlightedCurve->insert( curve );
      }
      ocit = temp;
    }
  }
  else
  { // highlight just the edge
    this->highlightedCurve->insert( this->removableHalfedge->curve( ) );
  }

  emit modelChanged( );
}

#endif // DELETE_CURVE_CALLBACK_H
