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

#ifndef ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
#define ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H

#include <CGAL/Qt/GraphicsViewInput.h>
//#include <CGAL/Qt/Converter.h>
//#include <QEvent>
//#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <iostream>

#include "GraphicsViewSegmentInput.h"
#include "Utils.h"

template < typename Arr_ >
class ArrangementSegmentInputCallback:
  public CGAL::Qt::GraphicsViewSegmentInput<
    typename ArrTraitsAdaptor< typename Arr_::Geometry_traits_2 >::Kernel >
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Vertex_iterator Vertex_iterator;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef CGAL::Qt::GraphicsViewSegmentInput< Kernel > Superclass;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Construct_x_monotone_curve_2
  Construct_x_monotone_curve_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::FT FT;

  ArrangementSegmentInputCallback( Arrangement* arrangement_, QObject* parent);
  void processInput( CGAL::Object );
  void setScene( QGraphicsScene* scene_ );

protected:
  Point_2 snapPoint( QGraphicsSceneMouseEvent* event );

  Arrangement* arrangement;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
  SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
  SnapToGridStrategy< Kernel > snapToGridStrategy;
}; // class ArrangementSegmentInputCallback

template < typename Arr_ >
ArrangementSegmentInputCallback< Arr_ >::
ArrangementSegmentInputCallback( Arrangement* arrangement_, QObject* parent ):
  Superclass( parent ),
  arrangement( arrangement_ )
{
  this->snapToVertexStrategy.setArrangement( arrangement_ );

  QObject::connect( this, SIGNAL( generate( CGAL::Object ) ),
                    this, SLOT( processInput( CGAL::Object ) ) );
}

template < typename Arr_ >
void
ArrangementSegmentInputCallback< Arr_ >::
processInput( CGAL::Object o )
{
  Segment_2 segment;
  if ( CGAL::assign( segment, o ) )
  {
    // insert a segment
    Point_2 p1 = segment.source( );
    Point_2 p2 = segment.target( );
    Curve_2 curve = this->construct_x_monotone_curve_2( p1, p2 );
    CGAL::insert( *( this->arrangement ), curve );
  }
    
  emit CGAL::Qt::GraphicsViewInput::modelChanged( );
}

template < typename Arr_ >
void
ArrangementSegmentInputCallback< Arr_ >::setScene( QGraphicsScene* scene_ )
{
  this->Superclass::setScene( scene_ );
  this->snapToVertexStrategy.setScene( scene_ );
  this->snapToGridStrategy.setScene( scene_ );
}

template < typename Arr_ >
typename ArrangementSegmentInputCallback< Arr_ >::Point_2
ArrangementSegmentInputCallback< Arr_ >::
snapPoint( QGraphicsSceneMouseEvent* event )
{
  if ( this->snapToGridEnabled )
  {
    return this->snapToGridStrategy.snapPoint( event );
  }
  else if ( this->snappingEnabled )
  {
    return this->snapToVertexStrategy.snapPoint( event );
  }
  else
  {
    return this->convert( event->scenePos( ) );
  }
}

#endif // ARRANGEMENT_SEGMENT_INPUT_CALLBACK_H
