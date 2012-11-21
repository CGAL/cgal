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

#ifndef CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
//#include <iostream>
//#include "Callback.h"
#include "ISnappable.h"

namespace CGAL {
namespace Qt {

class GraphicsViewSegmentInputBase :
    public GraphicsViewInput, public ISnappable
{
public:
  virtual void setScene( QGraphicsScene* scene_ );
  QGraphicsScene* getScene( ) const;

  void setSnappingEnabled( bool b );
  void setSnapToGridEnabled( bool b );

protected:
  GraphicsViewSegmentInputBase( QObject* parent_ );
  virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
  virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
  virtual bool eventFilter( QObject* obj, QEvent* event );

  QRectF viewportRect( ) const;

  QGraphicsScene* scene;
  bool snappingEnabled;
  bool snapToGridEnabled;

}; // class GraphicsViewSegmentInputBase

template < class K_ >
class GraphicsViewSegmentInput: public GraphicsViewSegmentInputBase
{
public:
  typedef K_                            Kernel;
  typedef typename Kernel::Point_2      Point_2;
  typedef typename Kernel::Segment_2    Segment_2;

  GraphicsViewSegmentInput( QObject* parent );

protected:
  void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
  void mousePressEvent( QGraphicsSceneMouseEvent* event );

  // override this to snap to the points you like
  virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event );

  Converter< Kernel > convert;
  Point_2 p1;
  Point_2 p2;
  bool second;

  QGraphicsLineItem segmentGuide;
}; // class GraphicsViewSegmentInput


template < class K_ >
GraphicsViewSegmentInput< K_ >::GraphicsViewSegmentInput( QObject* parent ):
  GraphicsViewSegmentInputBase( parent ),
  second( false )
{ }

template < class K_ >
typename GraphicsViewSegmentInput< K_ >::Point_2
GraphicsViewSegmentInput< K_ >::snapPoint( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->convert( event->scenePos( ) );
  return clickedPoint;
}

template < class K_ >
void
GraphicsViewSegmentInput< K_ >::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  if ( this->second )
  {
    Point_2 clickedPoint = this->snapPoint( event );
    Segment_2 segment( this->p1, clickedPoint );
    QLineF qSegment = this->convert( segment );
    segmentGuide.setLine( qSegment );
  }
}

template < class K_ >
void
GraphicsViewSegmentInput<K_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  if ( !this->second )
  {
    this->second = true;
    this->p1 = this->snapPoint( event );
    QPointF pt = this->convert( this->p1 );
    segmentGuide.setLine( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
    if ( this->scene != NULL )
    {
      this->scene->addItem( &(segmentGuide ) );
    }
  }
  else
  {
    this->second = false;
    this->p2 = this->snapPoint( event );
    if ( this->scene != NULL )
    {
      this->scene->removeItem( &(segmentGuide ) );
    }
    Segment_2 res( this->p1, this->p2 );
    emit generate( CGAL::make_object( res ) );
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
