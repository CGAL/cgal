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
  virtual void setScene( QGraphicsScene* scene_ );  //!< a virtual member setter function
  QGraphicsScene* getScene( ) const;                //!< a getter function.

  void setSnappingEnabled( bool b );                //!< To enable/disable snap
  void setSnapToGridEnabled( bool b );              //!< To enable/disable snaptoGrid

protected:
  GraphicsViewSegmentInputBase( QObject* parent_ );
  virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );       //!< a virtual function for mouse movement
  virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );      //!< a virtual function for mouse actions
  virtual bool eventFilter( QObject* obj, QEvent* event );              //!< a virtual function to filter the mouse actions

  QRectF viewportRect( ) const;                   						//!< returns the bbox of the viewport

  QGraphicsScene* scene;												//!< denoting the current scene
  bool snappingEnabled;													//!< boolean value to enable/disable
  bool snapToGridEnabled;												//!< boolean value to enable/disable

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
  Point_2 p1;												//!< CGAL Point_2 
  Point_2 p2;												//!< CGAL Point_2
  bool second;

  QGraphicsLineItem segmentGuide;
}; // class GraphicsViewSegmentInput


template < class K_ >
GraphicsViewSegmentInput< K_ >::GraphicsViewSegmentInput( QObject* parent ):
  GraphicsViewSegmentInputBase( parent ),
  second( false )
{ }

//! Temple type
//! overridden function of snap point
/*!
  \param event a QGraphicsSceneMouseEvent pointer to the class
  \return the position of the point snapped
*/
template < class K_ >
typename GraphicsViewSegmentInput< K_ >::Point_2
GraphicsViewSegmentInput< K_ >::snapPoint( QGraphicsSceneMouseEvent* event )
{
  Point_2 clickedPoint = this->convert( event->scenePos( ) );
  return clickedPoint;
}

//! Temple type
//! Drawing a line after clicking on the graph scene
/*!
  \param event a QGraphicsSceneMouseEvent pointer to the class
*/
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

//! Temple type
//! adding a segment onto the graph scene
/*!
  \param event a QGraphicsSceneMouseEvent pointer to the class
*/
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
	Q_EMIT generate( CGAL::make_object( res ) );
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
