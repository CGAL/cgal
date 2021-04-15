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

#include "ArrangementGraphicsItem.h"

namespace CGAL {
namespace Qt {

ArrangementGraphicsItemBase::ArrangementGraphicsItemBase( ) :
  bb( 0, 0, 0, 0 ),
  bb_initialized( false ),
  verticesPen( QPen( ::Qt::blue, 3. ) ),
  edgesPen( QPen( ::Qt::blue, 1. ) ),
  visible_edges( true ),
  visible_vertices( true ),
  scene( NULL ),
  backgroundColor( ::Qt::white )
{
  this->verticesPen.setCosmetic( true );
  this->verticesPen.setCapStyle( ::Qt::SquareCap );
  this->edgesPen.setCosmetic( true );
}

const QPen& ArrangementGraphicsItemBase::getVerticesPen( ) const
{
  return this->verticesPen;
}

const QPen& ArrangementGraphicsItemBase::getEdgesPen( ) const
{
  return this->edgesPen;
}

void ArrangementGraphicsItemBase::setVerticesPen( const QPen& pen )
{
  this->verticesPen = pen;
}

void ArrangementGraphicsItemBase::setEdgesPen( const QPen& pen )
{
  this->edgesPen = pen;
}

bool ArrangementGraphicsItemBase::visibleVertices( ) const
{
  return this->visible_vertices;
}

void ArrangementGraphicsItemBase::setVisibleVertices( const bool b )
{
  this->visible_vertices = b;
  this->update( );
}

bool ArrangementGraphicsItemBase::visibleEdges( ) const
{
  return this->visible_edges;
}

void ArrangementGraphicsItemBase::setVisibleEdges( const bool b )
{
  this->visible_edges = b;
  this->update( );
}

void ArrangementGraphicsItemBase::setBackgroundColor( QColor color )
{
  this->backgroundColor = color;
}

#if 0
void ArrangementGraphicsItemBase::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
}

QRectF ArrangementGraphicsItemBase::getViewportRect( ) const
{
  QRectF clipRect;
  if ( this->scene == NULL || this->scene->views( ).size( ) == 0 )
  {
    return clipRect;
  }

  QGraphicsView* view = this->scene->views( ).first( );
  QPointF p1 = view->mapToScene( 0, 0 );
  QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
  clipRect = QRectF( p1, p2 );

  return clipRect;
}

#endif

} // namespace Qt
} // namespace CGAL
