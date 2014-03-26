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

#include "ArrangementDemoGraphicsView.h"

#include <iostream>
#include <QVarLengthArray>
#include <QPen>

ArrangementDemoGraphicsView::ArrangementDemoGraphicsView( QWidget* parent ) :
  QGraphicsView( parent ),
  showGrid( false ),
  gridSize( 50 ),
  gridColor( ::Qt::black ),
  backgroundColor( ::Qt::white )
{
  QMatrix m( 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 );
  this->setMatrix( m );
  this->setBackgroundBrush( QBrush( backgroundColor ) );
}

void ArrangementDemoGraphicsView::setShowGrid( bool b )
{
  this->showGrid = b;
}

bool ArrangementDemoGraphicsView::getShowGrid( ) const
{
  return this->showGrid;
}

void ArrangementDemoGraphicsView::setGridSize( int size )
{
  this->gridSize = size;
}

int ArrangementDemoGraphicsView::getGridSize( ) const
{
  return this->gridSize;
}

void ArrangementDemoGraphicsView::setGridColor( QColor color )
{
  this->gridColor = color;
}

QColor ArrangementDemoGraphicsView::getGridColor( ) const
{
  return this->gridColor;
}

void ArrangementDemoGraphicsView::setBackgroundColor( QColor color )
{
  this->backgroundColor = color;
}

QColor ArrangementDemoGraphicsView::getBackgroundColor( ) const
{
  return this->backgroundColor;
}

void ArrangementDemoGraphicsView::drawForeground( QPainter* painter,
                                                  const QRectF& /* rect */)
{
  QRectF viewportRect = this->getViewportRect( );
  if ( this->showGrid )
  {
    // compute integer-spaced grid lines
    QVarLengthArray< QLineF, 100 > linesX;
    QVarLengthArray< QLineF, 100 > linesY;
    qreal left =
      int(viewportRect.left()) - (int(viewportRect.left()) % this->gridSize);
    qreal top =
      int(viewportRect.top()) - (int(viewportRect.top()) % this->gridSize);
    for ( qreal x = left; x < viewportRect.right( ); x += this->gridSize )
    {
      linesX.append( QLineF(x, viewportRect.top(), x, viewportRect.bottom()));
    }
    for ( qreal y = top; y < viewportRect.bottom( ); y += this->gridSize )
    {
      linesY.append(QLineF(viewportRect.left( ), y, viewportRect.right(), y));
    }

    // set up the painter
    QPen savePen = painter->pen( );
    QPen gridPen( savePen );
    gridPen.setColor( this->gridColor );
    painter->setPen( gridPen );

    // draw the grid
    painter->drawLines( linesX.data( ), linesX.size( ) );
    painter->drawLines( linesY.data( ), linesY.size( ) );

    // revert the painter
    painter->setPen( savePen );
  }
}

QRectF ArrangementDemoGraphicsView::getViewportRect( ) const
{
  QPointF p1 = this->mapToScene( 0, 0 );
  QPointF p2 = this->mapToScene( this->width( ), this->height( ) );

  double xmin = std::min( p1.x( ), p2.x( ) );
  double xmax = std::max( p1.x( ), p2.x( ) );
  double ymin = std::min( p1.y( ), p2.y( ) );
  double ymax = std::max( p1.y( ), p2.y( ) );

  QRectF res = QRectF( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );

  return res;
}
