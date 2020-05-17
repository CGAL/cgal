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

#include "ArrangementDemoGraphicsView.h"

#include <iostream>
#include <QVarLengthArray>
#include <QPen>
#include <QCoreApplication>
#include <QKeyEvent>

//! Member function to setup the viewport of the screen
/*!
  \param parent a Qwidget pointer to the class
*/
ArrangementDemoGraphicsView::ArrangementDemoGraphicsView( QWidget* parent ) :
  QGraphicsView( parent ),
  showGrid( false ),
  gridSize( 50 ),
  gridColor( ::Qt::black ),
  backgroundColor( ::Qt::white )
{
  QTransform m( 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 );
  this->setTransform( m );
  this->setBackgroundBrush( QBrush( backgroundColor ) );
  this->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
  this->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
}

void ArrangementDemoGraphicsView::wheelEvent(QWheelEvent* event)
{
  // std::cout<<"In ArrangementDemoGraphicsView wheelEvent\n";
  // this->centerOn(this->mapToScene(event->pos()));
}

//! Member function to display the grid
/*!
  \param b boolean value to toggle the grid view
*/
void ArrangementDemoGraphicsView::setShowGrid( bool b )
{
  this->showGrid = b;
}

//! Member function to get the grid
/*!
  \return the grid to be displayed
*/
bool ArrangementDemoGraphicsView::getShowGrid( ) const
{
  return this->showGrid;
}

//! Member function to set the size of grid
/*!
  \param size integer variable 
*/
void ArrangementDemoGraphicsView::setGridSize( int size )
{
  this->gridSize = size;
}

//! Member function to create the user screen
/*!
  \return The screen to be manipulated
*/
int ArrangementDemoGraphicsView::getGridSize( ) const
{
  return this->gridSize;
}

//! Member function to set the color of grid
/*!
  \param color a QColor object 
*/
void ArrangementDemoGraphicsView::setGridColor( QColor color )
{
  this->gridColor = color;
}

//! Member function to get the color of grid
/*!
  \return the grid color being set
*/
QColor ArrangementDemoGraphicsView::getGridColor( ) const
{
  return this->gridColor;
}

//! Member function to set color of background
/*!
  \param color a Qcolor object
*/
void ArrangementDemoGraphicsView::setBackgroundColor( QColor color )
{
  this->backgroundColor = color;
}

//! Member function to get back the color of background screen
/*!
  \return The of the viewport/screen
*/
QColor ArrangementDemoGraphicsView::getBackgroundColor( ) const
{
  return this->backgroundColor;
}

//! Member function to draw within the size of the viewport even when its changed
/*!
  \param painter a Qpainter pointer to its class 
*/
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

//! Member function to create the user screen
/*!
  \return The screen to be manipulated
*/
QRectF ArrangementDemoGraphicsView::getViewportRect( ) const
{
  QPointF p1 = this->mapToScene( 0, 0 );
  QPointF p2 = this->mapToScene( this->width( ), this->height( ) );

  double xmin = (std::min)( p1.x( ), p2.x( ) );
  double xmax = (std::max)( p1.x( ), p2.x( ) );
  double ymin = (std::min)( p1.y( ), p2.y( ) );
  double ymax = (std::max)( p1.y( ), p2.y( ) );

  QRectF res = QRectF( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );

  return res;
}
