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
  this->resetTransform();
  this->setBackgroundBrush( QBrush( backgroundColor ) );
  this->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
  this->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
  this->setMouseTracking( true );
  this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  this->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  // TODO: Make options menu work
  this->setRenderHint(QPainter::Antialiasing);
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
  QRectF viewportRect = this->viewportRect( );
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

QRectF ArrangementDemoGraphicsView::viewportRect() const
{
  return ArrangementDemoGraphicsView::viewportRect(this);
}

void ArrangementDemoGraphicsView::resetTransform()
{
  this->setTransform({1.0, 0.0, 0.0, -1.0, 0.0, 0.0});
}

QRectF
ArrangementDemoGraphicsView::viewportRect(const QGraphicsView* view)
{
  QPointF p1 = view->mapToScene(0, 0);
  QPointF p2 = view->mapToScene(view->width(), view->height());
  // we also need those because view might rotate
  // rotation of rectangle is a parallelogram
  QPointF p3 = view->mapToScene(view->width(), 0);
  QPointF p4 = view->mapToScene(0, view->height());

  double xmin =
    (std::min)((std::min)(p1.x(), p2.x()), (std::min)(p3.x(), p4.x()));
  double xmax =
    (std::max)((std::max)(p1.x(), p2.x()), (std::max)(p3.x(), p4.x()));
  double ymin =
    (std::min)((std::min)(p1.y(), p2.y()), (std::min)(p3.y(), p4.y()));
  double ymax =
    (std::max)((std::max)(p1.y(), p2.y()), (std::max)(p3.y(), p4.y()));

  QRectF res = QRectF(QPointF(xmin, ymin), QPointF(xmax, ymax));
  return res;
}
