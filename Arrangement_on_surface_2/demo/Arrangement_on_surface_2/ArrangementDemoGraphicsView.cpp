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
#include <cmath>
#include <QVarLengthArray>
#include <QPen>
#include <QCoreApplication>
#include <QKeyEvent>
#include <QFontMetrics>

//! Member function to setup the viewport of the screen
/*!
  \param parent a Qwidget pointer to the class
*/
ArrangementDemoGraphicsView::ArrangementDemoGraphicsView( QWidget* parent ) :
  QGraphicsView( parent ),
  maxScale( 1000000 ),
  minScale( 0.0002 ),
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
void ArrangementDemoGraphicsView::drawBackground(
  QPainter* painter, const QRectF& /* rect */)
{
}

void ArrangementDemoGraphicsView::paintEvent(QPaintEvent* event)
{
  qreal scale = std::sqrt(std::abs(this->transform().determinant()));
  if (scale > this->maxScale)
    this->scale(this->maxScale / scale, this->maxScale / scale);
  else if (scale < this->minScale)
    this->scale(this->minScale / scale, this->minScale / scale);
  QGraphicsView::paintEvent(event);
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
