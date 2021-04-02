// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

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
  maxScale( 500000 ),
  minScale( 0.0002 )
{
  this->resetTransform();
  this->setResizeAnchor(QGraphicsView::AnchorUnderMouse);
  this->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
  this->setMouseTracking( true );
  this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  this->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  // TODO: Make options menu work
  this->setRenderHint(QPainter::Antialiasing);
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

void ArrangementDemoGraphicsView::resetTransform()
{
  this->setTransform({1.0, 0.0, 0.0, -1.0, 0.0, 0.0});
}
