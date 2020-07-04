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

#ifndef ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
#define ARRANGEMENT_DEMO_GRAPHICS_VIEW_H

#include <QGraphicsView>
#include <QColor>

class QPaintEvent;

class ArrangementDemoGraphicsView : public QGraphicsView
{
public:
  ArrangementDemoGraphicsView( QWidget* parent = 0 );

  void setBackgroundColor( QColor color );
  QColor getBackgroundColor( ) const;
  QRectF viewportRect() const;
  static QRectF viewportRect(const QGraphicsView* view);
  void resetTransform();

protected:
  void drawBackground(QPainter* painter, const QRectF& rect) override;
  void paintEvent(QPaintEvent* event) override;

  qreal maxScale;
  qreal minScale;
  QColor backgroundColor;     /*!< color for the background */
};

#endif // ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
