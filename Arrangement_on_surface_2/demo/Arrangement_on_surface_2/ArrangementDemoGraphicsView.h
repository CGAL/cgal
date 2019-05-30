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

class ArrangementDemoGraphicsView : public QGraphicsView
{
public:
  ArrangementDemoGraphicsView( QWidget* parent = 0 );

  void setShowGrid( bool b );
  bool getShowGrid( ) const;
  void setGridSize( int size );
  int getGridSize( ) const;
  void setGridColor( QColor color );
  QColor getGridColor( ) const;
  void setBackgroundColor( QColor color );
  QColor getBackgroundColor( ) const;


protected:
  void wheelEvent(QWheelEvent* event);
  void drawForeground( QPainter* painter, const QRectF& rect );     //!< drawing on the screen
  QRectF getViewportRect( ) const;                                  //!< get current screen size

  bool showGrid;            /*!< displaying grid toggle */
  int gridSize;              /*!< an integer value for the grid size */
  QColor gridColor;           /*!< color for the grid */
  QColor backgroundColor;     /*!< color for the background */
};

#endif // ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
