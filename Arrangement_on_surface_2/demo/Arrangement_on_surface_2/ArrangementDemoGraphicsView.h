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
  void drawForeground( QPainter* painter, const QRectF& rect );
  QRectF getViewportRect( ) const;

  bool showGrid;
  int gridSize;
  QColor gridColor;
  QColor backgroundColor;
};

#endif // ARRANGEMENT_DEMO_GRAPHICS_VIEW_H
