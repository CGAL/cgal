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

#ifndef POINTS_GRAPHICS_ITEM_H
#define POINTS_GRAPHICS_ITEM_H

#include <vector>
#include <CGAL/Qt/GraphicsItem.h>
#include <QPen>

class QPainter;
class QPen;

/**
   Add a set of points to the QGraphicsScene.
*/
class PointsGraphicsItem: public CGAL::Qt::GraphicsItem
{
public:
  PointsGraphicsItem( );

  virtual void paint( QPainter* painter,
                      const QStyleOptionGraphicsItem* option, QWidget* widget );
  virtual QRectF boundingRect( ) const;

  template < class Point >
  void insert( const Point& point )
  {
    this->prepareGeometryChange( );

    double x = CGAL::to_double( point.x( ) );
    double y = CGAL::to_double( point.y( ) );
    this->points.push_back( QPointF( x, y ) );
  }

  void clear( );

  void setColor( QColor c );
  QColor getColor( ) const;

  void setPointRadius( double d );
  double getPointRadius( ) const;

public slots:
  virtual void modelChanged( );

protected:
  std::vector< QPointF > points;
  double pointRadius;
  QColor color;

}; // class PointsGraphicsItem

#endif // POINTS_GRAPHICS_ITEM_H
