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

#ifndef VERTICAL_RAY_GRAPHICS_ITEM_H
#define VERTICAL_RAY_GRAPHICS_ITEM_H

#include <CGAL/Qt/GraphicsItem.h>
#include <QPen>
#include <QColor>

/*
 * Represents a vertical ray in the scene. The ray doesn't necessarily extend
 * to infinity, but it's called a ray because it has an arrowhead at its target
 * end.
 *
 * If it is designated that the ray extends to infinity, we'll clip the ray to
 * the boundary of the visible viewport.
 */
class VerticalRayGraphicsItem : public CGAL::Qt::GraphicsItem
{
public:
  VerticalRayGraphicsItem( );

  virtual void paint( QPainter* painter,
					  const QStyleOptionGraphicsItem* option,
					  QWidget* widget );
  virtual QRectF boundingRect( ) const;				//!< a virtual member function for the bounding box

  const QPointF& source( ) const;
  void setSource( const QPointF& src );
  double targetY( ) const;
  void setTargetY( double y );						//!< direction of the arrow.
  bool isInfinite( ) const;							//!< if it is of infinite length
  void setIsInfinite( bool b );						//!< indefinite length arrow

  const QColor& color( ) const;							
  void setColor( const QColor& color );				//!< setting the color of the arrow.
  int width( ) const;
  void setWidth( int width );						//!< setting the thickness of the arrow.

  void reset( );

public Q_SLOTS:
  virtual void modelChanged( );

protected:
  QRectF viewportRect( ) const;
  void drawArrowhead( QPainter* painter, double targetY, bool isShootingUp );		//!< drawing the arrow on the viewport.

  QPointF m_source;                   				/*!< position of the arrow */
  double m_targetY;									
  bool m_isInfinite;								/*!< if the arrow is set to infinity */
  QColor m_color;									/*!< color of the arrow */
  int m_width;
}; // class VerticalRayGraphicsItem

#endif // VERTICAL_RAY_GRAPHICS_ITEM_H
