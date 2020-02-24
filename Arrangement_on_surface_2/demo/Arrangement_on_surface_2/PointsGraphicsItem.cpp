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

#include "PointsGraphicsItem.h"

#include <limits>
#include <QPen>
#include <QPainter>

PointsGraphicsItem::PointsGraphicsItem( ) :
  pointRadius( 3.0 ),
  color( ::Qt::blue )
{ }

void PointsGraphicsItem::paint(QPainter* painter,
                               const QStyleOptionGraphicsItem* /* option */,
                               QWidget* /* widget */)
{
  double scale = painter->worldTransform( ).m11( );
  double radius = this->pointRadius;
  radius /= scale;
  QPen savePen = painter->pen( );
  painter->setBrush( QBrush( this->color ) );

  for ( unsigned int i = 0; i < this->points.size( ); ++i )
  {
    QPointF pt = this->points[ i ];
    painter->drawEllipse( pt, radius, radius );
  }

  painter->setPen( savePen );
}

QRectF PointsGraphicsItem::boundingRect( ) const
{
  if ( this->points.size( ) == 0 )
  {
    return QRectF( );
  }
  double xmin = (std::numeric_limits< double >::max)( );
  double xmax = -(std::numeric_limits< double >::max)( );
  double ymin = (std::numeric_limits< double >::max)( );
  double ymax = -(std::numeric_limits< double >::max)( );
  for ( unsigned int i = 0; i < this->points.size(); ++i )
  {
    QPointF pt = this->points[ i ];
    double x = pt.x( );
    double y = pt.y( );
    xmin = (std::min)( xmin, x );
    xmax = (std::max)( xmax, x );
    ymin = (std::min)( ymin, y );
    ymax = (std::max)( ymax, y );
  }
  QRectF res( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );
  res.adjust( -5, -5, 5, 5 ); // pad the borders a bit
  return res;
}

void PointsGraphicsItem::clear( )
{
  this->prepareGeometryChange( );

  this->points.clear( );
}

void PointsGraphicsItem::setColor( QColor c )
{
  this->color = c;
}

QColor PointsGraphicsItem::getColor( ) const
{
  return this->color;
}

void PointsGraphicsItem::setPointRadius( double d )
{
  this->pointRadius = d;
}

double PointsGraphicsItem::getPointRadius( ) const
{
  return this->pointRadius;
}

void PointsGraphicsItem::modelChanged( )
{
  if ( this->points.size( ) == 0 )
  {
    this->hide( );
  }
  else
  {
    this->show( );
  }
  this->update( );
}
