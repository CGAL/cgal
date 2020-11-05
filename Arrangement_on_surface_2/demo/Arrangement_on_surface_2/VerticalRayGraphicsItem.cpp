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

#include "VerticalRayGraphicsItem.h"

#include <limits>
#include <QPainter>
#include <QGraphicsView>
#include <QScrollBar>

VerticalRayGraphicsItem::VerticalRayGraphicsItem( ) :
  m_source( QPointF( ) ), // null point ie. (+0.0, +0.0)
  m_targetY( 0.0 ),
  m_isInfinite( false ),
  m_color( ::Qt::green ),
  m_width( 1 )
{}

void VerticalRayGraphicsItem::paint( QPainter* painter,
                                     const QStyleOptionGraphicsItem* /* option*/,
                                     QWidget* /* widget */ )
{
  QPen rayPen( this->m_color, this->m_width );
  rayPen.setCosmetic(true);
  painter->setPen( rayPen );

  if ( this->m_source.isNull( ) && this->m_targetY == 0.0 )
  {
    return;
  }
  // -y towards the top
  bool isShootingUp = ( this->m_source.y( ) < this->m_targetY );
  if ( this->m_isInfinite )
  {
    QRectF clipRect = this->viewportRect( );
    double targetY = this->m_targetY;
    if ( isShootingUp && this->m_isInfinite )
    {
      targetY = clipRect.top( );
    }
    else if ( this->m_isInfinite )
    {
      targetY = clipRect.bottom( );
    }
    QPointF target( this->m_source.x( ), targetY );
    QLineF line( this->m_source, target );
    painter->drawLine( line );
    // TODO: draw arrowhead
    this->drawArrowhead( painter, targetY, isShootingUp );
    // std::cout << "drawing the ray to " << targetY << std::endl;
  }
  else
  {
    QPointF target( this->m_source.x( ), this->m_targetY );
    QLineF line( this->m_source, target );
    painter->drawLine( line );
    this->drawArrowhead( painter, this->m_targetY, isShootingUp );
  }
}

QRectF VerticalRayGraphicsItem::boundingRect( ) const
{
  if ( this->m_source.isNull( ) || // uninitialized
       this->m_source.y( ) == this->m_targetY ) // degenerate
  {
    return QRectF( );
  }

  double xmin = this->m_source.x( ) - 5;
  double xmax = this->m_source.x( ) + 5;
  // -y towards the top
  bool isShootingUp = ( this->m_source.y( ) > this->m_targetY );
  double ymin = ( isShootingUp )? this->m_targetY : this->m_source.y( );
  double ymax = ( isShootingUp )? this->m_source.y( ) : this->m_targetY;
  if ( this->m_isInfinite )
  {
    if ( isShootingUp )
    {
      //ymin = -std::numeric_limits< float >::max( );
      ymin = this->m_source.y( ) - 1e9;
    }
    else
    {
      //ymax = std::numeric_limits< float >::max( );
      ymax = this->m_source.y( ) + 1e9;
    }
  }
  // std::cout << ymin << " " << ymax << std::endl;

  QRectF res( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );
  return res;
}

const QPointF& VerticalRayGraphicsItem::source( ) const
{
  return this->m_source;
}

void VerticalRayGraphicsItem::setSource( const QPointF& src )
{
  this->prepareGeometryChange( );

  this->m_source = src;
}

double VerticalRayGraphicsItem::targetY( ) const
{
  return this->m_targetY;
}

void VerticalRayGraphicsItem::setTargetY( double y )
{
  this->prepareGeometryChange( );

  this->m_targetY = y;
}

bool VerticalRayGraphicsItem::isInfinite( ) const
{
  return this->m_isInfinite;
}

void VerticalRayGraphicsItem::setIsInfinite( bool b )
{
  this->prepareGeometryChange( );
  this->m_isInfinite = b;
}

const QColor& VerticalRayGraphicsItem::color( ) const
{
  return this->m_color;
}

void  VerticalRayGraphicsItem::setColor( const QColor& color )
{
  this->m_color = color;
}

int VerticalRayGraphicsItem::width( ) const
{
  return this->m_width;
}

void VerticalRayGraphicsItem::setWidth( int width )
{
  this->m_width = width;
}

void VerticalRayGraphicsItem::reset( )
{
  this->prepareGeometryChange( );

  this->m_source = QPointF( ); // null point ie. (+0.0, +0.0)
  this->m_targetY = ( 0.0 );
  this->m_isInfinite = ( false );
}

void VerticalRayGraphicsItem::modelChanged( )
{
  // std::cout << "ray changed" << std::endl;

  if ( this->m_source.isNull( ) || // uninitialized
       this->m_source.y( ) == this->m_targetY ) // degenerate
  {
    this->hide( );
  }
  else
  {
    this->show( );
  }
  this->update( );
}

// FIXME: there must be a way to get this info from just a QPainter object...
// doing it this way assumes we've only added this item to exactly one view...
QRectF VerticalRayGraphicsItem::viewportRect( ) const
{
  QRectF res;
  if ( this->scene( ) == NULL )
  {
    return res;
  }

  QList< QGraphicsView* > views = this->scene( )->views( );
  if ( views.size( ) == 0 )
  {
    return res;
  }
  // assumes the first view is the right one
  QGraphicsView* viewport = views.first( );
  QPointF p1 = viewport->mapToScene( 0, 0 );
  QPointF p2 = viewport->mapToScene( viewport->width( ), viewport->height( ) );
  res = QRectF( p1, p2 );

  return res;
}

void VerticalRayGraphicsItem::drawArrowhead( QPainter* painter,
                                             double targetY, bool isShootingUp )
{
  if ( this->scene( ) == 0 || this->scene( )->views( ).size( ) == 0 )
  {
    return;
  }
  QGraphicsView* view = this->scene( )->views( ).first( );
  QPointF arrowTip( this->m_source.x( ), targetY );
  QPoint pt = view->mapFromScene( arrowTip );
  if ( ! isShootingUp && this->m_isInfinite )
  {
    if (view->horizontalScrollBar( ) &&
        view->horizontalScrollBar( )->isVisible( ) )
    {
      // std::cout << view->horizontalScrollBar( )->height( ) << std::endl;
      pt.setY( pt.y( ) - view->horizontalScrollBar( )->height( ) - 5 );
      arrowTip = view->mapToScene( pt );
    }
    else
    {
      // std::cout << "no scroll bar " << std::endl;
      pt.setY( pt.y( ) - 5 );
      arrowTip = view->mapToScene( pt );
    }
  }

  int dy = -1;
  if ( isShootingUp )
  {
    dy = 1;
  }
  QPoint leftPt( pt.x( ) - 3, pt.y( ) + 3*dy );
  QPoint rightPt( pt.x( ) + 3, pt.y( ) + 3*dy );
  QPointF left = view->mapToScene( leftPt );
  QPointF right = view->mapToScene( rightPt );
  QLineF leftEdge( left, arrowTip );
  QLineF rightEdge( arrowTip, right );
  painter->drawLine( leftEdge );
  painter->drawLine( rightEdge );
}
