// Copyright (c) 2010  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/GraphicsView/include/CGAL/Qt/GpsCircleSegmentGraphicsItem.h $
// $Id: GpsCircleSegmentGraphicsItem.h 45725 2008-09-24 13:24:24Z lrineau $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
//                 Fernando Cacciola <Fernando.Cacciola@geometryfactory.com>

#ifndef CGAL_QT_GPS_CIRCLE_SEGMENT_GRAPHICS_ITEM_H
#define CGAL_QT_GPS_CIRCLE_SEGMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename CS>
class GpsCircleSegmentGraphicsItem : public GraphicsItem
{
  typedef CS Circle_segment_2 ;

  typedef Qt::Converter< Simple_cartesian<double> > Converter ;

    public:
      
  GpsCircleSegmentGraphicsItem();

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return this->vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    this->vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }
  
  void setCircleSegment(const Circle_segment_2& cs);

  Circle_segment_2 circleSegment() const
  {
    return cs_;
  }

protected:
  void updateBoundingBox();

  QRectF bounding_rect;

  QPen edges_pen;

  Circle_segment_2 cs_;
};


template <typename CS>
void 
GpsCircleSegmentGraphicsItem<CS>::setCircleSegment(const Circle_segment_2& cs)
{
  cs_ = cs;
  updateBoundingBox();
  update();
}

template <typename CS>
GpsCircleSegmentGraphicsItem<CS>::GpsCircleSegmentGraphicsItem()
{
  this->hide();
  setZValue(3);
}

template <typename CS>
QRectF 
GpsCircleSegmentGraphicsItem<CS>::boundingRect() const
{
  return bounding_rect;
}

template <typename CS>
void 
GpsCircleSegmentGraphicsItem<CS>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  painter->setPen(edgesPen());

  if ( !cs_.is_full() )
  {
    double sx = to_double(cs_.source().x());
    double sy = to_double(cs_.source().y());
    double tx = to_double(cs_.target().x());
    double ty = to_double(cs_.target().y());

    if( cs_.orientation() == COLLINEAR)
    {
      painter->drawLine(sx,sy,tx,ty);
    }
    else
    {
      double cx = to_double(cs_.supporting_circle().center().x());
      double cy = to_double(cs_.supporting_circle().center().y());

      double x0, y0, x1, y1 ;
      if(cs_.orientation() == CLOCKWISE)
      {
        x0 = sx ;
        y0 = sy ;
        x1 = tx ;
        y1 = ty ; 
      }
      else
      {
        x0 = tx ;
        y0 = ty ;
        x1 = sx ;
        y1 = sy ; 
      }
      double rad = std::sqrt(CGAL::to_double(cs_.supporting_circle().squared_radius()));

      double a   = std::atan2( y0 - cy, x0 - cx ) ;
      double a2p = std::atan2( y1 - cy, x1 - cx );

      if (a2p <= a)
        a2p += 2 * CGAL_PI;

      double alen2 = a2p - a;

      double to_deg = 180 / CGAL_PI ;

      painter->drawArc(cx - rad, cy - rad, 2 * rad, 2 * rad, - a * to_deg * 16, - alen2 * to_deg * 16 );
    }
  }
  else
  {
    if( cs_.orientation() != COLLINEAR)
    {
      double cx = to_double(cs_.supporting_circle().center().x());
      double cy = to_double(cs_.supporting_circle().center().y());
      double rad = std::sqrt(CGAL::to_double(cs_.supporting_circle().squared_radius()));
      painter->drawArc(cx -rad, cy - rad, 2 * rad, 2 * rad, 0, 360*16);
    }
  }


      
}

template <typename CS>
void 
GpsCircleSegmentGraphicsItem<CS>::updateBoundingBox()
{
  Converter convert;
      
  prepareGeometryChange();

  if(cs_.is_circular())
  {
    bounding_rect = convert(cs_.supporting_circle().bbox());
  }
  else
  {
    double x_min = to_double(cs_.source().x());
    double y_min = to_double(cs_.source().y()); 
    double x_max = to_double(cs_.target().x());   
    double y_max = to_double(cs_.target().y());
        
    if(x_min > x_max)
    {
      std::swap(x_min, x_max);
      std::swap(y_min, y_max);
    }

    if(y_min > y_max)
      std::swap(y_min, y_max);
        
    bounding_rect = QRectF(x_min,y_min,x_max,y_max) ;
   }
}


template <typename CS>
void 
GpsCircleSegmentGraphicsItem<CS>::modelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GPS_CIRCLE_SEGMENT_GRAPHICS_ITEM_H
