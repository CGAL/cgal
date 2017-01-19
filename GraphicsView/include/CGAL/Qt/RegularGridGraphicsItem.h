// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_REGULAR_GRID_GRAPHICS_ITEM_H
#define CGAL_QT_REGULAR_GRID_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

  template <typename K>
class RegularGridGraphicsItem : public GraphicsItem
{
  typedef K Geom_traits;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Segment_2 Segment_2;

public:
  RegularGridGraphicsItem(double dx, double dy);

  void modelChanged();

public:

  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

  bool visibleVertices() const
  {
    return visible_vertices;
  }

  void setVisibleVertices(const bool b)
  {
    visible_vertices = b;
    update();
  }

  bool visibleEdges() const
  {
    return visible_edges;
  }

  void setVisibleEdges(const bool b)
  {
    visible_edges = b;
    update();
  }

  void setDelta(double x, double y)
  {
    dx = x;
    dy = y;
    update();
  }

protected:
  void updateBoundingBox();
  double dx, dy;
  QPainter* m_painter;
  PainterOstream<Geom_traits> painterostream;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
};


  template <typename K>
  RegularGridGraphicsItem<K>::RegularGridGraphicsItem(double dx, double dy)
    :  dx(dx), dy(dy), painterostream(0),
     visible_edges(true), visible_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  setZValue(-1);
}

  template <typename K>
QRectF 
  RegularGridGraphicsItem<K>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}






  template <typename K>
void 
  RegularGridGraphicsItem<K>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  QRectF rect = boundingRect();
  double b = rect.bottom();
  double t = rect.top();
  double l = rect.left();
  double r = rect.right();

  if(b > t) std::swap(b,t); // because things are upside down in Qt
  
  painterostream = PainterOstream<Geom_traits>(painter);
  painter->setPen(this->edgesPen());

  double ll = l;
  ll = dx * static_cast<int>(ll/dx);
  
  for(; ll < r; ll += dx){
    painterostream << Segment_2(Point_2(ll,b),
                                Point_2(ll,t));
  }


  double bb = b;
  bb = dy * static_cast<int>(bb/dy);
  
  for(; bb < t; bb += dy){
    painterostream << Segment_2(Point_2(l,bb),
                                Point_2(r,bb));
  }


  /*
  painter->setPen(this->verticesPen());
  QMatrix matrix = painter->matrix();
  painter->resetMatrix();
  for(int i = 0; i < nw; i++){
    for(int j = 0; j < nh; j++){
      painter->drawPoint(matrix.map(QPointF(i*dw, j*dh)));
    }
  }
  */
}

  template <typename K>
void 
  RegularGridGraphicsItem<K>::updateBoundingBox()
{}


  template <typename K>
void 
  RegularGridGraphicsItem<K>::modelChanged()
{
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_REGULAR_GRID_GRAPHICS_ITEM_H
