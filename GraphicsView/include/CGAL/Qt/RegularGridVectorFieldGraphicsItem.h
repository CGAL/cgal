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

#ifndef CGAL_QT_REGULAR_GRID_VECTOR_FIELD_GRAPHICS_ITEM_H
#define CGAL_QT_REGULAR_GRID_VECTOR_FIELD_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

  template <typename T, typename K>
class RegularGridVectorFieldGraphicsItem : public GraphicsItem
{
  typedef typename T::Geom_traits Geom_traits;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Segment_2 Segment_2;

public:
  RegularGridVectorFieldGraphicsItem(T* t_);

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

protected:
  void updateBoundingBox();

  T * rg;
  QPainter* m_painter;
  PainterOstream<Geom_traits> painterostream;

  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
};


  template <typename T, typename K>
  RegularGridVectorFieldGraphicsItem<T,K>::RegularGridVectorFieldGraphicsItem(T * t_)
  :  rg(t_), painterostream(0),
     visible_edges(true), visible_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  updateBoundingBox();
  setZValue(3);
}

  template <typename T, typename K>
QRectF 
  RegularGridVectorFieldGraphicsItem<T,K>::boundingRect() const
{
  return bounding_rect;
}






  template <typename T, typename K>
void 
  RegularGridVectorFieldGraphicsItem<T,K>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{

  painterostream = PainterOstream<Geom_traits>(painter);
  painter->setPen(this->edgesPen());
  double w  = rg->get_size().first;
  double h  = rg->get_size().second;
  int nw =  rg->get_dimension().first;
  int nh =  rg->get_dimension().second;
  double dw = w/(nw-1);
  double dh = h/(nh-1);

  for(int i = 0; i < nw; i++){
    for(int j = 0; j < nh; j++){
      Vector_2 v = rg->get_field(i,j);
      v = (dw*0.45) * v/sqrt(v*v);
      painterostream << Segment_2(Point_2(i*dw,j*dh),
                                  Point_2(i*dw,j*dw)+v);
    }
  }
  painter->setPen(this->verticesPen());
  QMatrix matrix = painter->matrix();
  painter->resetMatrix();
  for(int i = 0; i < nw; i++){
    for(int j = 0; j < nh; j++){
      painter->drawPoint(matrix.map(QPointF(i*dw, j*dh)));
    }
  }

}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
  template <typename T, typename K>
void 
  RegularGridVectorFieldGraphicsItem<T,K>::updateBoundingBox()
{

  bounding_rect = QRectF(0,
                         0,
                         rg->get_size().first,
                         rg->get_size().second);
}


  template <typename T, typename K>
void 
  RegularGridVectorFieldGraphicsItem<T,K>::modelChanged()
{
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_REGULAR_GRID_VECTOR_FIELD_GRAPHICS_ITEM_H
