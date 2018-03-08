// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/GraphicsView/include/CGAL/Qt/TriangulationGraphicsItem.h $
// $Id: TriangulationGraphicsItem.h 61414 2011-02-24 16:36:04Z sloriot $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/apply_to_range.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>


//#define DRAW_OCTAGON_IMAGE ;

namespace CGAL {
namespace Qt {

template <typename T>
class TriangulationGraphicsItem : public GraphicsItem
{
  typedef typename T::Geometric_traits  Geom_traits;
  typedef typename Geom_traits::Point_2 Point_2;
public:
  TriangulationGraphicsItem(T* t_);

  void modelChanged();

public:

  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  virtual void operator()(typename T::Face_handle fh);

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

  void setVisibleOctagon(const bool b)
  {
    visible_octagon = b;
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

  void setMovingPoint(Point_2 p) {
    moving_point = p;
  }

  void setSource(Point_2 p) {
    source = p;
  }

  void setTarget(Point_2 p) {
    target = p;
  }


protected:
  virtual void drawAll(QPainter *painter);
  void paintVertices(QPainter *painter);
  void paintOneVertex(const typename T::Point& point);
  virtual void paintVertex(typename T::Vertex_handle vh);
  void updateBoundingBox();

  T * t;
  QPainter* m_painter;
  PainterOstream<Geom_traits> painterostream;

  typename T::Vertex_handle vh;
  typename T::Point p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
  bool visible_octagon;

  Point_2 moving_point;
  Point_2 source;
  Point_2 target;
};


template <typename T>
TriangulationGraphicsItem<T>::TriangulationGraphicsItem(T * t_)
  :  t(t_), painterostream(0),
     bb(0,0,0,0), bb_initialized(false),
     visible_edges(true), visible_vertices(true), visible_octagon(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(t->number_of_vertices() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}

template <typename T>
QRectF 
TriangulationGraphicsItem<T>::boundingRect() const
{
  return bounding_rect;
}


template <typename T>
void 
TriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if(visible_edges) {
    for (int i=0; i<3; i++) {
      if (fh < fh->neighbor(i)) {
        m_painter->setPen(this->edgesPen());
        painterostream << t->segment(fh, i);
      }
    }
  }
  if(visible_vertices) {
    for (int i=0; i<3; i++) {
      paintVertex(fh->vertex(i));
    }
  }
}

template <typename T>
void 
TriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  //delete
  QPen temp = painter->pen();
  QPen old = temp;
  
  temp.setWidthF(0.0125);
  temp.setColor(::Qt::black);
  painter->setPen(temp);
  painterostream = PainterOstream<Geom_traits>(painter);
  
  painter->drawEllipse(QRectF(-1.0, -1.0, 2.0, 2.0));

  if (visible_octagon) {
    typedef typename Geom_traits::FT        FT;
    typedef typename Geom_traits::Point_2   Point_2;
    FT ep = CGAL::sqrt(FT(2)+CGAL::sqrt(FT(2)));
    FT em = CGAL::sqrt(FT(2)-CGAL::sqrt(FT(2)));
    FT p14( CGAL::sqrt(CGAL::sqrt(FT(2))) );
    FT p34( p14*p14*p14 );

    Point_2 v0(ep*p34/FT(4), -em*p34/FT(4));
    Point_2 v1(p14*(ep + em)/FT(4), p14*(ep - em)/FT(4));
    Point_2 v2(em*p34/FT(4), ep*p34/FT(4));
    Point_2 v3(-p14*(ep - em)/FT(4), p14*(ep + em)/FT(4));
    Point_2 v4(-ep*p34/FT(4), em*p34/FT(4));
    Point_2 v5(-p14*(ep + em)/FT(4), -p14*(ep - em)/FT(4));
    Point_2 v6(-em*p34/FT(4), -ep*p34/FT(4));
    Point_2 v7(p14*(ep - em)/FT(4), -p14*(ep + em)/FT(4));

    painterostream << t->segment(v0, v1);
    painterostream << t->segment(v1, v2);
    painterostream << t->segment(v2, v3);
    painterostream << t->segment(v3, v4);
    painterostream << t->segment(v4, v5);
    painterostream << t->segment(v5, v6);
    painterostream << t->segment(v6, v7);
    painterostream << t->segment(v7, v0);
  }


  temp.setWidthF(0.01);
  temp.setColor(::Qt::darkGreen);
  painter->setPen(temp);
  painterostream = PainterOstream<Geom_traits>(painter);

  // for (typename T::Face_iterator fit = t->faces_begin();
  //      fit != t->faces_end(); fit++) {

  //   for (int k = 0; k < 3; k++) {
  //     painterostream << t->segment(fit, k);
  //   }

  // }

  paintVertices(painter);

}


template <typename T>
void 
TriangulationGraphicsItem<T>::paintVertices(QPainter *painter)
{

    Converter<Geom_traits> convert;

    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
      
    // delete
    QPen temp = painter->pen();
    QPen old = temp;
    temp.setWidth(9);
    
    double px = to_double(moving_point.x());
    double py = to_double(moving_point.y());
    double dist = px*px + py*py;

    double width = 10.*(exp(.85) - exp(dist));
    temp.setWidthF(width);
    temp.setColor(::Qt::blue);
    painter->setPen(temp);
    
    QPointF point = matrix.map(convert(moving_point));
    painter->drawPoint(point);

    //--------------------------------------

    double sx = to_double(source.x());
    double sy = to_double(source.y());
    double sdist = sx*sx + sy*sy;

    width = 10;
    temp.setWidthF(width);
    temp.setColor(::Qt::red);
    painter->setPen(temp);

    QPointF spoint = matrix.map(convert(source));
    painter->drawPoint(spoint);

    //--------------------------------------

    double tx = to_double(target.x());
    double ty = to_double(target.y());
    double tdist = tx*tx + ty*ty;

    QPointF tpoint = matrix.map(convert(target));
    painter->drawPoint(tpoint);
    
    painter->setPen(old);
  
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(point)));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;
  m_painter->setPen(this->verticesPen());
  
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(vh->point())));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
//   painter->drawRect(boundingRect());
  if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
    drawAll(painter);
  } else {
    m_painter = painter;
    painterostream = PainterOstream<Geom_traits>(painter);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
TriangulationGraphicsItem<T>::updateBoundingBox()
{
  prepareGeometryChange();
  bounding_rect = QRectF(-1., -1., 2., 2.);
}


template <typename T>
void 
TriangulationGraphicsItem<T>::modelChanged()
{
  if((t->number_of_vertices() == 0) ){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
