// Copyright (c) 2008  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/Hyperbolic_octagon_translation.h>

//#define DRAW_OCTAGON_IMAGE ;

namespace CGAL {
namespace Qt {

template <typename T>
class TriangulationGraphicsItem : public GraphicsItem
{
  typedef typename T::Geom_traits Geom_traits;
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

  void setVisibleConflictZone(const bool b) {
    visible_cz = b;
  }

  void setVisibleDemo(const bool val) {
    visible_demo = val;
  }

  void setVisibleCopies(const bool val) {
    if(val)
      initializeTranslations();
    visible_copies = val;
  }

  bool visibleCopies() {
    return visible_copies;
  }

  bool visibleConflictZone() const {
    return visible_cz;
  }

  void setSource(typename T::Point src) {
    source = src;
  }


  void setTarget(typename T::Point tgt) {
    target = tgt;
  }

  void setMovingPoint(typename T::Point mp) {
    this->moving_point = mp;
    cfaces.clear();
    t->find_conflicts(moving_point, std::back_inserter(cfaces));
  }

  void initializeTranslations() {
    std::vector<Matrix> gens;
    T::Hyperbolic_translation::generators(gens);
    std::vector<Matrix> tmp;
    for(unsigned int j=0; j<gens.size(); ++j) {
      tmp.push_back(gens[j]);
    }

    for(int j=0; j<3; ++j) {
      std::size_t N = tmp.size();
      for(unsigned int i=0; i<N; ++i) {
        for(unsigned int k=0; k< gens.size(); k++) {
          tmp.push_back(gens[k]*tmp[i]);
        }
      }
    }

    for(unsigned int k=0; k< tmp.size(); k++) {
      if(std::find(trans.begin(), trans.end(), tmp[k]) == trans.end()) {
        trans.push_back(tmp[k]);
      }
    }
    std::cout << "translations ready; tmp has " << tmp.size() << " elems, trans has " << trans.size() << " elems" << std::endl;
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
  typename T::Point moving_point;
  typename T::Point source, target;
  CGAL::Bbox_2 bb;
  bool bb_initialized;
  QRectF bounding_rect;
  Converter<Geom_traits> convert;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
  bool visible_octagon;
  bool visible_cz;
  bool visible_demo;
  bool visible_copies;

  std::list<typename T::Face_handle> cfaces;

  typedef typename T::Hyperbolic_translation Matrix;
  std::vector<Matrix> trans;

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
    for(int i=0; i<3; ++i) {
      if(fh < fh->neighbor(i)) {
        m_painter->setPen(this->edgesPen());
        painterostream << t->construct_hyperbolic_segment(fh, i);
      }
    }
  }
  if(visible_vertices) {
    for(int i=0; i<3; ++i) {
      paintVertex(fh->vertex(i));
    }
  }
}

template <typename T>
void
TriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  QPen temp = painter->pen();
  QPen old = temp;
  temp.setWidthF(0.0125);
  temp.setColor(::Qt::black);

  QColor bcol = ::Qt::green;
  bcol.setAlpha(120);
  QBrush brush;
  brush.setColor(bcol);

  painter->setPen(temp);
  painter->setBrush(brush);

  painterostream = PainterOstream<Geom_traits>(painter);

  painter->drawEllipse(QRectF(-1.0, -1.0, 2.0, 2.0));

  if(visible_octagon) {
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

    painterostream << t->construct_hyperbolic_segment(v0, v1);
    painterostream << t->construct_hyperbolic_segment(v1, v2);
    painterostream << t->construct_hyperbolic_segment(v2, v3);
    painterostream << t->construct_hyperbolic_segment(v3, v4);
    painterostream << t->construct_hyperbolic_segment(v4, v5);
    painterostream << t->construct_hyperbolic_segment(v5, v6);
    painterostream << t->construct_hyperbolic_segment(v6, v7);
    painterostream << t->construct_hyperbolic_segment(v7, v0);
  }


  if(visible_cz) {
    painter->setBrush(QColor(122, 20, 39, 100));
    QPen oldpen = painter->pen();
    QPen npen(QColor(122, 20, 39, 0));
    npen.setWidthF(0.0);
    painter->setPen(npen);
    //painter->pen().setColor(::Qt::white);
    painterostream = PainterOstream<Geom_traits>(painter);
    for(typename std::list<typename T::Face_handle>::iterator it = cfaces.begin(); it != cfaces.end(); it++) {
      painterostream << t->construct_hyperbolic_triangle(*it);
    }
    painter->setPen(oldpen);
  }

  if(visible_edges) {
    //cout << "painting edges" << std::endl;
    temp.setWidthF(0.01);
    temp.setColor(::Qt::darkGreen);
    painter->setPen(temp);
    painterostream = PainterOstream<Geom_traits>(painter);

    typedef typename Geom_traits::Construct_hyperbolic_point_2 CP2;

    for(typename T::Face_iterator fit = t->faces_begin(); fit != t->faces_end(); fit++) {
      typename Geom_traits::Point_2 pts[] = { CP2()(fit->vertex(0)->point(), fit->translation(0)),
                                              CP2()(fit->vertex(1)->point(), fit->translation(1)),
                                              CP2()(fit->vertex(2)->point(), fit->translation(2)) } ;

      painterostream << t->construct_hyperbolic_segment(pts[0], pts[1]);
      painterostream << t->construct_hyperbolic_segment(pts[1], pts[2]);
      painterostream << t->construct_hyperbolic_segment(pts[2], pts[0]);
      //cout << " original painted" << std::endl;
      if(visible_copies) {
      //   std::cout << "painting copies" << std::endl;
      //   typename Geom_traits::Point_2 pts[] = {fit->translation(Triangulation_cw_ccw_2::ccw(0)).apply(fit->vertex(Triangulation_cw_ccw_2::ccw(k))->point());
      //   typename Geom_traits::Point_2 tgt = fit->translation(Triangulation_cw_ccw_2::cw(k)).apply(fit->vertex(Triangulation_cw_ccw_2::cw(k))->point());
        for(unsigned int j=0; j<trans.size(); ++j) {
          painterostream << t->construct_hyperbolic_segment( CP2()(pts[0], trans[j]), CP2()(pts[1], trans[j]) );
          painterostream << t->construct_hyperbolic_segment( CP2()(pts[1], trans[j]), CP2()(pts[2], trans[j]) );
          painterostream << t->construct_hyperbolic_segment( CP2()(pts[2], trans[j]), CP2()(pts[0], trans[j]) );
        }
      }
      //cout << "   copies painted" << std::endl;
    }

    if(visible_demo) {
      temp.setColor(::Qt::red);
      painter->setPen(temp);
      painterostream = PainterOstream<Geom_traits>(painter);
      painterostream << t->construct_hyperbolic_segment(source, target);
    }
  }

  //delete
  painter->setPen(old);
  //

  //cout << "painting vertices" << std::endl;
  paintVertices(painter);
}

template <typename T>
void
TriangulationGraphicsItem<T>::paintVertices(QPainter *painter)
{
  if(visibleVertices()) {
    Converter<Geom_traits> convert;

    painter->setPen(verticesPen());
    QTransform matrix = painter->worldTransform();
    //QMatrix tr(1, 0, 0, -1, 0, 0);
    //matrix = tr*matrix;
    painter->resetTransform();

    QPen temp = painter->pen();
    QPen old = temp;

    double width;

    for(typename T::Vertex_iterator it = t->vertices_begin();
        it != t->vertices_end();
        it++){

      double px = to_double(it->point().x());
      double py = to_double(it->point().y());
      double dist = px*px + py*py;

      width = 20.*(exp(.85) - exp(dist));
      temp.setWidthF(width);

      if(t->is_dummy_vertex(it)) {
        temp.setColor(::Qt::green);
      } else {
        temp.setColor(::Qt::red);
      }

      painter->setPen(temp);
      QPointF point = matrix.map(convert(it->point()));
      painter->drawPoint(point);

      if(visible_copies) {
        for(unsigned int k=0; k< trans.size(); k++) {
          typedef typename Geom_traits::Construct_hyperbolic_point_2 CP2;
          typename Geom_traits::Point_2 img = CP2()(it->point(), trans[k]);

          px = to_double(img.x());
          py = to_double(img.y());
          dist = px*px + py*py;

          width = 20.*(exp(1.0) - exp(dist));
          temp.setWidthF(width);
          painter->setPen(temp);
          painter->drawPoint(matrix.map(convert(img)));
        }
      }

    }

    if(visible_demo) {

      width = 10;
      temp.setWidthF(width);
      temp.setColor(::Qt::red);
      painter->setPen(temp);

      QPointF spoint = matrix.map(convert(source));
      painter->drawPoint(spoint);

      //--------------------------------------

      QPointF tpoint = matrix.map(convert(target));
      painter->drawPoint(tpoint);

      //--------------------------------------

      double px = to_double(moving_point.x());
      double py = to_double(moving_point.y());
      double dist = px*px + py*py;

      double width = 10.*(exp(.85) - exp(dist));
      temp.setWidthF(width);
      temp.setColor(::Qt::blue);
      painter->setPen(temp);

      QPointF point = matrix.map(convert(moving_point));
      painter->drawPoint(point);
    }

    painter->setPen(old);
  }
}

template <typename T>
void
TriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QTransform matrix = m_painter->worldTransform();
  m_painter->resetTransform();
  m_painter->drawPoint(matrix.map(convert(point)));
  m_painter->setWorldTransform(matrix);
}

template <typename T>
void
TriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;
  m_painter->setPen(this->verticesPen());

  QTransform matrix = m_painter->worldTransform();
  m_painter->resetTransform();
  m_painter->drawPoint(matrix.map(convert(vh->point())));
  m_painter->setWorldTransform(matrix);
}

template <typename T>
void
TriangulationGraphicsItem<T>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
  if( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
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
