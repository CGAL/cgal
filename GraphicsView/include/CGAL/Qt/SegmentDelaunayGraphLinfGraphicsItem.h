// Copyright (c) 2015  Università della Svizzera italiana.
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_QT_SEGMENT_DELAUNAY_GRAPH_LINF_GRAPHICS_ITEM_H
#define CGAL_QT_SEGMENT_DELAUNAY_GRAPH_LINF_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/utility.h>
//#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/Polychain_2.h>


namespace CGAL {

namespace Qt {

template <typename T>
class SegmentDelaunayGraphLinfGraphicsItem : public GraphicsItem
{ 
  typedef typename T::Geom_traits Geom_traits;
  typedef typename T::Point_2 Point_2;
  typedef typename Kernel_traits<Point_2> ::Kernel Kern;

  T* t;
  QPainter* m_painter;
  PainterOstream<Kern> painterostream;

  QPen vertices_pen, segment_pen, voronoi_pen ;
  Bbox_2 bb;
  bool bb_initialized;

public:
  SegmentDelaunayGraphLinfGraphicsItem(T  * t_)
    : t(t_), painterostream(0),
      segment_pen(::Qt::blue),
      voronoi_pen(::Qt::blue)
  {
  }
  

  void updateBoundingBox();

  void modelChanged();


  QRectF boundingRect() const;
 
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
 
  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }


  const QPen& segmentPen() const
  {
    return segment_pen;
  }

  void setSegmentPen(const QPen& pen)
  {
    segment_pen = pen;
  }



  const QPen& voronoiPen() const
  {
    return voronoi_pen;
  }

  void setVoronoiPen(const QPen& pen)
  {
    voronoi_pen = pen;
  }


protected:
  void drawAll(QPainter *painter, const QStyleOptionGraphicsItem *option);
  void drawDualEdge(QPainter *painter, typename T::Edge e);



};

template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::drawDualEdge(QPainter * /*painter*/, typename T::Edge e)
{
  CGAL_precondition( ! t->is_infinite(e) );

  CGAL::Polychainline_2<Geom_traits>      pcl;
  CGAL::Polychainray_2<Geom_traits>       pcr;
  CGAL::Polychainsegment_2<Geom_traits>   pcs;
//  CGAL::Polychainline_2<Kern>         pcl;
//  CGAL::Polychainray_2<Kern>          pcr;
//  CGAL::Polychainsegment_2<Kern>      pcs;

//  std::cout << "debug drawDE pcl type= " 
//            << typeid(pcl).name() << std::endl;
//  std::cout << "debug drawDE pcr type= " 
//            << typeid(pcr).name() << std::endl;
//  std::cout << "debug drawDE pcs type= " 
//            << typeid(pcs).name() << std::endl;

  Object o = t->primal(e);

  //std::cout << "debug drawDE trying to assign object at "
  //          << &o << std::endl;

  if (CGAL::assign(pcl, o)) { 
    /* m_painter->setPen(::Qt::cyan); 
       std::cerr << "line " << std::endl; */ 
//    std::cout << "debug drawDE pcl-assigned object at "
//              << &o << std::endl;
    //painterostream << pcl; 
    pcl.draw(painterostream);      
  }
  else if (CGAL::assign(pcr, o)) { 
    /* m_painter->setPen(::Qt::magenta); 
       std::cerr << "segment " << std::endl; */ 
//    std::cout << "debug drawDE pcr-assigned object at "
//              << &o << std::endl;
    //painterostream << pcr;
    pcr.draw(painterostream);      
  }
  else if (CGAL::assign(pcs, o))  { 
    /* m_painter->setPen(::Qt::darkMagenta);  
       std::cerr << "ray " << r << std::endl;  */ 
//    std::cout << "debug drawDE pcs-assigned object at "
//              << &o << std::endl;
    //painterostream << pcs; 
    pcs.draw(painterostream);      
  }
  else { 
    std::cerr << "error: drawDualEdge unknown curve" << std::endl; 
  }

  /* m_painter->setPen(::Qt::black); */
}


template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::drawAll(QPainter *painter, const QStyleOptionGraphicsItem *option)
{
  QRectF rect = option->exposedRect;
  m_painter = painter;
  painterostream = PainterOstream<Kern>(m_painter, rect);
  m_painter->setPen(this->voronoiPen());
  typename T::Finite_edges_iterator eit = t->finite_edges_begin();
  for (; eit != t->finite_edges_end(); ++eit) {
    drawDualEdge(m_painter, *eit);
  }
  {
    m_painter->setPen(this->segmentPen());
      typename T::Finite_vertices_iterator vit;
      for (vit = t->finite_vertices_begin();
	   vit != t->finite_vertices_end(); ++vit) {
	typename T::Site_2 s = vit->site();
	if ( s.is_segment() ) {
	  painterostream << s.segment();
	}
      }
    }
    {
    m_painter->setPen(this->verticesPen());
    QMatrix matrix = m_painter->matrix();
    m_painter->resetMatrix();
    Converter<Kern> convert;
      typename T::Finite_vertices_iterator vit;
      for (vit = t->finite_vertices_begin();
	   vit != t->finite_vertices_end(); ++vit) {
	typename T::Site_2 s = vit->site();
	if ( s.is_input() ) {
	  //*widget << CGAL::RED;
	} else {
	  //*widget << CGAL::YELLOW;
	}
	if ( s.is_point() ) {
          QPointF point = matrix.map(convert(s.point()));
          m_painter->drawPoint(point);
	}
      }
    }

}

/*
template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if(visibleFacesInDomain()) {
    if(fh->is_in_domain()){
      this->painterostream = PainterOstream<typename T::Geom_traits>(this->m_painter);
      this->m_painter->setBrush(facesInDomainBrush());
      this->m_painter->setPen(::Qt::NoPen) ;
      this->painterostream << this->t->triangle(fh);
    }
  }
  Base::operator()(fh);
}
*/


template <typename T>
QRectF 
SegmentDelaunayGraphLinfGraphicsItem<T>::boundingRect() const
{

  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::modelChanged()
{
  if((t->number_of_vertices() == 0) ){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    this->show();
  }
  update();
}

template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::updateBoundingBox()
{
  prepareGeometryChange();
  if(t->number_of_vertices() == 0){
    bb = Bbox_2(0,0,0,0);
    bb_initialized = false;
    return;
  } else if(! bb_initialized){
    //    bb = t->finite_vertices_begin()->point().bbox();
    bb_initialized = true;
  }
  /*
  if(t->dimension() <2){
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
	it != t->finite_vertices_end();
	++it){
      bb = bb + it->point().bbox();
    }
  } else {
    typename T::Vertex_handle inf = t->infinite_vertex();
    typename T::Vertex_circulator vc = t->incident_vertices(inf), done(vc);
    do {
      bb = bb + vc->point().bbox();
      ++vc;
    } while(vc != done);
  }
  bounding_rect = QRectF(bb.xmin(),
                         bb.ymin(),
                         bb.xmax()-bb.xmin(),
                         bb.ymax()-bb.ymin());
  */
}



template <typename T>
void 
SegmentDelaunayGraphLinfGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{

//   painter->drawRect(boundingRect());
  //  if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {

    drawAll(painter, option);
    //  } else {
    //    std::cerr << "else" << std::endl;
    //  }
}



} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_SEGMENT_DELAUNAY_GRAPH_LINF_GRAPHICS_ITEM_H
