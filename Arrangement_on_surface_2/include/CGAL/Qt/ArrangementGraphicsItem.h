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

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
//#include <CGAL/apply_to_range.h>
#include <CGAL/Kernel/global_functions.h> // TODO: should be included in PainterOstream.h
#include <CGAL/Qt/ArrangementPainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename TArr >
class ArrangementGraphicsItem : public GraphicsItem
{
  typedef typename TArr::Geometry_traits_2 Traits;
  typedef typename TArr::Vertex_iterator Vertex_iterator;
  typedef typename TArr::Edge_iterator Edge_iterator;
  typedef typename Traits::Kernel Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;

public:
  ArrangementGraphicsItem(TArr* t_);

public slots:
  void modelChanged();

public:

  QRectF boundingRect() const;
  
  virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  const QPen& getVerticesPen() const
  {
    return this->verticesPen;
  }

  const QPen& getEdgesPen() const
  {
    return this->edgesPen;
  }

  void setVerticesPen(const QPen& pen)
  {
    this->verticesPen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    this->edgesPen = pen;
  }

  bool visibleVertices() const
  {
    return this->visible_vertices;
  }

  void setVisibleVertices(const bool b)
  {
    this->visible_vertices = b;
    this->update();
  }

  bool visibleEdges() const
  {
    return this->visible_edges;
  }

  void setVisibleEdges(const bool b)
  {
    this->visible_edges = b;
    this->update();
  }

protected:
  void updateBoundingBox();

  TArr* arr;
  QPainter* m_painter;
  ArrangementPainterOstream< Kernel > painterostream;

  //typename Traits::Point_2 p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;

  QPen verticesPen;
  QPen edgesPen;
  bool visible_edges;
  bool visible_vertices;
  CGAL::Qt::Converter< Traits > convert;
};


template < typename TArr >
ArrangementGraphicsItem< TArr >::ArrangementGraphicsItem( TArr * arr_ )
  :  arr( arr_ ), painterostream( 0 ),
     bb( 0, 0, 0, 0 ), bb_initialized( false ),
     visible_edges( true ), visible_vertices( true )
{
  this->setVerticesPen( QPen( ::Qt::black, 3. ) );
  this->setEdgesPen( QPen( ::Qt::black, 1. ) );
  if ( this->arr->number_of_vertices() == 0 ) {
    this->hide( );
  }
  //this->updateBoundingBox( );
  this->setZValue( 3 );
}

template <typename TArr>
QRectF 
ArrangementGraphicsItem< TArr >::boundingRect() const
{
    QRectF rect = this->convert( this->bb );
    return rect;
}

template <typename TArr>
void 
ArrangementGraphicsItem< TArr >::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
    std::cout << "ArrangementGraphicsItem::paint stub" << std::endl;
    painter->drawRect( this->boundingRect( ) );

    painter->setPen( this->verticesPen );
    this->painterostream = ArrangementPainterOstream< Kernel >( painter, this->boundingRect( ) );
    for ( Vertex_iterator it = this->arr->vertices_begin( ); it != this->arr->vertices_end( ); ++it )
    {
        this->painterostream << it->point( );
    }
    painter->setPen( this->edgesPen );
    for ( Edge_iterator it = this->arr->edges_begin( ); it != this->arr->edges_end( ); ++it )
    {
        Point_2 p1 = it->source( )->point( );
        Point_2 p2 = it->target( )->point( );
        Segment_2 edge( p1, p2 );
        this->painterostream << edge;
    }
#if 0
  painter->setPen(this->edgesPen());
//   painter->drawRect(boundingRect());
  if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
    drawAll(painter);
  } else {
    m_painter = painter;
    painterostream = PainterOstream<Traits>(painter);
    CGAL::apply_to_range (*t, 
                          typename T::Point(option->exposedRect.left(),
                                            option->exposedRect.bottom()), 
                          typename T::Point(option->exposedRect.right(),
                                            option->exposedRect.top()), 
                          *this);
  }
#endif
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename TArr>
void 
ArrangementGraphicsItem< TArr >::updateBoundingBox()
{
    this->prepareGeometryChange( );
    if ( this->arr->number_of_vertices( ) == 0 )
    {
        this->bb = Bbox_2( 0, 0, 0, 0 );
        this->bb_initialized = false;
        return;
    }
    else
    {
        this->bb = this->arr->vertices_begin( )->point( ).bbox( );
        this->bb_initialized = true;
    }

    for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( );
        ++it )
    {
        this->bb = this->bb + it->point( ).bbox( );
    }
#if 0
  prepareGeometryChange();
  if(t->number_of_vertices() == 0){
    bb = Bbox_2(0,0,0,0);
    bb_initialized = false;
    return;
  } else if(! bb_initialized){
    bb = t->finite_vertices_begin()->point().bbox();
    bb_initialized = true;
  }
  
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
#endif
}


template <typename TArr>
void 
ArrangementGraphicsItem< TArr >::modelChanged()
{
  std::cout << "ArrangementGraphicsItem modelChanged stub" << std::endl;
  if ( this->arr->is_empty( ) )
  {
      this->hide( );
  }
  else
  {
      this->show( );
  }
  this->updateBoundingBox();
  this->update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
