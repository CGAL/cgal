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
// Author(s)     : Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_DELAUNAY_MESH_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_DELAUNAY_MESH_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Qt/ConstrainedTriangulationGraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <QBrush>
#include <QPen>
#include <list>

namespace CGAL {

namespace Qt {

template <typename T>
class DelaunayMeshTriangulationGraphicsItem
  : public ConstrainedTriangulationGraphicsItem<T>
{
  typedef ConstrainedTriangulationGraphicsItem<T> Base;
  typedef typename T::Geom_traits Geom_traits;
  typedef typename T::Point Point;

public:
  DelaunayMeshTriangulationGraphicsItem(T  * t_)
    : Base(t_)
    , visible_in_domain(true)
    , visible_blind_faces(false)
    , visible_voronoi(false)
    , visible_seeds(false)
    , visible_inside_edges(false)
    , seeds_begin()
    , seeds_end()
    , in_domain_brush(::Qt::blue)
    , blind_brush(::Qt::gray)
  {
    setSeedsPen(QPen(::Qt::black, 10.));
  }
  
  void operator()(typename T::Face_handle fh);

  const QBrush& facesInDomainBrush() const
  {
    return in_domain_brush;
  }

  void setFacesInDomainBrush(const QBrush& brush)
  {
    in_domain_brush = brush;
  }

  const QBrush& blindFacesBrush() const
  {
    return blind_brush;
  }

  void setBlindFacesBrush(const QBrush& brush)
  {
    blind_brush = brush;
  }

  bool visibleFacesInDomain() const
  {
    return visible_in_domain;
  }

  void setVisibleFacesInDomain(const bool b)
  {
    visible_in_domain = b;
    this->update();
  }

  bool visibleBlindFaces() const
  {
    return visible_blind_faces;
  }

  void setVisibleBlindFaces(const bool b)
  {
    visible_blind_faces = b;
    this->update();
  }

  bool visibleSeeds() const
  {
    return visible_seeds;
  }

  void setVisibleSeeds(const bool b,
    typename std::list<Point>::iterator begin,
    typename std::list<Point>::iterator end)
  {
    visible_seeds = b;
    seeds_begin = begin;
    seeds_end = end;
    this->update();
  }

  const QPen& voronoiPen() const
  {
    return voronoi_pen;
  }

  void setVoronoiPen(const QPen& pen)
  {
    voronoi_pen = pen;
  }

  const QPen& seedsPen() const
  {
    return seeds_pen;
  }

  void setSeedsPen(const QPen& pen)
  {
    seeds_pen = pen;
  }

  bool visibleVoronoiEdges() const
  {
    return visible_voronoi;
  }

  void setVisibleVoronoiEdges(const bool b)
  {
    visible_voronoi = b;
    this->update();
  }

  bool visibleInsideEdges() const
  {
    return visible_inside_edges;
  }

  void setVisibleInsideEdges(const bool b)
  {
    visible_inside_edges = b;
    this->update();
  }


protected:
  void drawAll(QPainter *painter);
  void paintSeeds(QPainter *painter);

  bool visible_in_domain;
  bool visible_blind_faces;
  bool visible_voronoi;
  bool visible_seeds;
  bool visible_inside_edges;
  typename std::list<Point>::iterator seeds_begin, seeds_end;
  
  QBrush in_domain_brush;
  QBrush blind_brush;
  QPen voronoi_pen;
  QPen seeds_pen;
};

template <typename T>
void 
DelaunayMeshTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  if(visibleFacesInDomain())
  {
    this->painterostream = PainterOstream<typename T::Geom_traits>(painter);
    painter->setBrush(facesInDomainBrush());
    painter->setPen(::Qt::NoPen);
    for(typename T::Finite_faces_iterator fit = this->t->finite_faces_begin();
	fit != this->t->finite_faces_end();
	++fit){
      if(fit->is_in_domain()){
	this->painterostream << this->t->triangle(fit);
      }
    }
    painter->setBrush(::Qt::NoBrush);
  }
  if(visibleInsideEdges())
  {
    this->painterostream = PainterOstream<typename T::Geom_traits>(painter);
    painter->setBrush(::Qt::NoBrush);
    painter->setPen(this->edgesPen());
    for(typename T::Finite_faces_iterator fit = this->t->finite_faces_begin();
        fit != this->t->finite_faces_end();
        ++fit){
      if(fit->is_in_domain()){
        this->painterostream << this->t->triangle(fit);
      }
    }
    painter->setBrush(::Qt::NoBrush);
  }
  if(visibleBlindFaces())
  {
    this->painterostream = PainterOstream<typename T::Geom_traits>(painter);
    painter->setBrush(blindFacesBrush());
    painter->setPen(::Qt::NoPen);
    for(typename T::Finite_faces_iterator fit = this->t->finite_faces_begin();
	fit != this->t->finite_faces_end();
	++fit){
      if(fit->is_blind()){
	this->painterostream << this->t->triangle(fit);
      }
    }
    painter->setBrush(::Qt::NoBrush);
  }
  if(this->visibleVoronoiEdges())
  {
    painter->setBrush(::Qt::NoBrush);
    painter->setPen(::Qt::darkGreen);
    //painter->setPen(this->voronoiPen());

    this->painterostream = PainterOstream<typename T::Geom_traits>(painter);
    typedef CGAL::Cvd_cell_2<T> Cvd_cell;
    for(typename T::Finite_vertices_iterator
          vit = this->t->finite_vertices_begin();
          vit != this->t->finite_vertices_end();
          ++vit)
    {
      Cvd_cell cell = CGAL::dual(*(this->t), vit);
      for(typename Cvd_cell::segment_iterator sit = cell.segments_begin();
        sit != cell.segments_end();
        ++sit)
      {
        this->painterostream << *sit;
      }
      for(typename Cvd_cell::ray_iterator rit = cell.rays_begin();
        rit != cell.rays_end();
        ++rit)
      {
        this->painterostream << *rit;
      }
    }
  }
  paintSeeds(painter);

  Base::drawAll(painter);
}

template <typename T>
void
DelaunayMeshTriangulationGraphicsItem<T>::paintSeeds(QPainter *painter)
{
  if(this->visibleSeeds())
  {
    Converter<Geom_traits> convert;
    painter->setPen(this->seedsPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();

    typename std::list<Point>::iterator sit;
    for(sit = seeds_begin; sit != seeds_end; ++sit)
      painter->drawPoint(matrix.map(convert(*sit)));

    painter->setMatrix(matrix);
  }
}

template <typename T>
void 
DelaunayMeshTriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if(visibleFacesInDomain()) {
    if(fh->is_in_domain()){
      this->painterostream = PainterOstream<typename T::Geom_traits>(this->m_painter);
      this->m_painter->setBrush(facesInDomainBrush());
      this->m_painter->setPen(::Qt::NoPen) ;
      this->painterostream << this->t->triangle(fh);
    }
  }
  if(visibleBlindFaces()) {
    if(fh->is_blind()){
      this->painterostream = PainterOstream<typename T::Geom_traits>(this->m_painter);
      this->m_painter->setBrush(blindFacesBrush());
      this->m_painter->setPen(::Qt::NoPen) ;
      this->painterostream << this->t->triangle(fh);
    }
  }
  Base::operator()(fh);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_DELAUNAY_MESH_TRIANGULATION_GRAPHICS_ITEM_H
