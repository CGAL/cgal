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
#include <QBrush>

namespace CGAL {

namespace Qt {

template <typename T>
class DelaunayMeshTriangulationGraphicsItem : public ConstrainedTriangulationGraphicsItem<T>
{
  typedef ConstrainedTriangulationGraphicsItem<T> Base;

public:
  DelaunayMeshTriangulationGraphicsItem(T  * t_)
    : Base(t_),
      visible_in_domain(true),
      in_domain_brush(::Qt::blue)
  {
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

  bool visibleFacesInDomain() const
  {
    return visible_in_domain;
  }

  void setVisibleFacesInDomain(const bool b)
  {
    visible_in_domain = b;
    this->update();
  }

protected:
  void drawAll(QPainter *painter);

  bool visible_in_domain;
  QBrush in_domain_brush;
};

template <typename T>
void 
DelaunayMeshTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  if(visibleFacesInDomain()) {
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
  }
  Base::drawAll(painter);
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
  Base::operator()(fh);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_DELAUNAY_MESH_TRIANGULATION_GRAPHICS_ITEM_H
