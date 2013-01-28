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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_TRIANGULATION_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_TRIANGULATION_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <QPen>

namespace CGAL {

namespace Qt {

template <typename T>
class TriangulationArrangementGraphicsItem : public TriangulationGraphicsItem<T>
{
  typedef typename T::Geom_traits Geom_traits;
public:
  TriangulationArrangementGraphicsItem(T  * t_)
    : TriangulationGraphicsItem<T>(t_), visible_constraints(true)
  {
    constraints_pen = this->edgesPen();
    constraints_pen.setColor(::Qt::red);
    
    fixed_vertices_pen = this->verticesPen();
    fixed_vertices_pen.setColor(::Qt::blue);
  }
  
  void operator()(typename T::Face_handle fh);

  const QPen& constraintsPen() const
  {
    return constraints_pen;
  }

  void setConstraintsPen(const QPen& pen)
  {
    constraints_pen = pen;
  }

  const QPen& fixedVerticesPen() const
  {
    return fixed_vertices_pen;
  }

  void setFixedVerticesPen(const QPen& pen)
  {
    fixed_vertices_pen = pen;
  }

  bool visibleConstraints() const
  {
    return visible_constraints;
  }

  void setVisibleConstraints(const bool b)
  {
    visible_constraints = b;
    this->update();
  }

protected:

  void drawAll(QPainter *painter);
  void paintVertices(QPainter *painter);
  virtual void paintVertex(typename T::Vertex_handle vh);

  QPen constraints_pen;
  QPen fixed_vertices_pen;

private:
  bool visible_constraints;

};

template <typename T>
void 
TriangulationArrangementGraphicsItem<T>::drawAll(QPainter *painter)
{
  this->painterostream = PainterOstream<Geom_traits>(painter);
  
  for(typename T::Finite_edges_iterator eit = this->t->finite_edges_begin();
      eit != this->t->finite_edges_end();
      ++eit)
  {
     
    if (   !eit->first->vertex(this->t->ccw(eit->second))->is_corner()
        || !eit->first->vertex(this->t->cw (eit->second))->is_corner()
       )
    {
     if(this->visibleConstraints() && this->t->is_constrained(*eit))
     {
       painter->setPen(constraintsPen());
       this->painterostream << this->t->segment(*eit);
     } 
     else if( this->visibleEdges() )
     {
       painter->setPen(this->edgesPen());
       this->painterostream << this->t->segment(*eit);
     }
    }
  }
  
  this->paintVertices(painter);
}

template <typename T>
void 
TriangulationArrangementGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  for (int i=0; i<3; i++) {
    if ( fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i)) )  {
      if(this->visibleConstraints() && this->t->is_constrained(typename T::Edge(fh,i)) && (! this->t->is_infinite(fh->neighbor(i)))){
        this->m_painter->setPen(constraintsPen());
	this->painterostream << this->t->segment(fh,i);
      } else if( this->visibleEdges() ){
	this->m_painter->setPen(this->edgesPen());
	this->painterostream << this->t->segment(fh,i);
      }
    }
    if(this->visibleVertices()) {
      paintVertex(fh->vertex(i));
    }
  }
}


template <typename T>
void 
TriangulationArrangementGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  if ( ! vh->is_corner() ) {
    Converter<Geom_traits> convert;
    
    if ( vh->is_fixed() ) {     
      this->m_painter->setPen(this->fixedVerticesPen());
    } else {
      this->m_painter->setPen(this->verticesPen());
    }
    QMatrix matrix = this->m_painter->matrix();
    this->m_painter->resetMatrix();
    this->m_painter->drawPoint(matrix.map(convert(vh->point())));
    this->m_painter->setMatrix(matrix);
  }
}


template <typename T>
void 
TriangulationArrangementGraphicsItem<T>::paintVertices(QPainter *painter)
{
  if(this->visibleVertices()) 
  {
    Converter<Geom_traits> convert;

    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename T::Finite_vertices_iterator it = this->t->finite_vertices_begin();
        it != this->t->finite_vertices_end();
        it++)
    {
      if ( ! it->is_corner() )
      {
        QPointF point = matrix.map(convert(it->point()));

        if ( it->is_fixed() )       
             painter->setPen(this->fixedVerticesPen());
        else painter->setPen(this->verticesPen());
        
        painter->drawPoint(point);
      }
    }
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_H
