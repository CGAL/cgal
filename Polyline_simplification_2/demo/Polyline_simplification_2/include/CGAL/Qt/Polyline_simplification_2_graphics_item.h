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

#ifndef CGAL_QT_POLYLINE_SIMPLIFICATION_2_GRAPHICS_ITEM_H
#define CGAL_QT_POLYLINE_SIMPLIFICATION_2_GRAPHICS_ITEM_H

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <QPen>

namespace CGAL {

namespace Qt {

template <typename PCT>
class PolylineSimplificationGraphicsItem : public TriangulationGraphicsItem<PCT >
{
  typedef typename PCT::Geom_traits Geom_traits;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Constraint_iterator Constraint_iterator;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  //typedef typename PCT::Points_in_constraint_iterator Points_in_constraint_iterator;

  Vertices_in_constraint_iterator
  decrement(Vertices_in_constraint_iterator it)
  {
    do{ 
      --it;
    } while(it->removed);
    return it;
  }
  
  Vertices_in_constraint_iterator
  increment(Vertices_in_constraint_iterator it)
  {
    do{ 
      ++it;
    } while(it->removed);
    return it;
  }
  

  
public:

  PolylineSimplificationGraphicsItem(PCT  * pct)
    : TriangulationGraphicsItem<PCT>(pct), visible_constraints(true)
  {
    constraints_pen = this->edgesPen();
    constraints_pen.setColor(::Qt::red);
    
    unremovable_vertices_pen = this->verticesPen();
    unremovable_vertices_pen.setColor(::Qt::blue);
  }
  
  void operator()(typename PCT::Face_handle fh);

  const QPen& constraintsPen() const
  {
    return constraints_pen;
  }

  void setConstraintsPen(const QPen& pen)
  {
    constraints_pen = pen;
  }

  const QPen& unremovableVerticesPen() const
  {
    return unremovable_vertices_pen;
  }

  void setUnremovableVerticesPen(const QPen& pen)
  {
    unremovable_vertices_pen = pen;
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
  virtual void paintVertex(typename PCT::Vertex_handle vh);

  QPen constraints_pen;
  QPen unremovable_vertices_pen;

private:
  bool visible_constraints;

};

template <typename PCT>
void 
PolylineSimplificationGraphicsItem<PCT>::drawAll(QPainter *painter)
{
  this->painterostream = PainterOstream<Geom_traits>(painter);
  
  for(typename PCT::Finite_edges_iterator eit = this->t->finite_edges_begin();
      eit != this->t->finite_edges_end();
      ++eit)
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
  
  this->paintVertices(painter);
}

template <typename PCT>
void 
PolylineSimplificationGraphicsItem<PCT>::operator()( typename PCT::Face_handle fh )
{
  for (int i=0; i<3; i++) {
    if ( fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i)) )  {
//      if(this->visibleConstraints() && this->t->is_constrained(typename PCT::Edge(fh,i)) && (! this->t->is_infinite(fh->neighbor(i)))){
      if(this->visibleConstraints() && this->t->is_constrained(typename PCT::Edge(fh,i)) ){
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


template <typename PCT>
void 
PolylineSimplificationGraphicsItem<PCT>::paintVertex( typename PCT::Vertex_handle vh )
{
  Converter<Geom_traits> convert;
  
  //  if ( vh->is_unremovable() || vh->is_shared() ) {     
  //    this->m_painter->setPen(this->unremovableVerticesPen());
  // } else {
    this->m_painter->setPen(this->verticesPen());
    //  }
  QMatrix matrix = this->m_painter->matrix();
  this->m_painter->resetMatrix();
  this->m_painter->drawPoint(matrix.map(convert(vh->point())));
  this->m_painter->setMatrix(matrix);
}


template <typename PCT>
void 
PolylineSimplificationGraphicsItem<PCT>::paintVertices(QPainter *painter)
{
  if(this->visibleVertices()) 
  {
    Converter<Geom_traits> convert;

    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(Constraint_iterator cit = this->t->constraints_begin();
        cit != this->t->constraints_end();
        ++cit){
      for(Vertices_in_constraint_iterator it = this->t->vertices_in_constraint_begin(*cit);
        it != this->t->vertices_in_constraint_end(*cit);
        it++){
        QPointF point = matrix.map(convert((*it)->point()));  
        if ( (*it)->is_removable() )       
          painter->setPen(this->verticesPen());
        else 
          painter->setPen(this->unremovableVerticesPen());

        painter->drawPoint(point);
      }
    }
  }
}
  
} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POLYLINE_SIMPLIFICATION_2_GRAPHICS_ITEM_H
