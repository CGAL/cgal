
#ifndef CGAL_QT_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <QPen>

class QGraphicsSceneMouseEvent;

namespace CGAL {

namespace Qt {

template <typename T>
class ConstrainedTriangulationGraphicsItem : public TriangulationGraphicsItem<T>
{
  typedef typename T::Geom_traits Geom_traits;
public:
  ConstrainedTriangulationGraphicsItem(T  * t_)
    : TriangulationGraphicsItem<T>(t_), visible_constraints(true)
  {
    constraints_pen = this->edgesPen();
    constraints_pen.setColor(::Qt::red);
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

  bool visibleConstraints() const
  {
    return visible_constraints;
  }

  void setvisibleConstraints(const bool b)
  {
    visible_constraints = b;
    this->update();
  }

protected:
  void drawAll(QPainter *painter);

  QPen constraints_pen;

private:
  bool visible_constraints;

};

template <typename T>
void 
ConstrainedTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  this->painterostream = PainterOstream<Geom_traits>(painter);
  for(typename T::Finite_edges_iterator eit = this->t->finite_edges_begin();
      eit != this->t->finite_edges_end();
      ++eit){
    if(this->t->is_constrained(*eit)){
      painter->setPen(constraintsPen());
    } else {
      painter->setPen(this->edgesPen());
    }
    if(this->visibleEdges() || this->visibleConstraints()){
      this->painterostream << this->t->segment(*eit);
    }
  }
  
  this->paintVertices(painter);
}

template <typename T>
void 
ConstrainedTriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  for (int i=0; i<3; i++) {
    if ( fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i)) )  {
      if(this->t->is_constrained(typename T::Edge(fh,i))){
        this->m_painter->setPen(constraintsPen());
      } else {
	this->m_painter->setPen(this->edgesPen());
      }
      if(this->visibleEdges() || this->visibleConstraints()){
	this->painterostream << this->t->segment(fh,i);
      }   
    }
    if(this->visibleVertices()) {
      paintOneVertex(fh->vertex(i)->point());
    }
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_H
