
#ifndef CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
#define CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H

#include "QTriangulationGraphicsItem_2.h"
#include <QPen>

class QGraphicsSceneMouseEvent;

namespace CGAL {


template <typename T>
class QConstrainedTriangulationGraphicsItem_2 : public QTriangulationGraphicsItem_2<T>
{
public:
  QConstrainedTriangulationGraphicsItem_2(T  * t_)
    : QTriangulationGraphicsItem_2<T>(t_)
  {
    constraints_pen = this->pen();
    constraints_pen.setColor(Qt::red);
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

protected:
  void drawAll(QPainter *painter);

  QPen constraints_pen;
};

template <typename T>
void 
QConstrainedTriangulationGraphicsItem_2<T>::drawAll(QPainter *painter)
{
  for(typename T::Finite_edges_iterator eit = this->t->finite_edges_begin();
      eit != this->t->finite_edges_end();
      ++eit){
    if(this->t->is_constrained(*eit)){
      painter->setPen(constraintsPen());
    } else {
      painter->setPen(this->pen());
    }
    (*painter) << this->t->segment(*eit);
  }
  this->paintVertices(painter);
}

template <typename T>
void 
QConstrainedTriangulationGraphicsItem_2<T>::operator()(typename T::Face_handle fh)
{
  for (int i=0; i<3; i++) {
    if (fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i))){
      if(this->t->is_constrained(typename T::Edge(fh,i))){
	this->m_painter->setPen(constraintsPen());
      } else {
	this->m_painter->setPen(this->pen());
      }
      (*this->m_painter) << this->t->segment(fh,i);      
    }
    if(this->drawVertices()) {
      paintOneVertex(fh->vertex(i)->point());
    }
  }
}

} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
