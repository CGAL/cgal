
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
public:
  ConstrainedTriangulationGraphicsItem(T  * t_)
    : TriangulationGraphicsItem<T>(t_)
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

protected:
  void drawAll(QPainter *painter);

  QPen constraints_pen;
};

template <typename T>
void 
ConstrainedTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  if(this->drawEdges()) {
    for(typename T::Finite_edges_iterator eit = this->t->finite_edges_begin();
        eit != this->t->finite_edges_end();
        ++eit){
      if(this->t->is_constrained(*eit)){
        painter->setPen(constraintsPen());
      } else {
        painter->setPen(this->edgesPen());
      }
      (*painter) << this->t->segment(*eit);
    }
  }
  this->paintVertices(painter);
}

template <typename T>
void 
ConstrainedTriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  for (int i=0; i<3; i++) {
    if (this->drawEdges() &&
        ( fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i)) ) ) {
      if(this->t->is_constrained(typename T::Edge(fh,i))){
        this->m_painter->setPen(constraintsPen());
      } else {
	this->m_painter->setPen(this->edgesPen());
      }
      (*this->m_painter) << this->t->segment(fh,i);      
    }
    if(this->drawVertices()) {
      paintOneVertex(fh->vertex(i)->point());
    }
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_H
