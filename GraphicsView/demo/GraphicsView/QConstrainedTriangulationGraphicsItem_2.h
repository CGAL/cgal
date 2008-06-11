
#ifndef CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
#define CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H

#include "QTriangulationGraphicsItem_2.h"



class QGraphicsSceneMouseEvent;

namespace CGAL {


template <typename T>
class QConstrainedTriangulationGraphicsItem_2 : public QTriangulationGraphicsItem_2<T>
{
public:
  QConstrainedTriangulationGraphicsItem_2(T  * t_)
    : QTriangulationGraphicsItem_2<T>(t_)
  {}

  
  virtual void operator()(typename T::Face_handle fh);
  
};



template <typename T>
void 
QConstrainedTriangulationGraphicsItem_2<T>::operator()(typename T::Face_handle fh)
{
  QPen blackpen(Qt::black, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  QPen redpen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

  for (int i=0; i<3; i++)
    if (fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i))){
      if(this->t->is_constrained(typename T::Edge(fh,i))){
	this->m_painter->setPen(redpen);
      } else {
	this->m_painter->setPen(blackpen);
      }
      (*this->m_painter) << this->t->segment(fh,i);      
    }
}


} // namespace CGAL

#endif // CGAL_Q_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
