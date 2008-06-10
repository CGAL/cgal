
#ifndef CGAL_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H

#include "TriangulationGraphicsItem_2.h"



class QGraphicsSceneMouseEvent;

namespace CGAL {


template <typename T>
class ConstrainedTriangulationGraphicsItem_2 : public TriangulationGraphicsItem_2<T>
{
public:
  ConstrainedTriangulationGraphicsItem_2(T  * t_)
    : TriangulationGraphicsItem_2<T>(t_)
  {}

  
  virtual void operator()(typename T::Face_handle fh);
  
};



template <typename T>
 void 
ConstrainedTriangulationGraphicsItem_2<T>::operator()(typename T::Face_handle fh)
{
  QPen blackpen(Qt::black, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  QPen redpen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

  for (int i=0; i<3; i++)
    if (fh < fh->neighbor(i) || t->is_infinite(fh->neighbor(i))){
      if(t->is_constrained(typename T::Edge(fh,i))){
	m_painter->setPen(redpen);
      } else {
	m_painter->setPen(blackpen);
      }
      (*m_painter) << t->segment(fh,i);      
    }
}


} // namespace CGAL

#endif // CGAL_CONSTRAINED_TRIANGULATION_GRAPHICS_ITEM_2_H
