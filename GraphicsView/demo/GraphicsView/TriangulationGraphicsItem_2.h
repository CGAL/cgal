
#ifndef CGAL_TRIANGULATION_GRAPHICS_ITEM_2_H
#define CGAL_TRIANGULATION_GRAPHICS_ITEM_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/apply_to_range.h>
#include "QPainterOstream.h"
#include "GraphicsItem_2.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>




class QGraphicsSceneMouseEvent;

namespace CGAL {


template <typename T>
class TriangulationGraphicsItem_2 : public GraphicsItem_2
{
public:
  TriangulationGraphicsItem_2(T* t_);

  enum { Type = UserType + 4 };
  int type() const { return Type; }
  
  
  virtual void operator()(typename T::Face_handle fh);
  
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void vModelChanged();

protected:
  T * t;
  QPainter* m_painter;

private:
  void updateBoundingBox();

  typename T::Vertex_handle vh;
  typename T::Point p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
};



template <typename T>
TriangulationGraphicsItem_2<T>::TriangulationGraphicsItem_2(T * t_)
  :  t(t_), bb_initialized(false)
{
  updateBoundingBox();
  setZValue(3);

}

template <typename T>
 void 
TriangulationGraphicsItem_2<T>::operator()(typename T::Face_handle fh)
{
  m_painter->setPen(this->pen());
  for (int i=0; i<3; i++)
    if (fh < fh->neighbor(i) || t->is_infinite(fh->neighbor(i))){
      (*m_painter) << t->segment(fh,i);
    }
}

template <typename T>
QRectF TriangulationGraphicsItem_2<T>::boundingRect() const
{
  return QRectF(bb.xmin()-1, bb.ymin()-1,
		(bb.xmax()-bb.xmin())+2, (bb.ymax()-bb.ymin())+2); // width height
}




template <typename T>
void TriangulationGraphicsItem_2<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *)
{
  if (t->dimension()<2) {
    for(typename T::Finite_edges_iterator eit = t->finite_edges_begin();
	eit != t->finite_edges_end();
	++eit){
      (*painter) << t->segment(*eit);
    }
  } else {
    m_painter = painter;
    CGAL::apply_to_range(*t, 
			 typename T::Point(option->exposedRect.x(),
					    option->exposedRect.y()+option->exposedRect.height()), 
			 typename T::Point(option->exposedRect.x()+option->exposedRect.width(), 
					    option->exposedRect.y()), 
			 *this);
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void TriangulationGraphicsItem_2<T>::updateBoundingBox()
{
  if(t->number_of_vertices() == 0){
    bb = Bbox_2();
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
}


template <typename T>
void TriangulationGraphicsItem_2<T>::vModelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace CGAL

#endif // CGAL_TRIANGULATION_GRAPHICS_ITEM_2_H
