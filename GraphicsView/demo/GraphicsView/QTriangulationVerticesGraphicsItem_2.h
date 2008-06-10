
#ifndef CGAL_TRIANGULATION_VERTICES_GRAPHICS_ITEM_2_H
#define CGAL_TRIANGULATION_VERTICES_GRAPHICS_ITEM_2_H

#include <CGAL/Bbox_2.h>
#include "QPainterOstream.h"
#include "QGraphicsItem_2.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>




class QGraphicsSceneMouseEvent;

namespace CGAL {


template <typename T>
class QTriangulationVerticesGraphicsItem_2 : public QGraphicsItem_2
{
public:
  QTriangulationVerticesGraphicsItem_2(T  * t_);

  void 
  operator()(typename T::Face_handle fh);
  
  QRectF 
  boundingRect() const;
  
  void 
  paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void 
  modelChanged();



private:
  void 
  updateBoundingBox();

  T * t;
  typename T::Vertex_handle vh;
  typename T::Point p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QPainter* m_painter;
};



template <typename T>
QTriangulationVerticesGraphicsItem_2<T>::QTriangulationVerticesGraphicsItem_2(T * t_)
  :  t(t_), bb(0,0,0,0), bb_initialized(false)
{
  if(t->number_of_vertices() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(4);

}

template <typename T>
QRectF 
QTriangulationVerticesGraphicsItem_2<T>::boundingRect() const
{
  return QRectF(bb.xmin()-1, bb.ymin()-1,
		(bb.xmax()-bb.xmin())+2, (bb.ymax()-bb.ymin())+2); // width height
}




template <typename T>
void 
QTriangulationVerticesGraphicsItem_2<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *)
{
  painter->setPen(pen());
  for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
      it != t->finite_vertices_end();
      it++){
    (*painter) << it->point();
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
QTriangulationVerticesGraphicsItem_2<T>::updateBoundingBox()
{
  if(t->number_of_vertices() == 0){
    bb = Bbox_2(0,0,0,0);
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
void 
QTriangulationVerticesGraphicsItem_2<T>::modelChanged()
{
  if((t->number_of_vertices() == 0) ){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    updateBoundingBox();
    this->show();
    update();
  }
}


} // namespace CGAL

#endif // CGAL_Q_TRIANGULATION_VERTICES_GRAPHICS_ITEM_2_H
