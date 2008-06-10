#ifndef CGAL_TRIANGULATION_CIRCUMCENTER_2
#define CGAL_TRIANGULATION_CIRCUMCENTER_2

#include "Input.h"
#include <QGraphicsSceneMouseEvent> 
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QEvent>
#include <QPen>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


namespace CGAL {

template <typename DT>
class TriangulationCircumcenter_2 : public Input
{


public:
  TriangulationCircumcenter_2(QGraphicsScene* s, DT  * dt_);
 
  void setPen(const QPen& pen);

  void show();
  void hide();
  
protected:

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  DT * dt;
  typename DT::Face_handle fh;
  QGraphicsScene *scene_;
  QGraphicsEllipseItem circle;
};


template <typename T>
TriangulationCircumcenter_2<T>::TriangulationCircumcenter_2(QGraphicsScene* s, T * dt_)
  :  dt(dt_), scene_(s)
{
  circle.hide();
  scene_->addItem(&circle);
}


template <typename T>
void
TriangulationCircumcenter_2<T>::setPen(const QPen& pen)
{
  circle.setPen(pen);
}

template <typename T>
void
TriangulationCircumcenter_2<T>::show()
{
  circle.show();
}

template <typename T>
void
TriangulationCircumcenter_2<T>::hide()
{
  circle.hide();
}


template <typename T>
void TriangulationCircumcenter_2<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->dimension() != 2){
    circle.hide();
    return;
  }
  typename T::Point p = typename T::Point(event->scenePos().x(), event->scenePos().y());
  typename T::Face_handle fh = dt->locate(p);
  if(!dt->is_infinite(fh)){
    CGAL::Exact_predicates_inexact_constructions_kernel::Circle_2 c(fh->vertex(0)->point(), 
								    fh->vertex(1)->point(), 
								    fh->vertex(2)->point());
    CGAL::Bbox_2 bb = c.bbox();
    circle.setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
    circle.show();
  } else {
    circle.hide();
  }
}


template <typename T>
bool TriangulationCircumcenter_2<T>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // don't consume the event
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
} 


} // namespace CGAL

#endif // CGAL_TRIANGULATION_CIRCUMCENTER_2
