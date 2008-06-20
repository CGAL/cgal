#ifndef CGAL_Q_TRIANGULATION_CIRCUMCENTER_2
#define CGAL_Q_TRIANGULATION_CIRCUMCENTER_2

#include <CGAL/IO/QtInput.h>
#include <QGraphicsSceneMouseEvent> 
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QEvent>
#include <QPen>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


namespace CGAL {

template <typename DT>
class QTriangulationCircumcenter_2 : public QtInput
{
public:
  QTriangulationCircumcenter_2(QGraphicsScene* s, DT  * dt_, QObject* parent);
  ~QTriangulationCircumcenter_2();
 
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
  QGraphicsEllipseItem* circle;
};


template <typename T>
QTriangulationCircumcenter_2<T>::QTriangulationCircumcenter_2(QGraphicsScene* s,
                                                              T * dt_
                                                              QObject* parent)
  :  QtInput(parent), dt(dt_), scene_(s)
{
  circle = new QGraphicsEllipseItem();
  circle->hide();
  scene_->addItem(circle);
}

template <typename T>
QTriangulationCircumcenter_2<T>::~QTriangulationCircumcenter_2()
{
  delete circle;
}

template <typename T>
void
QTriangulationCircumcenter_2<T>::setPen(const QPen& pen)
{
  circle->setPen(pen);
}

template <typename T>
void
QTriangulationCircumcenter_2<T>::show()
{
  circle->show();
}

template <typename T>
void
QTriangulationCircumcenter_2<T>::hide()
{
  circle->hide();
}


template <typename T>
void 
QTriangulationCircumcenter_2<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->dimension() != 2){
    circle->hide();
    return;
  }
  typename T::Point p = typename T::Point(event->scenePos().x(), event->scenePos().y());
  fh = dt->locate(p, fh); // fh is used as a hint
  if(!dt->is_infinite(fh)){
    typename T::Geom_traits::Circle_2 c(fh->vertex(0)->point(), 
                                        fh->vertex(1)->point(), 
                                        fh->vertex(2)->point());
    CGAL::Bbox_2 bb = c.bbox();
    circle->setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
    circle->show();
  } else {
    circle->hide();
  }
}


template <typename T>
bool 
QTriangulationCircumcenter_2<T>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_Q_TRIANGULATION_CIRCUMCENTER_2
