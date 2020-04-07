#ifndef CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H
#define CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H

#include <QGraphicsSceneMouseEvent>
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QEvent>
#include <QPen>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Qt/GraphicsViewInput.h>

namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationCircumcircle : public GraphicsViewInput
{
public:
  TriangulationCircumcircle(QGraphicsScene* s, DT  * dt_, QObject* parent);
  ~TriangulationCircumcircle();

  void setPen(const QPen& pen);

  void show();
  void hide();

protected:

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  DT * dt;
  typedef typename DT::Vertex_handle Vertex_handle;
  typename DT::Vertex_handle hint;
  typename DT::Face_handle fh;
  QGraphicsScene *scene_;
  QGraphicsEllipseItem* circle;
};


template <typename T>
TriangulationCircumcircle<T>::TriangulationCircumcircle(QGraphicsScene* s,
                                                              T * dt_,
                                                              QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{
  hint = dt->infinite_vertex();
  circle = new QGraphicsEllipseItem();
  circle->hide();
  scene_->addItem(circle);
}


template <typename T>
TriangulationCircumcircle<T>::~TriangulationCircumcircle()
{
}


template <typename T>
void
TriangulationCircumcircle<T>::setPen(const QPen& pen)
{
  circle->setPen(pen);
}


template <typename T>
void
TriangulationCircumcircle<T>::show()
{
  circle->show();
}


template <typename T>
void
TriangulationCircumcircle<T>::hide()
{
  circle->hide();
}


template <typename T>
void
TriangulationCircumcircle<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->dimension() != 2){
    circle->hide();
    hint = Vertex_handle();
    return;
  }
  if (hint == Vertex_handle()){
    hint = dt->infinite_vertex();
  }
  typename T::Point p = typename T::Point(event->scenePos().x(), event->scenePos().y());
  fh = dt->locate(p, hint->face());
  hint = fh->vertex(0);
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
TriangulationCircumcircle<T>::eventFilter(QObject *obj, QEvent *event)
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

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H
