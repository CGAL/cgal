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
  typename DT::Face_handle fh;
  QGraphicsScene *scene_;
  QGraphicsEllipseItem* circle;
};


template <typename DT>
TriangulationCircumcircle<DT>::TriangulationCircumcircle(QGraphicsScene* s,
                                                         DT * dt_,
                                                         QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{
  circle = new QGraphicsEllipseItem();
  circle->hide();
  scene_->addItem(circle);
}


template <typename DT>
TriangulationCircumcircle<DT>::~TriangulationCircumcircle()
{
}


template <typename DT>
void
TriangulationCircumcircle<DT>::setPen(const QPen& pen)
{
  circle->setPen(pen);
}


template <typename DT>
void
TriangulationCircumcircle<DT>::show()
{
  circle->show();
}


template <typename DT>
void
TriangulationCircumcircle<DT>::hide()
{
  circle->hide();
}


template <typename DT>
void
TriangulationCircumcircle<DT>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
    if(dt->dimension() != 2){
        circle->hide();
        return;
    }
  typename DT::Point p = typename DT::Point(event->scenePos().x(), event->scenePos().y());

  double dx = dt->domain().xmax() - dt->domain().xmin();
  double dy = dt->domain().ymax() - dt->domain().ymin();
  p = typename DT::Point(p.x()- std::floor(p.x()/dx), p.y()- std::floor(p.y()/dy));

  fh = dt->locate(p);

  typename DT::Triangle triangle = dt->triangle(dt->periodic_triangle(fh));
  typename DT::Geom_traits::Circle_2 c(triangle[0],
                                       triangle[1],
                                       triangle[2]);
  CGAL::Bbox_2 bb = c.bbox();
  circle->setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
  circle->show();
}


template <typename DT>
bool
TriangulationCircumcircle<DT>::eventFilter(QObject *obj, QEvent *event)
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
