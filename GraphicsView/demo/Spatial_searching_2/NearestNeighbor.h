#ifndef CGAL_QT_NEAREST_NEIGHBOR_H
#define CGAL_QT_NEAREST_NEIGHBOR_H

#include <vector>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QEvent>
#include <QPen>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Qt/GraphicsViewInput.h>

#include <CGAL/Qt/PointsGraphicsItem.h>

namespace CGAL {
namespace Qt {

template <typename T>
class NearestNeighbor : public GraphicsViewInput
{
  typedef typename T::Point_d Point_2;
public:
  NearestNeighbor(QGraphicsScene* s, typename T::Tree  * t_, QObject* parent, int N);
  ~NearestNeighbor();

  void setPen(const QPen& pen);

  void show();
  void hide();

  void setN(int n)
  {
    N = n;
  }

protected:

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  int N;
  typename T::Tree * t;
  QGraphicsScene *scene_;
  QGraphicsEllipseItem* circle;
  std::vector<Point_2> points;
  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * pgi;
};


template <typename T>
NearestNeighbor<T>::NearestNeighbor(QGraphicsScene* s,
                                    typename T::Tree * t,
                                    QObject* parent,
                                    int N)
  :  GraphicsViewInput(parent), N(N), t(t), scene_(s)
{
  pgi = new PointsGraphicsItem<std::vector<Point_2> >(&points);
  pgi->setVerticesPen(QPen(::Qt::red, 3, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
  pgi->setZValue(10);
  scene_->addItem(pgi);

  circle = new QGraphicsEllipseItem();
  circle->hide();
  scene_->addItem(circle);
}


template <typename T>
NearestNeighbor<T>::~NearestNeighbor()
{
}


template <typename T>
void
NearestNeighbor<T>::setPen(const QPen& pen)
{
  circle->setPen(pen);
}


template <typename T>
void
NearestNeighbor<T>::show()
{
  circle->show();
}


template <typename T>
void
NearestNeighbor<T>::hide()
{
  circle->hide();
}


template <typename T>
void
NearestNeighbor<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(t->empty()){
    points.clear();
    circle->hide();
    return;
  }
  Point_2 query(event->scenePos().x(), event->scenePos().y());
  T search(*t, query, N);

  points.clear();
  double sr = 0;
  for(typename T::iterator it = search.begin(); it != search.end(); ++it){
    if(sr < it->second){
      points.push_back(it->first);
      sr = it->second;
    }
  }
  pgi->modelChanged();
  typename CGAL::Kernel_traits<Point_2>::Kernel::Circle_2 c(query, sr);
  CGAL::Bbox_2 bb = c.bbox();
  circle->setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
  circle->show();
}


template <typename T>
bool
NearestNeighbor<T>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_QT_NEAREST_NEIGHBOR_H
