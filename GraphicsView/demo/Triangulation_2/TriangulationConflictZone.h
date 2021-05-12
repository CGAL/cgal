
#ifndef CGAL_QT_TRIANGULATION_CONFLICT_ZONE
#define CGAL_QT_TRIANGULATION_CONFLICT_ZONE

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationConflictZone : public GraphicsViewInput
{
public:
  typedef typename DT::Geom_traits K;
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Point Point;

  TriangulationConflictZone(QGraphicsScene* s, DT  * dt_, QObject* parent);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  std::list<Face_handle> faces;
  std::list<QGraphicsPolygonItem*> qfaces;
  DT * dt;
  Converter<K> convert;
  QGraphicsScene *scene_;
  bool animate;
  Face_handle hint;
};


template <typename T>
TriangulationConflictZone<T>::TriangulationConflictZone(QGraphicsScene* s,
                                                        T * dt_,
                                                        QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s), animate(false)
{}


template <typename T>
void
TriangulationConflictZone<T>::localize_and_insert_point(QPointF qt_point)
{
  Point p(convert(qt_point));

  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    delete *it;
  }
  qfaces.clear();
  hint = dt->locate(p, hint);
  dt->find_conflicts(p, faces, hint);
  for(typename std::list<Face_handle>::iterator it = faces.begin();
      it != faces.end();
      ++it){
    if(! dt->is_infinite(*it)){
      QGraphicsPolygonItem *item = new QGraphicsPolygonItem(convert(dt->triangle(*it)));
      QColor color(::Qt::blue);
      color.setAlpha(150);
      item->setBrush(color);
      scene_->addItem(item);
      qfaces.push_back(item);
    }
  }
}



template <typename T>
void
TriangulationConflictZone<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->number_of_vertices() == 0 ||
     event->modifiers() != 0 ||
     event->button() != ::Qt::LeftButton) {
    return;
  }
  hint = dt->locate(convert(event->scenePos()));
  localize_and_insert_point(event->scenePos());
  animate = true;
}


template <typename T>
void
TriangulationConflictZone<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(animate){
    localize_and_insert_point(event->scenePos());
  }
}


template <typename T>
void
TriangulationConflictZone<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent * /*event*/)
{
  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    delete *it;
  }
  qfaces.clear();
  animate = false;
}



template <typename T>
bool
TriangulationConflictZone<T>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // do not eat move event!
  } else if (event->type() == QEvent::GraphicsSceneMouseRelease) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseReleaseEvent(mouseEvent);
    return true;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_CONFLICT_ZONE
