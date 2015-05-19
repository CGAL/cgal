
#ifndef CGAL_QT_TRIANGULATION_POINT_INPUT_AND_CONFLICT_ZONE
#define CGAL_QT_TRIANGULATION_POINT_INPUT_AND_CONFLICT_ZONE

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {
namespace Qt {

template <typename PT>
class TriangulationPointInputAndConflictZone : public GraphicsViewInput
{
public:
  typedef typename PT::Geom_traits K;
  typedef typename PT::Face_handle Face_handle;
  typedef typename PT::Point Point;

  TriangulationPointInputAndConflictZone(QGraphicsScene* s, PT  * dt_, QObject* parent);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  std::list<Face_handle> faces;
  std::list<QGraphicsPolygonItem*> qfaces;
  PT * dt;
  Converter<K> convert;
  QGraphicsScene *scene_;
  Point p;
};


template <typename T>
TriangulationPointInputAndConflictZone<T>::TriangulationPointInputAndConflictZone(QGraphicsScene* s,
							T * dt_,
							QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{}




template <typename T>
void 
TriangulationPointInputAndConflictZone<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  p = convert(event->scenePos());
  double dx = dt->domain().xmax() - dt->domain().xmin();
  double dy = dt->domain().ymax() - dt->domain().ymin();
  p = Point(p.x()- std::floor(p.x()/dx), p.y()- std::floor(p.y()/dy));

  if(dt->dimension() < 2 ||
     event->modifiers() != 0 ||
     event->button() != ::Qt::LeftButton) {
    return;
  }
  

  dt->get_conflicts(p, std::back_inserter(faces));
  for(typename std::list<Face_handle>::iterator it = faces.begin();
      it != faces.end();
      ++it) {
      QGraphicsPolygonItem *item = new QGraphicsPolygonItem(convert(dt->triangle(*it)));
      QColor color(::Qt::blue);
      color.setAlpha(150);
      item->setBrush(color);
      scene_->addItem(item);
      qfaces.push_back(item);
  }
}


template <typename T>
void 
TriangulationPointInputAndConflictZone<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *)
{
  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    scene_->removeItem(*it);
    delete *it;
  }
  qfaces.clear();
  Q_EMIT( generate(CGAL::make_object(p)));
}



template <typename T>
bool 
TriangulationPointInputAndConflictZone<T>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
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

#endif // CGAL_QT_TRIANGULATION_POINT_INPUT_AND_CONFLICT_ZONE
