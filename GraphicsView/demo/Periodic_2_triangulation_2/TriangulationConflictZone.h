#ifndef CGAL_QT_PERIODIC_TRIANGULATION_CONFLICT_ZONE
#define CGAL_QT_PERIODIC_TRIANGULATION_CONFLICT_ZONE

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsPolygonItem>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {
namespace Qt {

template <typename PT>
class TriangulationConflictZone : public GraphicsViewInput
{
public:
  typedef PT                                       Periodic_triangulation;
  typedef typename PT::Geom_traits                 K;
  typedef typename PT::Face_handle                 Face_handle;
  typedef typename PT::Periodic_triangle_iterator  Periodic_triangle_iterator;
  typedef typename PT::Point                       Point;
  typedef typename PT::Iso_rectangle               Iso_rectangle;
  typedef typename PT::Geom_traits::Vector_2       Vector;
  typedef typename PT::Triangle                    Triangle;

  TriangulationConflictZone(QGraphicsScene* s, PT  * tr_, QObject* parent);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  Face_handle m_hint;
  std::list<Face_handle> faces;
  std::list<QGraphicsPolygonItem*> qfaces;
  Periodic_triangulation *m_tr;
  Face_handle m_containing_face;
  Converter<K> m_convert;
  QGraphicsScene *m_scene;
  QGraphicsPolygonItem *m_triangle;
  bool m_animate;
};


template <typename T>
TriangulationConflictZone<T>::TriangulationConflictZone(QGraphicsScene* s,
                                                            T * tr_,
                                                            QObject* parent)
    :  GraphicsViewInput(parent), m_hint(NULL), m_tr(tr_), m_containing_face(Face_handle()), m_scene(s), m_triangle(NULL), m_animate(false)
{
}

template <typename T>
void
TriangulationConflictZone<T>::localize_and_insert_point(QPointF qt_point)
{
  Point p(m_convert(qt_point));
  double dx = m_tr->domain().xmax() - m_tr->domain().xmin();
  double dy = m_tr->domain().ymax() - m_tr->domain().ymin();
  p = Point(p.x()- std::floor(p.x()/dx), p.y()- std::floor(p.y()/dy));

  if (m_hint == NULL) {
      m_hint = m_tr->faces_begin();
  }

  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    delete *it;
  }
  qfaces.clear();
  m_tr->get_conflicts(p, std::back_inserter(faces), m_hint);
  for(typename std::list<Face_handle>::iterator it = faces.begin();
      it != faces.end();
      ++it){
    if(! m_tr->is_infinite(*it)){
      QGraphicsPolygonItem *item = new QGraphicsPolygonItem(m_convert(m_tr->triangle(*it)));
      QColor color = ::Qt::blue;
      color.setAlpha(150);
      item->setBrush(QBrush(color));
      item->setPen(QPen(::Qt::black, .01));
      m_scene->addItem(item);
      qfaces.push_back(item);
    }
  }
}



template <typename T>
void
TriangulationConflictZone<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(m_tr->number_of_vertices() == 0 ||
     event->modifiers() != 0 ||
     event->button() != ::Qt::LeftButton) {
    return;
  }
  localize_and_insert_point(event->scenePos());
  m_animate = true;
}


template <typename T>
void
TriangulationConflictZone<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(m_animate){
    localize_and_insert_point(event->scenePos());
  }
}


template <typename T>
void
TriangulationConflictZone<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *)
{
  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    delete *it;
  }
  qfaces.clear();
  m_animate = false;
  m_hint = NULL;
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

#endif // CGAL_QT_PERIODIC_TRIANGULATION_CONFLICT_ZONE
