
#ifndef CGAL_QT_PERIODIC_TRIANGULATION_CONFLICT_ZONE
#define CGAL_QT_PERIODIC_TRIANGULATION_CONFLICT_ZONE

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <queue>



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

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  Periodic_triangulation *m_tr;
  Face_handle m_containing_face;
  Converter<K> m_convert;
  QGraphicsScene *m_scene;
  QGraphicsPolygonItem *m_triangle;
};


template <typename T>
TriangulationConflictZone<T>::TriangulationConflictZone(QGraphicsScene* s,
                                                            T * tr_,
                                                            QObject* parent)
  :  GraphicsViewInput(parent), m_tr(tr_), m_containing_face(Face_handle()), m_scene(s), m_triangle(NULL)
{
}


template <typename T>
void 
TriangulationConflictZone<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
}

template <typename T>
void 
TriangulationConflictZone<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if (m_triangle) {
    m_scene->removeItem(m_triangle);
    m_triangle = NULL;
  }

  /// Do a point location
  Point p = m_convert(event->scenePos());
  const Iso_rectangle & domain = m_tr->domain();
  const Point min_coord = domain.min();
  const Vector domain_size = domain.max() - min_coord;
  p = Point(p.x() - std::floor((p.x()-min_coord.x())/domain_size.x()),
            p.y() - std::floor((p.y()-min_coord.y())/domain_size.y()));

  m_containing_face = m_tr->locate(p);

  CGAL_assertion(m_containing_face != typename Periodic_triangulation::Face_handle());
  if (m_containing_face != typename Periodic_triangulation::Face_handle()) {
    // Display the triangle
    m_triangle = new QGraphicsPolygonItem(m_convert(m_tr->triangle(m_containing_face)));
    m_triangle->setBrush(QColor(::Qt::green));
    
    m_scene->addItem(m_triangle);
  }
}


template <typename T>
void 
TriangulationConflictZone<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
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
