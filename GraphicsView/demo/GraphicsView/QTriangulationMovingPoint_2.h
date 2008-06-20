
#ifndef CGAL_Q_TRIANGULATION_MOVING_POINT_2
#define CGAL_Q_TRIANGULATION_MOVING_POINT_2

#include <CGAL/IO/QtInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {

template <typename DT>
class QTriangulationMovingPoint_2 : public QtInput
{
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;

  QTriangulationMovingPoint_2(DT  * dt_);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  DT * dt;
  Vertex_handle vh;
  Face_handle fh;
  bool movePointToInsert;
};


template <typename T>
QTriangulationMovingPoint_2<T>::QTriangulationMovingPoint_2(T * dt_,
                                                            QObject* parent)
  :  QtInput(parent), dt(dt_), movePointToInsert(false), fh()
{}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::localize_and_insert_point(QPointF qt_point)
{
  Point p(qt_point.x(), qt_point.y());
  typename T::Locate_type lt;
  int li;
  fh = dt->locate(p, lt, li, fh); // fh serves as a hint

  vh = dt->insert(p, lt, fh, li);
  fh = vh->face(); // update the hint fh after the insertion
  emit(generate(CGAL::Object()));
}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->number_of_vertices() == 0){
    return;
  }
  movePointToInsert = true;
  localize_and_insert_point(event->scenePos());
}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{

  if(! movePointToInsert) return;

  // fh will be destroyed by the removal of vh.
  // Let us take a neighbor that is not in the star of vh.
  fh = fh->neighbor(fh->index(vh));
  dt->remove(vh);
  localize_and_insert_point(event->scenePos());
}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  if(! movePointToInsert) return;

  dt->remove(vh);
  fh = Face_handle();
  
  emit(generate(CGAL::Object()));
 
  movePointToInsert = false;
}



template <typename T>
bool 
QTriangulationMovingPoint_2<T>::eventFilter(QObject *obj, QEvent *event)
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



} // namespace CGAL

#endif // CGAL_Q_TRIANGULATION_MOVING_POINT_2
