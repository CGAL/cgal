
#ifndef CGAL_QT_TRIANGULATION_MOVING_POINT
#define CGAL_QT_TRIANGULATION_MOVING_POINT

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationMovingPoint : public GraphicsViewInput
{
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;

  TriangulationMovingPoint(DT  * dt_, QObject* parent);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  DT * dt;
  Vertex_handle vh;
  bool movePointToInsert;
  bool insertedPoint;
};


template <typename T>
TriangulationMovingPoint<T>::TriangulationMovingPoint(T * dt_,
							  QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), vh(), movePointToInsert(false), insertedPoint(false)
{}


template <typename T>
void 
TriangulationMovingPoint<T>::localize_and_insert_point(QPointF qt_point)
{
  Point p(qt_point.x(), qt_point.y());
  typename T::Locate_type lt;
  int li;
  Face_handle fh = (vh == Vertex_handle()) ? Face_handle() : vh->face();
  fh = dt->locate(p, lt, li, fh); // fh serves as a hint
  if(lt != T::VERTEX){
    vh = dt->insert(p, lt, fh, li);
    insertedPoint = true;
    emit(modelChanged());
  } else {
    vh = fh->vertex(0);
    insertedPoint = false;
  }
}


template <typename T>
void 
TriangulationMovingPoint<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->number_of_vertices() == 0 ||
     event->modifiers() != 0 ||
     event->button() != ::Qt::LeftButton) {
    return;
  }
  movePointToInsert = true;
  localize_and_insert_point(event->scenePos());
}


template <typename T>
void 
TriangulationMovingPoint<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{

  if(! movePointToInsert) return;

  // fh will be destroyed by the removal of vh.
  // Let us take a neighbor that is not in the star of vh.
  const Face_handle fh = vh->face();
  Vertex_handle next_hint = fh->vertex((fh->index(vh)+1)&3);
  if(insertedPoint){
    dt->remove(vh);
  }
  vh = next_hint;
  localize_and_insert_point(event->scenePos());
}


template <typename T>
void 
TriangulationMovingPoint<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  if(! movePointToInsert ||
     event->button() != ::Qt::LeftButton) {
    return;
  }

  if(insertedPoint){
    dt->remove(vh);
  }
  vh = Vertex_handle();
  
  emit(modelChanged());
 
  movePointToInsert = false;
}



template <typename T>
bool 
TriangulationMovingPoint<T>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_QT_TRIANGULATION_MOVING_POINT
