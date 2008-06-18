
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
  QTriangulationMovingPoint_2(DT  * dt_);

  void operator()(typename DT::Face_handle fh);
 
protected:

  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  DT * dt;
  typename DT::Vertex_handle vh;
  typename DT::Point p;
  bool movePointToInsert;
  CGAL::Bbox_2 bb;  
};




template <typename T>
QTriangulationMovingPoint_2<T>::QTriangulationMovingPoint_2(T * dt_)
  :  dt(dt_), movePointToInsert(false)
{}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->number_of_vertices() == 0){
    return;
  }
  p = typename T::Point(event->scenePos().x(), event->scenePos().y());
  movePointToInsert = true;
  typename T::Locate_type lt;
  int li;
  typename T::Face_handle fh = dt->locate(p, lt, li);

  vh = dt->insert(p, lt, fh, li);
  emit(generate(CGAL::Object()));
}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{

  if(! movePointToInsert) return;

  dt->remove(vh);
  p = typename T::Point(event->scenePos().x(), event->scenePos().y());
  typename T::Locate_type lt;
  int li;
  typename T::Face_handle fh = dt->locate(p, lt, li);
 
  vh = dt->insert(p, lt, fh, li);
  emit(generate(CGAL::Object()));

}


template <typename T>
void 
QTriangulationMovingPoint_2<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  if(! movePointToInsert) return;
  dt->remove(vh);
  
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



} // namespace CGAL

#endif // CGAL_Q_TRIANGULATION_MOVING_POINT_2
