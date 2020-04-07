
#ifndef CGAL_QT_REGULAR_TRIANGULATION_REMOVE_VERTEX_H
#define CGAL_QT_REGULAR_TRIANGULATION_REMOVE_VERTEX_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>
#include <CGAL/Qt/Converter.h>



namespace CGAL {
namespace Qt {

template <typename DT>
class RegularTriangulationRemoveVertex : public GraphicsViewInput
{
public:
  typedef typename DT::Vertex_handle Vertex_handle;

  RegularTriangulationRemoveVertex(DT  * dt_, QObject* parent);

protected:

  void mousePressEvent(QGraphicsSceneMouseEvent *event);

  bool eventFilter(QObject *obj, QEvent *event);

  DT * dt;
};


template <typename T>
RegularTriangulationRemoveVertex<T>::RegularTriangulationRemoveVertex(T * dt_,
                                                          QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_)
{}



template <typename T>
void
RegularTriangulationRemoveVertex<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  if((event->modifiers()  & ::Qt::ShiftModifier)
     && (! (event->modifiers() & ::Qt::ControlModifier))){
    if(dt->number_of_vertices() == 0){
      dt->clear();
    }else {
      typedef typename Kernel_traits<typename T::Bare_point>::Kernel K;
      Converter<K> convert;
      typename T::Vertex_handle selected_vertex = dt->nearest_power_vertex(convert(event->scenePos()));
      dt->remove(selected_vertex);
    }
    Q_EMIT( modelChanged());
  }
}



template <typename T>
bool
RegularTriangulationRemoveVertex<T>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return false;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_REGULAR_TRIANGULATION_REMOVE_VERTEX_H
