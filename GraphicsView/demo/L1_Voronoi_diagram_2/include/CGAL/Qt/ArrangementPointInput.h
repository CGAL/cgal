
#ifndef CGAL_QT_ARRANGEMENT_POINT_INPUT_H
#define CGAL_QT_ARRANGEMENT_POINT_INPUT_H

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>

namespace CGAL {
namespace Qt {

template <typename Arrangement>
class ArrangementPointInput : public GraphicsViewInput
{
public:
  typedef typename Arrangement::Geometry_traits_2::Kernel Kernel;
  typedef typename Kernel::Point_2                        Point_2;

  ArrangementPointInput(QObject* parent);
  
protected:
  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  Converter<Kernel> m_convert;
  Point_2 m_p;
};


template <typename T>
ArrangementPointInput<T>::ArrangementPointInput(QObject* parent)
  :  GraphicsViewInput(parent)
{}

template <typename T>
void 
ArrangementPointInput<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  m_p = m_convert(event->scenePos());
}


template <typename T>
void 
ArrangementPointInput<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent* )
{
  emit (generate(CGAL::make_object(m_p)));
}



template <typename T>
bool 
ArrangementPointInput<T>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_QT_ARRANGEMENT_POINT_INPUT_H
