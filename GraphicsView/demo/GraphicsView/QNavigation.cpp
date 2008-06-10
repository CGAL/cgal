
#include "Navigation.h"
#include <cmath>
#include <iostream>

namespace CGAL {

  bool 
  Navigation::eventFilter(QObject *obj, QEvent *event)
  {
    if (event->type() == QEvent::Wheel) {
      QWheelEvent *wheelEvent = static_cast<QWheelEvent*>(event);
      scaleView(pow((double)2, -wheelEvent->delta() / 240.0));
      return true;
    } else if (event->type() == QEvent::MouseMove) {
      QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
      QPointF pos = v->mapToScene(mouseEvent->pos());
      QString xy = QString(" ") + QString::number(pos.x(),'g', 6) + " , " + QString::number(pos.y(),'g', 6) + " ";
      emit mouseCoordinates(xy);
      return false; // QObject::eventFilter(obj, event);
    } else{
      // standard event processing
      return false; // QObject::eventFilter(obj, event);
    }
  }


  void 
  Navigation::scaleView(qreal scaleFactor)
  {
    qreal factor = v->matrix().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    //if (factor < 0.001 || factor > 2000)
    //    return;

    v->scale(scaleFactor, scaleFactor);
}

#include "Navigation.moc"

} // namespace CGAL

