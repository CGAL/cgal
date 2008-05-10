
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
      return QObject::eventFilter(obj, event);
    } else if (event->type() == QEvent::KeyPress) {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      keyPressEvent(keyEvent);
      return true;
    } else{
      // standard event processing
      return QObject::eventFilter(obj, event);
    }
  }



  void 
Navigation::keyPressEvent(QKeyEvent *event)
  {
          std::cout << "in keyPressEvent" << std::endl;
    /*
    QRectF rect;
    QRect vprect = this->viewport()->rect();
    QPoint tl = vprect.topLeft();
    QPoint br = vprect.bottomRight();
    QPointF tlf = this->mapToScene(tl);
    QPointF brf = this->mapToScene(br);
    rect = QRectF(tlf, brf);
    */

    switch (event->key()) {
    case Qt::Key_Up:
      std::cout << "up" << std::endl;

      v->translate(0, -10);  // rect.height()
      break;
    case Qt::Key_Down:
      v->translate(0, 10); // rect.height()
      break;
    case Qt::Key_Left:
      v->translate(10 ,0);  // -rect.width()
      break;
    case Qt::Key_Right:
      v->translate(10,0);// rect.width()
      break;
    case Qt::Key_Plus:
      scaleView(1.2);
      break;
    case Qt::Key_Minus:
      scaleView(1 / 1.2);
      break;
    case Qt::Key_Space:
    case Qt::Key_Enter:

      break;
    default:
      /// what should be the default??? 
      break;
    }
  }

  void 
  Navigation::scaleView(qreal scaleFactor)
  {
    qreal factor = v->matrix().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    if (factor < 0.001 || factor > 2000)
        return;

    v->scale(scaleFactor, scaleFactor);
}

#include "Navigation.moc"

} // namespace CGAL

