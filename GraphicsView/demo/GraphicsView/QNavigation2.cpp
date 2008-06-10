
#include "Navigation2.h"
#include <cmath>
#include <iostream>

namespace CGAL {

  bool 
  Navigation2::eventFilter(QObject *obj, QEvent *event)
  {
    if (event->type() == QEvent::KeyPress) {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      return keyPressEvent(keyEvent);
    } else{
      // standard event processing
      return false; //return QObject::eventFilter(obj, event);
    }
  }



  bool
  Navigation2::keyPressEvent(QKeyEvent *event)
  {    
    std::cout << "Navigation2::keyPressEvent" << std::endl;
    QRectF rect;
    QRect vprect = v->viewport()->rect();
    QPoint tl = vprect.topLeft();
    QPoint br = vprect.bottomRight();
    QPointF tlf = v->mapToScene(tl);
    QPointF brf = v->mapToScene(br);
    rect = QRectF(tlf, brf);
    
    switch (event->key()) {
    case Qt::Key_Up:
      v->translate(0, -rect.height()/2);
      break;
    case Qt::Key_Down:
      v->translate(0, rect.height()/2);
      break;
    case Qt::Key_Left:
      v->translate(-rect.width()/2 ,0);
      break;
    case Qt::Key_Right:
      v->translate(rect.width()/2,0);
      break;
    case Qt::Key_Plus:
      scaleView(1.2);
      break;
    case Qt::Key_Minus:
      scaleView(1 / 1.2);
      break;
    default:
      return false;
    }
    return true;
  }


  void 
  Navigation2::scaleView(qreal scaleFactor)
  {
    v->scale(scaleFactor, scaleFactor);
  }

#include "Navigation2.moc"

} // namespace CGAL

