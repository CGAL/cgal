#ifndef CGAL_Q_NAVIGATION2_H
#define CGAL_Q_NAVIGATION2_H

#include <QObject>
#include <QEvent>
#include <QMouseEvent>
#include <QPointF>
#include <QString>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QGraphicsView>

namespace CGAL {

class QNavigation2: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  QNavigation2(QGraphicsView* v_)
    : v(v_)
  {}
  
  bool eventFilter(QObject *obj, QEvent *event);

private:

  bool keyPressEvent(QKeyEvent *event);

  void scaleView(qreal scaleFactor);


  QGraphicsView* v;

};


} // namespace CGAL

#endif // CGAL_Q_NAVIGATION2_H
