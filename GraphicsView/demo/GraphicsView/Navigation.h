#ifndef CGAL_NAVIGATION_H
#define CGAL_NAVIGATION_H

#include <QObject>
#include <QEvent>
#include <QMouseEvent>
#include <QPointF>
#include <QString>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QGraphicsView>

namespace CGAL {

class Navigation: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  Navigation(QGraphicsView* v_)
    : v(v_)
  {}
  
  bool eventFilter(QObject *obj, QEvent *event);

private:

  void keyPressEvent(QKeyEvent *event);

  void scaleView(qreal scaleFactor);


  QGraphicsView* v;

};


} // namespace CGAL

#endif // CGAL_NAVIGATION_H
