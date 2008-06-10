#ifndef CGAL_Q_NAVIGATION_H
#define CGAL_Q_NAVIGATION_H

#include <QObject>
#include <QEvent>
#include <QMouseEvent>
#include <QPointF>
#include <QString>
#include <QWheelEvent>
#include <QGraphicsView>

namespace CGAL {

class QNavigation: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  QNavigation(QGraphicsView* v_)
    : v(v_)
  {}
  
  bool eventFilter(QObject *obj, QEvent *event);

private:

  void scaleView(qreal scaleFactor);


  QGraphicsView* v;

};


} // namespace CGAL

#endif // CGAL_Q_NAVIGATION_H
