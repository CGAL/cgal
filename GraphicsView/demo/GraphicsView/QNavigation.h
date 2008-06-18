#ifndef CGAL_Q_NAVIGATION_H
#define CGAL_Q_NAVIGATION_H

#include <QObject>
#include <QPointF>
#include <QString>
#include <QCursor>
#include <QRect>
#include <QRectF>

class QGraphicsView;
class QEvent;
class QGraphicsRectItem;

namespace CGAL {

class QNavigation: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  QNavigation(QGraphicsView* v_);
  ~QNavigation();
  
  bool eventFilter(QObject *obj, QEvent *event);

private:

  void scaleView(qreal scaleFactor);
  void translateView(int dx,  int dy);
  void drag_to(QPoint new_pos);
  QRectF mapToScene(QRect rect) const;
  void display_parameters();

  QGraphicsView* v;
  QGraphicsRectItem* rectItem;
  QPointF rect_first_point;
  bool dragging;
  QPointF dragging_start;
  QCursor cursor_backup;
};


} // namespace CGAL

#endif // CGAL_Q_NAVIGATION_H
