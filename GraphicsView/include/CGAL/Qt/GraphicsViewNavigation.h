#ifndef CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
#define CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H

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
namespace Qt {

class GraphicsViewNavigation: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  GraphicsViewNavigation(QGraphicsView* v_);
  ~GraphicsViewNavigation();
  
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

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
