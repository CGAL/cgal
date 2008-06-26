#ifndef CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
#define CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H

#include <QObject>
#include <QPointF>
#include <QString>
#include <QCursor>
#include <QRect>
#include <QRectF>

class QGraphicsView;
class QGraphicsScene;
class QEvent;
class QGraphicsRectItem;

namespace CGAL {
namespace Qt {

class GraphicsViewNavigation: public QObject {

  Q_OBJECT

  signals:
  void mouseCoordinates(QString);

public:
  GraphicsViewNavigation();
  ~GraphicsViewNavigation();
  
  bool eventFilter(QObject *obj, QEvent *event);

private:

  void scaleView(QGraphicsView*, qreal scaleFactor);
  void translateView(QGraphicsView*, int dx,  int dy);
  void drag_to(QGraphicsView*, QPoint new_pos);
  QRectF mapToScene(QGraphicsView*, QRect rect) const;
  void display_parameters(QGraphicsView*);

  QGraphicsRectItem* rectItem;
  QPointF rect_first_point;
  bool dragging;
  QPointF dragging_start;
  QCursor cursor_backup;
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
