
#ifndef CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H

#include <QPolygonF>
#include <QPointF>

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>

class QGraphicsScene;
class QGraphicsSceneMouseEvent;
class QGraphicsItem;
class QGraphicsPathItem;
class QKeyEvent;
class QEvent;
class QObject;

namespace CGAL {

class GraphicsViewPolylineInput_non_templated_base : public GraphicsViewInput
{
public:
  void setNumberOfVertices(int n)
  {
    n_ = n;
  }

protected:
  // protected constructor
  GraphicsViewPolylineInput_non_templated_base(QObject* parent, 
                                     QGraphicsScene* s,
                                     int n = 0,
                                     bool closed = true);


  // mousePressEvent returns true iff the event is consummed
  bool mousePressEvent(QGraphicsSceneMouseEvent *event);

  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

  // keyPressEvent returns true iff the event is consummed
  bool keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  
  void rubberbands(const QPointF& p);

  virtual void generate_polygon() = 0;

protected:
  QPolygonF polygon;

private:
  QGraphicsPathItem *path_item;
  QGraphicsLineItem *b, *e;
  bool closed_;
  int n_;
  QPointF sp;
  QGraphicsScene *scene_;
}; // end class GraphicsViewPolylineInput_non_templated_base

template <typename K>
class GraphicsViewPolylineInput : public GraphicsViewPolylineInput_non_templated_base
{
public:
  GraphicsViewPolylineInput(QObject* parent, QGraphicsScene* s, int n = 0, bool closed = true)
    : GraphicsViewPolylineInput_non_templated_base(parent, s, n, closed)
  {
  }

protected:
  void generate_polygon() {
    std::list<typename K::Point_2> points;
    Converter<K> convert;
    convert(points, this->polygon); 
    emit(generate(CGAL::make_object(points)));
  }
}; // end class GraphicsViewPolylineInput

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H
