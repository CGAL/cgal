
#ifndef CGAL_Q_POLYLINE_INPUT_2_H
#define CGAL_Q_POLYLINE_INPUT_2_H

#include <QPolygonF>
#include <QPointF>

#include <CGAL/IO/QtInput.h>
#include <CGAL/IO/QtConverter.h>

class QGraphicsScene;
class QGraphicsSceneMouseEvent;
class QGraphicsItem;
class QGraphicsPathItem;
class QKeyEvent;
class QEvent;
class QObject;

namespace CGAL {

class QtPolylineInput_non_templated_base : public QtInput
{
public:
  void setNumberOfVertices(int n)
  {
    n_ = n;
  }

protected:
  // protected constructor
  QtPolylineInput_non_templated_base(QObject* parent, 
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
}; // end class QtPolylineInput_non_templated_base

template <typename K>
class QtPolylineInput : public QtPolylineInput_non_templated_base
{
public:
  QtPolylineInput(QObject* parent, QGraphicsScene* s, int n = 0, bool closed = true)
    : QtPolylineInput_non_templated_base(parent, s, n, closed)
  {
  }

protected:
  void generate_polygon() {
    std::list<typename K::Point_2> points;
    QtConverter<K> convert;
    convert(points, this->polygon); 
    emit(generate(CGAL::make_object(points)));
  }
}; // end class QtPolylineInput

} // namespace CGAL

#endif // CGAL_Q_POLYLINE_INPUT_2_H
