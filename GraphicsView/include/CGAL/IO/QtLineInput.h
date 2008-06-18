
#ifndef CGAL_Q_LINE_INPUT_2_H
#define CGAL_Q_LINE_INPUT_2_H

#include <QGraphicsView>
#include <QRectF>
#include <QPointF>
#include <QGraphicsItem>
#include <QGraphicsLineItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>
#include <CGAL/IO/QtConverter.h>

#include <CGAL/IO/QtInput.h>

namespace CGAL {


template <typename K>
class QtLineInput : public QtInput
{
public:
  QtLineInput(QGraphicsScene* s);

protected:
    
  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  

  

private:

  QRectF boundingRect();
  QLineF qlinef();

  bool second;
  QGraphicsLineItem line;
  QPointF qsp, qtp;
  typename K::Point_2 sp, tp;
  typename K::Line_2 l;
  QGraphicsScene *scene_;  
  QtConverter<K> convert;
};

template <typename K>
QRectF
LineInput_2<K>::boundingRect()
{
  QRectF rect;
  QList<QGraphicsView *>  views = scene_->views();
  for (int i = 0; i < views.size(); ++i) {
    QGraphicsView *view = views.at(i);
    QRect vprect = view->viewport()->rect();
    QPoint tl = vprect.topLeft();
    QPoint br = vprect.bottomRight();
    QPointF tlf = view->mapToScene(tl);
    QPointF brf = view->mapToScene(br);
    rect = QRectF(tlf, brf);
  }
  return rect;
}


template <typename K>
QtLineInput<K>::QtLineInput(QGraphicsScene* s)
  : second(false), scene_(s)
{}


template <typename K>
QLineF
QtLineInput<K>::qlinef()
{

  sp = convert(qsp);
  tp = convert(qtp);
  typename K::Line_2  l(sp,tp);
  QRectF qrect(boundingRect());
  typename K::Iso_rectangle_2 rect;
  rect = convert(qrect);
  Object o = intersection(l,rect);
  typename K::Segment_2 s;
  assign(s, o);
  return convert(s);
}

template <typename K>
void 
QtLineInput<K>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{ 
  if(second){
      qtp = event->scenePos();
      sp = convert(qsp);
      tp = convert(qtp);
      scene_->removeItem(&line);
      emit generate(CGAL::make_object(typename K::Line_2(sp,tp)));
  } else {
    qsp = event->scenePos();
    qtp = QPointF(qsp.x()+1, qsp.y());
    line.setLine(qlinef());
    scene_->addItem(&line);
  }
  second = !second;
}




template <typename K>
void 
QtLineInput<K>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  qtp = event->scenePos();
  if(qtp == qsp){
    qtp = QPointF(qsp.x()+1, qsp.y());
  } 
  line.setLine(qlinef());
}


template <typename K>
void 
QtLineInput<K>::keyPressEvent ( QKeyEvent * event ) 
{
  if(event->key() != Qt::Key_Delete){
    return;
  }
  if(second){
    scene_->removeItem(&line);
    second = false;
  }
}



template <typename K>
bool 
QtLineInput<K>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::KeyPress) {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
    keyPressEvent(keyEvent);
    return true;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
} 

} // namespace CGAL

#endif // CGAL_Q_LINE_INPUT_H
