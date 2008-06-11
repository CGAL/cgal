
#ifndef CGAL_Q_POLYLINE_INPUT_2_H
#define CGAL_Q_POLYLINE_INPUT_2_H

#include <QGraphicsItem>
#include <QGraphicsPathItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>
#include <QPolygonF>
#include <QPainterPath>
#include "QConverter.h"

#include "QInput.h"

namespace CGAL {


template <typename K>
class QPolylineInput_2 : public QInput
{
public:
  QPolylineInput_2(QGraphicsScene* s, int n = 0, bool closed = true);


  void setNumberOfVertices(int n)
  {
    n_ = n;
  }

protected:

  // mousePressEvent returns true iff the event is consummed
  bool mousePressEvent(QGraphicsSceneMouseEvent *event);

  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

  // keyPressEvent returns true iff the event is consummed
  bool keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  
  void rubberbands(const QPointF& p);
  

private:
  bool first;
  QPolygonF polygon;
  QGraphicsPathItem *path_item;
  QGraphicsLineItem *b, *e;
  bool closed_;
  int n_;
  QPointF sp;
  QGraphicsScene *scene_;
};





template <typename K>
QPolylineInput_2<K>::QPolylineInput_2(QGraphicsScene* s, int n, bool closed)
  : path_item(NULL), b(NULL), e(NULL), n_(n), closed_(closed) , scene_(s)
{}


template <typename K>
bool
QPolylineInput_2<K>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{ 
  if( event->modifiers() )
    return false;
  if( event->button() != Qt::RightButton
      && event->button() != Qt::LeftButton )
    return false;
  polygon.push_back(event->scenePos());
  if(path_item){
    scene_->removeItem(path_item);
    delete path_item;
    path_item = NULL;
  }
  if( (n_ == 0 && event->button() == Qt::RightButton)
      || ( event->button() == Qt::LeftButton && (polygon.size() == n_)) ){
    std::list<typename K::Point_2> points;
    QConverter<K> convert;
    convert(points, polygon); 
    emit(generate(CGAL::make_object(points)));
    polygon.clear();
    if(b){
      scene_->removeItem(b);
      delete b;
      b = NULL;
    }
    if(e){
      scene_->removeItem(e);
      delete e;
      e = NULL;
    }
    return true;
  }
  if(event->button() == Qt::LeftButton){
    QPainterPath qpp;
    qpp.addPolygon(polygon);
    path_item = new QGraphicsPathItem(qpp);
    path_item->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    scene_->addItem(path_item);
    return true;
  }
  return false;
}


template <typename K>
void 
QPolylineInput_2<K>::rubberbands(const QPointF& p)
{
  if(polygon.empty()){
    return;
  }
  if(!b && closed_ ){
    b = new QGraphicsLineItem();
    b->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    scene_->addItem(b);
  }
  if( !e){
    e = new QGraphicsLineItem();    
    e->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    scene_->addItem(e);
  }
  if(closed_){
    QLineF bLine(polygon.front(), p);
    b->setLine(bLine);
  }
  QLineF eLine(polygon.back(), p);
  e->setLine(eLine); 
}


template <typename K>
void 
QPolylineInput_2<K>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  sp = event->scenePos();
  rubberbands(sp);
}


template <typename K>
bool
QPolylineInput_2<K>::keyPressEvent ( QKeyEvent * event ) 
{
  if( event->modifiers() )
    return false;

  switch(event->key())
  {
  case Qt::Key_Delete:
  case Qt::Key_Escape:
  case Qt::Key_Backspace:
    break;
  default:
    return false;
  }
  if(polygon.empty()){
    return true;
  }
  polygon.pop_back();
  if(polygon.empty()){
    if(b){
      scene_->removeItem(b);
      delete b;
      b = NULL;
    }
    if(e){
      scene_->removeItem(e);
      delete e;
      e = NULL;
    }
    return true;
  }
  if(path_item){
    scene_->removeItem(path_item);
    delete path_item;
    path_item = NULL;
  }
  QPainterPath qpp;
  qpp.addPolygon(polygon);
  path_item = new QGraphicsPathItem(qpp);
  path_item->setPen(QPen(Qt::red, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene_->addItem(path_item);
  rubberbands(sp);
  return true;
}



template <typename K>
bool 
QPolylineInput_2<K>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    if(!mousePressEvent(mouseEvent)) {
      return QObject::eventFilter(obj, event);
    }
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return QObject::eventFilter(obj, event);
  } else if (event->type() == QEvent::KeyPress) {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
    if(!keyPressEvent(keyEvent)) {
      return QObject::eventFilter(obj, event);
    }
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
} 

} // namespace CGAL

#endif // CGAL_Q_POLYLINE_INPUT_H
