
#ifndef CGAL_Q_CIRCLE_INPUT_2_H
#define CGAL_Q_CIRCLE_INPUT_2_H

#include <QGraphicsView>
#include <QRectF>
#include <QPointF>
#include <QGraphicsItem>
#include <QGraphicsEllipseItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>
#include <CGAL/IO/QtConverter.h>

#include <CGAL/IO/QtInput.h>

namespace CGAL {


template <typename K>
class QtCircleInput : public QtInput
{
public:
  QtCircleInput(QGraphicsScene* s, int pointsOnCircle=1); 

protected:
    
  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  

  

private:

  typedef typename K::Point_2 Point_2;
  int m_pointsOnCircle; // 1, 2 or 3
  int count;
  QGraphicsEllipseItem *qcircle;
  QPointF qp, qq, qr;
  Point_2 p, q, r;
  QGraphicsScene *scene_;  
  QtConverter<K> convert;
};


template <typename K>
QtCircleInput<K>::QtCircleInput(QGraphicsScene* s, int pointsOnCircle)
  : qcircle(new QGraphicsEllipseItem()), scene_(s), m_pointsOnCircle(pointsOnCircle), count(0)
{
  qcircle->hide();
  s->addItem(qcircle);
}


template <typename K>
void 
QtCircleInput<K>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{ 
  if(m_pointsOnCircle < 3){
    if(count == 0){
      qp = event->scenePos();
      p = convert(qp);
      count = 1;
    } else {
      qq = event->scenePos();
      if(qp != qq){
	qcircle->hide();
	q = convert(qq);
	if(m_pointsOnCircle == 1){
	  K::FT sd = squared_distance(p,q);
	  emit generate(CGAL::make_object(std::make_pair(p, sd)));
	} else {
	  emit generate(CGAL::make_object(std::make_pair(p, q)));
	}
	count = 0;
      } 
    }
  } else {
    if(count == 0){
      qp = event->scenePos();
      p = convert(qp);
      count = 1;
    } else if(count == 1){
      qq = event->scenePos();
      if(qp != qq){
	q = convert(qq);
	count = 2;
      }
    } else { // count == 2
      qr  = event->scenePos();
      r = convert(qr);
      typename K::Collinear_2 collinear;
      if(! collinear(p,q,r)){
	emit generate(CGAL::make_object(CGAL::Triple<Point_2, Point_2, Point_2>(p,q,r)));
	count = 0;
      }
    }
  }
}


template <typename K>
void 
QtCircleInput<K>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  CGAL::Bbox_2 bb;
  typename K::Construct_circle_2 construct_circle;
  if(count == 0){
    qcircle->hide();
    return;
  } else if(count == 1) {
    qq = event->scenePos();
    q = convert(qq);
    if(qp == qq){
      qcircle->hide();
      return;
    } else {
      if(m_pointsOnCircle == 1){
	K::FT sd = squared_distance(p,q);
	bb = construct_circle(p, sd).bbox();
      } else {
	bb = construct_circle(p, q).bbox();
      }
    }
  } else { // count == 2
    qr = event->scenePos();
    r = convert(qr);
    typename K::Collinear_2 collinear;
    if(collinear(p,q,r)){
      qcircle->hide();
      return;
    } else {
      bb = construct_circle(p, q, r).bbox();
    }
  }
  qcircle->setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
  qcircle->show();
}


template <typename K>
void 
QtCircleInput<K>::keyPressEvent ( QKeyEvent * event ) 
{
  if(event->key() == Qt::Key_Delete){
    if(count>0){
      --count;
    }
  }
}



template <typename K>
bool 
QtCircleInput<K>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_Q_CIRCLE_INPUT_H
