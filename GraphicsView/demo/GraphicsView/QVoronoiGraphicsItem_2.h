
#ifndef CGAL_Q_VORONOI_GRAPHICS_ITEM_2_H
#define CGAL_Q_VORONOI_GRAPHICS_ITEM_2_H



#include "QGraphicsItem_2.h"
#include "QPainterOstream.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/intersection_2.h>

class QGraphicsSceneMouseEvent;


namespace CGAL {

template <typename DT>
class QVoronoiGraphicsItem_2 : public QGraphicsItem_2
{
public:
  QVoronoiGraphicsItem_2(DT  * dt_);


  QRectF 
  boundingRect() const;
  
  void 
  paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void 
  modelChanged();

private:
  DT * dt;
};



template <typename DT>
QVoronoiGraphicsItem_2<DT>::QVoronoiGraphicsItem_2(DT * dt_)
  :  dt(dt_)
{
  setZValue(3);
}

template <typename DT>
QRectF 
QVoronoiGraphicsItem_2<DT>::boundingRect() const
{
  QRectF rect;
  QList<QGraphicsView *>  views = scene()->views();
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


template <typename DT>
void 
QVoronoiGraphicsItem_2<DT>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *w)
{
  QConverter<typename DT::Geom_traits> convert;
  QRectF rect = boundingRect();
  
  painter->setPen(pen());
  for(typename DT::Finite_edges_iterator eit = dt->finite_edges_begin();
      eit != dt->finite_edges_end();
      eit++){
    CGAL::Object o = dt->dual(eit);
    typename DT::Segment s;
    typename DT::Geom_traits::Ray_2 r;
    typename DT::Geom_traits::Line_2 l;
    if(CGAL::assign(s,o)){
      (*painter) << s;
    } else if(CGAL::assign(r,o)) {
      typename DT::Geom_traits::Iso_rectangle_2 ir;
      ir = convert(rect);
      o = CGAL::intersection(r, ir);
      if(CGAL::assign(s,o)){
	(*painter) << s;
      }
    }else if(CGAL::assign(l,o)) {
      typename DT::Geom_traits::Iso_rectangle_2 ir;
      ir = convert(rect);
      o = CGAL::intersection(l, ir);
      if(CGAL::assign(s,o)){
	(*painter) << s;
      }
    } 
  }
}


template <typename T>
void 
QVoronoiGraphicsItem_2<T>::modelChanged()
{
  update();
}

} // namespace CGAL

#endif // CGAL_Q_VORONOI_GRAPHICS_ITEM_2_H
