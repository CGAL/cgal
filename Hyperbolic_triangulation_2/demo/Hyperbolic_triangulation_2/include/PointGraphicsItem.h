#ifndef CGAL_QT_POINT_GRAPHICS_ITEM_H
#define CGAL_QT_POINT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename T>
class PointGraphicsItem : public GraphicsItem
{
public:
  typedef typename T::Point Point;
  typedef typename T::Geom_traits Geom_traits;

  PointGraphicsItem(Point p = Point());
  
  void modelChanged();
  
public:
  
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void setPoint(const Point& p)
  {
    m_p = p;
    
    updateBoundingBox();
    update();
  }
  
  Point getPoint() const
  {
    return m_p;
  }
  
  const QPen& vertexPen() const
  {
    return vertex_pen;
  }

  void setVertexPen(const QPen& pen)
  {
    vertex_pen = pen;
  }
    
  void paintOneVertex(const Point& point) const;

protected:
  void updateBoundingBox();
  void paintVertex(QPainter* painter, const Point& point) const;
  
  QPainter* m_painter;

  Point m_p;
  
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;
  
  QPen vertex_pen;
};


template <typename T>
PointGraphicsItem<T>::PointGraphicsItem(Point p)
: m_painter(0), m_p(p), bb(0,0,0,0), bb_initialized(false)
{
  setVertexPen(QPen(::Qt::green, 3.));
  
  updateBoundingBox();
  setZValue(3);
}
  

template <typename T>
QRectF 
PointGraphicsItem<T>::boundingRect() const
{
  return bounding_rect;
}
  
  
template <typename T>
void 
PointGraphicsItem<T>::modelChanged()
{
  updateBoundingBox();
  update();
}

  
template <typename T>
void 
PointGraphicsItem<T>::paintOneVertex(const Point& p) const
{
  pointVertex(m_painter, p);
}
  
  
template <typename T>
void 
PointGraphicsItem<T>::paintVertex(QPainter *painter, const Point& p) const
{
  Converter<Geom_traits> convert;
  
  painter->setPen(this->vertexPen());
  QMatrix matrix = painter->matrix();
  painter->resetMatrix();
  painter->drawPoint(matrix.map(convert(p)));
  painter->setMatrix(matrix);
}

template <typename T>
void 
PointGraphicsItem<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                            QWidget * /*widget*/)
{
  painter->setPen(this->vertexPen());
    
  paintVertex(painter, m_p);
}
  

template <typename T>
void 
PointGraphicsItem<T>::updateBoundingBox()
{
  bb = m_p.bbox();
  bb_initialized = true;
 
  bounding_rect = QRectF(bb.xmin(),
                         bb.ymin(),
                         bb.xmax()-bb.xmin(),
                         bb.ymax()-bb.ymin());
}
  
} // namespace Qt
  
} // namespace CGAL
  
#endif // CGAL_QT_POINT_GRAPHICS_ITEM_H
  
