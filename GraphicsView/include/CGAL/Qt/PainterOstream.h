#ifndef CGAL_QT_PAINTER_OSTREAM_H
#define CGAL_QT_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <CGAL/Qt/Converter.h>

namespace CGAL {
namespace Qt {

template <typename K>
class PainterOstream {

private:
  QPainter* qp;
  Converter<K> convert;
  
public:
  PainterOstream(QPainter* p, QRectF rect = QRectF())
    : qp(p), convert(rect)
  {}

  PainterOstream& operator<<(const Point_2<K>& p)
  {
    qp->drawPoint(convert(p));
    return *this;
  }
  
  PainterOstream& operator<<(const Segment_2<K>& s)
  {
    qp->drawLine(convert(s.source()), convert(s.target()));
    return *this;
  }
  
  
  PainterOstream& operator<<(const Ray_2<K>& r)
  {
    qp->drawLine(convert(r));
    return *this;
  }

  
  PainterOstream& operator<<(const Line_2<K>& l)
  {
    qp->drawLine(convert(l));
    return *this;
  }


  PainterOstream& operator<<(const Triangle_2<K>& t)
  {
    qp->drawPolygon(convert(t));
    return *this;
  }

  PainterOstream& operator<<(const Iso_rectangle_2<K>& r)
  {
    qp->drawRect(convert(r));
    return *this;
  }

  PainterOstream& operator<<(const Circle_2<K>& c)
  {
    qp->drawRect(convert(c.bbox()));
    return *this;
  }
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PAINTER_OSTREAM_H
