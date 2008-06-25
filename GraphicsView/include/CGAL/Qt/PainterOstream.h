#ifndef CGAL_QT_PAINTER_OSTREAM_H
#define CGAL_QT_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <CGAL/Qt/Converter.h>

namespace CGAL {
namesapce Qt {

template <typename K>
QPainter& operator<<(QPainter& qp, const Point_2<K>& p)
{
  QtConverter<K> convert;
  qp.drawPoint(convert(p));
  return qp;
}

template <typename K>
QPainter& operator<<(QPainter& qp, const Segment_2<K>& s)
{
  QtConverter<K> convert;
  qp.drawLine(convert(s.source()), convert(s.target()));
  return qp;
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PAINTER_OSTREAM_H
