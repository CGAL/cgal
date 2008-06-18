#ifndef CGAL_Q_PAINTER_OSTREAM_H
#define CGAL_Q_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <CGAL/IO/QtConverter.h>

namespace CGAL {

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

} // namespace CGAL

#endif // CGAL_Q_PAINTER_OSTREAM_H
