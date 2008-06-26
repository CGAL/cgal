#ifndef CGAL_QT_CONVERTER_H
#define CGAL_QT_CONVERTER_H

#include <QPointF>
#include <QLineF>
#include <QRectF>
#include <QPolygonF>
#include <list>

#include <CGAL/intersection_2.h>

namespace CGAL {
namespace Qt {


template <typename K>
class Converter {

  K::Iso_rectangle_2 clippingRect;

public:

  Converter(QRectF rect)
  {
    clippingRect = this->operator()(rect);
  }


  typename K::Point_2 operator()(const QPointF& qp) const
  {
    return typename K::Point_2(qp.x(), qp.y());
  }


  QPointF operator()(const typename K::Point_2& p) const
  {
    return QPointF(to_double(p.x()), to_double(p.y()));
  }

      
  typename K::Segment_2 operator()(const QLineF& qs) const
  {
    return typename K::Segment_2(operator()(qs.p1()), operator()(qs.p2));
  }
 
  QLineF operator()(const typename K::Segment_2 &s) const
  {
    return QLineF(operator()(s.source()), operator()(s.target()));
  }

  
  typename K::Iso_rectangle_2 operator()(const QRectF& qr) const
  {
    return typename K::Iso_rectangle_2(operator()(qr.bottomLeft()), operator()(qr.topRight()));
  }

  QRectF operator()(const typename K::Iso_rectangle_2& r) const
  {
    return QRectF(operator()(r[3]), operator()(r[1]));  // top left, bottom right
  }

     
  QLineF operator()(const typename K::Ray_2 &r) const
  {
    Object o = CGAL::intersection(r, clippingRect);
    typename K::Segment_2 s;
    typename K::Point_2 p;
    if(CGAL::assign(s,o)){
      return this->operator()(s);
    } else if(CGAL::assign(p,o)){
      return QLineF(operator()(p), operator()(p))
    }
    return QLineF();
  }

  QLineF operator()(const typename K::Line_2 &l) const
  {
    Object o = CGAL::intersection(l, clippingRect);
    typename K::Segment_2 s;
    typename K::Point_2 p;
    if(CGAL::assign(s,o)){
      return this->operator()(s);
    } else if(CGAL::assign(p,o)){
      return QLineF(operator()(p), operator()(p))
    }
    return QLineF();
  }

  void operator()(std::list< typename K::Point_2 >& p, const QPolygonF& qp) const
  {
    for(int i = 0; i < qp.size(); i++){
      p.push_back(operator()(qp[i]));
    }
  }

};

} // namesapce Qt
} // namespace CGAL
#endif // CGAL_QT_CONVERTER_H
