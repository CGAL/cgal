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

  bool clippingRectIsInitialized;
  typename K::Iso_rectangle_2 clippingRect;

public:

  Converter()
    : clippingRectIsInitialized(false)
  {
  }

  Converter(QRectF rect)
  {
    if(rect.isValid()) {
      clippingRect = this->operator()(rect);
      clippingRectIsInitialized = true;
    }
    else
      clippingRectIsInitialized = false;
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


  QRectF operator()(const CGAL::Bbox_2& bb)
  {
    return QRectF(bb.xmin(),
		  bb.ymin(),
		  bb.xmax()-bb.xmin(),
		  bb.ymax()-bb.ymin());
  }

     
  QLineF operator()(const typename K::Ray_2 &r) const
  {
    CGAL_assertion(clippingRectIsInitialized);
    Object o = CGAL::intersection(r, clippingRect);
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Point_2 Point_2;
    if(const Segment_2 *s = CGAL::object_cast<Segment_2>(&o)){
      return this->operator()(*s);
    } else if(const Point_2 *p = CGAL::object_cast<Point_2>(&o)){
      return QLineF(operator()(*p), operator()(*p));
    }
    return QLineF();
  }

  QLineF operator()(const typename K::Line_2 &l) const
  {
    CGAL_assertion(clippingRectIsInitialized);
    Object o = CGAL::intersection(l, clippingRect);
    typedef typename K::Segment_2 Segment_2;
    typedef typename K::Point_2 Point_2;
    if(const Segment_2 *s = CGAL::object_cast<Segment_2>(&o)){
      return this->operator()(*s);
    } else if(const Point_2 *p = CGAL::object_cast<Point_2>(&o)){
      return QLineF(operator()(*p), operator()(*p));
    }
    return QLineF();
  }

  QPolygonF operator()(const typename K::Triangle_2 &t)
  {
    QPolygonF qp;
    qp << operator()(t.vertex(0)) << operator()(t.vertex(1)) << operator()(t.vertex(2));
    return qp;
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
