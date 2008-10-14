// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_PAINTER_OSTREAM_H
#define CGAL_QT_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <QPainterPath>
#include <CGAL/Qt/Converter.h>

namespace CGAL {

  template <class CK>
  class  Circular_arc_point_2;

  template <class CK>
  class  Circular_arc_2;

  template <class CK>
  class  Line_arc_2;

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
    qp->drawEllipse(convert(c.bbox()));
    return *this;
  }


  PainterOstream& operator<<(const Circular_arc_point_2<K>& p)
  {
    typedef typename K::Point_2   Point_2;
    (*this) << Point_2(to_double(p.x()), to_double(p.y()));
    return *this;
  }


  PainterOstream& operator<<(const Circular_arc_2<K>& arc)
  {
    const typename K::Circle_2 & circ = arc.supporting_circle();
    const typename K::Point_2 & center = circ.center();
    const typename K::Circular_arc_point_2 & source = arc.source();
    const typename K::Circular_arc_point_2 & target = arc.target();

    double asource = std::atan2( -to_double(source.y() - center.y()),
				 to_double(source.x() - center.x())); 
    double atarget = std::atan2( -to_double(target.y() - center.y()),
				 to_double(target.x() - center.x()));

    std::swap(asource, atarget);
    double aspan = atarget - asource;

    if(aspan < 0.)
      aspan += 2 * CGAL_PI;

    const double coeff = 180*16/CGAL_PI;
    qp->drawArc(convert(circ.bbox()), 
		(int)(asource * coeff), 
	         (int)(aspan * coeff));
    return *this;
  }

  PainterOstream& operator<<(const Line_arc_2<K>& arc)
  {
    (*this) << Segment_2<K>(Point_2<K>(to_double(arc.source().x()), to_double(arc.source().y())),
			    Point_2<K>(to_double(arc.target().x()), to_double(arc.target().y())));
     return *this;
  }

  void draw_parabola_segment(const  Point_2<K>& center, const Line_2<K>& line, 
			     const  Point_2<K>& source, const Point_2<K>& target)
  {
    const Point_2<K> proj_source = line.projection(source);
    const Point_2<K> proj_target = line.projection(target);
    const Point_2<K> intersection = circumcenter(proj_source,
						 proj_target,
						 center);
    // Property: "intersection" is the intersection of the two tangent
    // lines in source and target.
    QPainterPath path;
    path.moveTo(convert(source));
    path.quadTo(convert(intersection), convert(target));
    qp->drawPath(path);
  }
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PAINTER_OSTREAM_H
