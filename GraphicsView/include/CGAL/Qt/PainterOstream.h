// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

namespace internal {

  template < class K >
  struct Is_circular_kernel : public has_Circular_arc_point_2< K >
  { };

  template < class K, bool b = Is_circular_kernel< K >::value >
  struct Circular_kernel_2_types
  {
      struct Circular_arc_point_2 { };
      struct Circular_arc_2 { };
      struct Line_arc_2 { };
  };

  template < class K >
  struct Circular_kernel_2_types< K, true >
  {
      typedef typename K::Circular_arc_point_2 Circular_arc_point_2;
      typedef typename K::Circular_arc_2 Circular_arc_2;
      typedef typename K::Line_arc_2 Line_arc_2;
  };
}

namespace Qt {

template <typename K>
class PainterOstream {

private:
  QPainter* qp;
  Converter<K> convert;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Ray_2 Ray_2;
  typedef typename K::Line_2 Line_2;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename K::Circle_2 Circle_2;
  typedef typename internal::Circular_kernel_2_types<K>::Circular_arc_point_2 
    Circular_arc_point_2;
  typedef typename internal::Circular_kernel_2_types<K>::Circular_arc_2 
    Circular_arc_2;
  typedef typename internal::Circular_kernel_2_types<K>::Line_arc_2 
    Line_arc_2;
  
public:
  PainterOstream(QPainter* p, QRectF rect = QRectF())
    : qp(p), convert(rect)
  {}

  PainterOstream& operator<<(const Point_2& p)
  {
    qp->drawPoint(convert(p));
    return *this;
  }
  
  PainterOstream& operator<<(const Segment_2& s)
  {
    qp->drawLine(convert(s.source()), convert(s.target()));
    return *this;
  }
  
  
  PainterOstream& operator<<(const Ray_2& r)
  {
    qp->drawLine(convert(r));
    return *this;
  }

  
  PainterOstream& operator<<(const Line_2& l)
  {
    qp->drawLine(convert(l));
    return *this;
  }


  PainterOstream& operator<<(const Triangle_2& t)
  {
    qp->drawPolygon(convert(t));
    return *this;
  }

  PainterOstream& operator<<(const Iso_rectangle_2& r)
  {
    qp->drawRect(convert(r));
    return *this;
  }

  PainterOstream& operator<<(const Bbox_2& bb)
  {
    qp->drawRect(convert(bb));
    return *this;
  }

  PainterOstream& operator<<(const Circle_2& c)
  {
    qp->drawEllipse(convert(c.bbox()));
    return *this;
  }


  PainterOstream& operator<<(const Circular_arc_point_2& p)
  {
    typedef typename K::Point_2   Point_2;
    (*this) << Point_2(to_double(p.x()), to_double(p.y()));
    return *this;
  }


  PainterOstream& operator<<(const Circular_arc_2& arc)
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

  PainterOstream& operator<<(const Line_arc_2& arc)
  {
    (*this) << Segment_2(Point_2(to_double(arc.source().x()), to_double(arc.source().y())),
			 Point_2(to_double(arc.target().x()), to_double(arc.target().y())));
     return *this;
  }

  void draw_parabola_segment(const  Point_2& center, const Line_2& line, 
			     const  Point_2& source, const Point_2& target)
  {
    if (CGAL::collinear(source,target,center))
      qp->drawLine(convert(source), convert(target));      
    else
    {
      const Point_2 proj_source = line.projection(source);
      const Point_2 proj_target = line.projection(target);      
      const Point_2 intersection = circumcenter(proj_source,
                                                proj_target,
                                                center);
      // Property: "intersection" is the intersection of the two tangent
      // lines in source and target.
      QPainterPath path;
      path.moveTo(convert(source));
      path.quadTo(convert(intersection), convert(target));
      qp->drawPath(path);
    }
  }
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PAINTER_OSTREAM_H
