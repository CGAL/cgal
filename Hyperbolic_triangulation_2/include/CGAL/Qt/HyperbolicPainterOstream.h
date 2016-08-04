// Copyright (c) 2011-2016   INRIA Sophia Antipolis, INRIA Nancy (France).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_PAINTER_OSTREAM_H
#define CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

#include <CGAL/Qt/PainterOstream.h>

#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

namespace CGAL{

namespace Qt {
  
  template <typename K>
  class PainterOstream<Hyperbolic_Delaunay_triangulation_traits_2<K> > : public PainterOstream<K> {
  public:
    typedef PainterOstream<K> Base;
    typedef PainterOstream<Hyperbolic_Delaunay_triangulation_traits_2<K> > Self;
    
    typedef Hyperbolic_Delaunay_triangulation_traits_2<K> Gt;
    
  private:
    typedef typename Gt::Segment_2      Segment_2;
    typedef typename Gt::Line_segment_2 Line_segment_2;
    typedef typename Gt::Arc_2          Arc_2;
    //typedef typename Gt::Line_2         Line_2;
    
    typedef typename K::Point_2    Point_2;
    typedef typename K::Circle_2   Circle_2;
    
  public:
    PainterOstream(QPainter* p, QRectF rect = QRectF())
    : Base(p, rect), qp(p), convert(rect) {}
    
    using Base::operator <<;
    
    PainterOstream& operator << (const Segment_2& s)
    {
      if(const Line_segment_2* line_segment = boost::get<Line_segment_2>(&s)){
        operator << (*line_segment);
        return *this;
      }
      if(const Arc_2* arc = boost::get<Arc_2>(&s)){
        
        const Circle_2& circle = get<0>(*arc);
        const Point_2& center = circle.center();
        const Point_2& source = get<1>(*arc);
        const Point_2& target = get<2>(*arc);
        
        if (circle.squared_radius() > 10) {
          const Line_segment_2 seg(source,target);
          operator << (seg);
          return *this;
        }
        
        double asource = std::atan2( -to_double(source.y() - center.y()),
                                    to_double(source.x() - center.x())); 
        double atarget = std::atan2( -to_double(target.y() - center.y()),
                                    to_double(target.x() - center.x()));
        
        std::swap(asource, atarget);
        double aspan = atarget - asource;
        
        if(aspan < 0.)
          aspan += 2 * CGAL_PI;
        
        const double coeff = 180*16/CGAL_PI;
        this->qp->drawArc(this->convert(circle.bbox()), 
                          (int)(asource * coeff), 
                          (int)(aspan * coeff));
      }
      return *this;
    }
    
  private:
    // ToDo: These objects must be deleted
    // Copies of these objects are in the base class.
    // We need access to the copies in the base class.
    QPainter* qp;
    Converter<K> convert;      
  };
  
} //namespace Qt
  
} //namespace CGAL

#endif // CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

