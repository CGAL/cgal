// Copyright (c) 2011   INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Triangulation_2/include/CGAL/Triangulation_hyperbolic_traits_2.h $
// $Id: Triangulation_hyperbolic_traits_2.h 57323 2010-07-05 10:07:39Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_PAINTER_OSTREAM_H
#define CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

#include <CGAL/Qt/PainterOstream.h>

namespace CGAL{

namespace Qt {
  
  template <typename K>
  class PainterOstream<Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<K> > : public PainterOstream<K> {
  public:
    typedef PainterOstream<K> Base;
    typedef PainterOstream<Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<K> > Self;
    typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<K> Gt;
    
    typedef typename Gt::Hyperbolic_segment_2      Hyperbolic_segment_2;    
    typedef typename K::Point_2     Point_2;
    
    typedef Line_arc_2<K>           Line_arc_2;
    typedef Circular_arc_2<K>       Circular_arc_2;
    typedef Circular_arc_point_2<K> Circular_arc_point_2;
    
  public:
    PainterOstream(QPainter* p, QRectF rect = QRectF())
    : Base(p, rect), qp(p), convert(rect) {}
    
    using Base::operator <<;
    
    PainterOstream& operator << (Hyperbolic_segment_2 s) {
      if (const Line_arc_2* seg = boost::get<Line_arc_2>(&s)) {
        operator << (*seg);
        return *this;
      }
      
      Circular_arc_2* arc = boost::get<Circular_arc_2>(&s);

      if (arc->squared_radius() > 10 )
        // due to rounding, the arc drawn does not look like it 
        // passes through the endpoints
        // so we replace the arc by a line segment
        {
          Point_2 source(to_double(arc->source().x()),to_double(arc->source().y()));
          Point_2 target(to_double(arc->target().x()),to_double(arc->target().y()));      
          const Line_arc_2 seg(source,target);
          operator << (seg);
          return *this;
        }
          
      operator << (*arc);
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

