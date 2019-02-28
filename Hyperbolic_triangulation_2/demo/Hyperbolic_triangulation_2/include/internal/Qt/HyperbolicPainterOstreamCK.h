// Copyright (c) 2011-2016 INRIA Sophia Antipolis, INRIA Nancy (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     :  Monique Teillaud <Monique.Teillaud@inria.fr>
//                  Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_PAINTER_OSTREAM_H
#define CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

#include <CGAL/Qt/PainterOstream.h>

namespace CGAL{

namespace Qt {
  
  template <typename K>
  class PainterOstream<Hyperbolic_Delaunay_triangulation_CK_traits_2<K> > 
    : public PainterOstream<K> 
  {
  public:
    typedef PainterOstream<K> Base;
    
    typedef Hyperbolic_Delaunay_triangulation_CK_traits_2<K> Gt;
    typedef PainterOstream<Gt> Self;
        
    typedef typename Gt::Hyperbolic_segment_2      Hyperbolic_segment_2;

    typedef typename Gt::Point_2    Point_2;
    typedef Line_arc_2<K>           Line_arc;
    typedef Circular_arc_2<K>       Circular_arc;
    typedef Circular_arc_point_2<K> Circular_arc_point;
    
    PainterOstream(QPainter* p, QRectF rect = QRectF())
      : Base(p, rect), qp(p), convert(rect) 
    {}
    
    using Base::operator <<;
    
    PainterOstream& operator << (Hyperbolic_segment_2 s)
      {
	if(const Line_arc* seg = boost::get<Line_arc>(&s)) {
	  operator << (*seg);
	  return *this;
	}
      
	Circular_arc* arc = boost::get<Circular_arc>(&s);

	if(arc->squared_radius() > 10)
	  // due to rounding, the arc drawn does not look like it 
	  // passes through the endpoints
	  // so we replace the arc by a line segment
	  {
	    Point_2 source(to_double(arc->source().x()),to_double(arc->source().y()));
	    Point_2 target(to_double(arc->target().x()),to_double(arc->target().y()));	    
	    const Line_arc seg(source,target);
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

