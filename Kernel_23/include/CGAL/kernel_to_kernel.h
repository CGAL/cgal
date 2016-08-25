// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_KERNEL_TO_KERNEL_H
#define CGAL_KERNEL_TO_KERNEL_H

#include <CGAL/config.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/LEDA_basic.h>
#include <CGAL/leda_integer.h>
#include <LEDA/geo/rat_point.h>
#endif

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>

namespace CGAL {

template <class NumberType>
struct Cartesian_double_to_Homogeneous
{
  typedef Point_2< Homogeneous< NumberType> >    Point2;
  typedef Segment_2< Homogeneous< NumberType> >  Segment;

  Cartesian_double_to_Homogeneous() {}

  Point2
  operator()(  const Point_2<Cartesian<double> >& p) const
  { return Point2( NumberType(p.x()), NumberType(p.y()) ); }

  Segment
  operator()(  const Segment_2<Cartesian<double> >& s) const
  {
    return Segment( Point2( NumberType(s.source().x()),
                            NumberType(s.source().y()) ),
                    Point2( NumberType(s.target().x()),
                            NumberType(s.target().y()) ) );
  }
};

#ifdef CGAL_USE_LEDA
struct Cartesian_double_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point2;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< double> >& s) const
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return Segment(
      Point2(CGAL_LEDA_SCOPE::to_double(rs.X()),
             CGAL_LEDA_SCOPE::to_double(rs.Y()),
             CGAL_LEDA_SCOPE::to_double(rs.W())),
      Point2(CGAL_LEDA_SCOPE::to_double(rt.X()),
             CGAL_LEDA_SCOPE::to_double(rt.Y()),
             CGAL_LEDA_SCOPE::to_double(rt.W())) );
  }
};

struct Cartesian_float_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point2;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< float> >& s) const
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return Segment(
      Point2(CGAL_LEDA_SCOPE::to_double(rs.X()),
             CGAL_LEDA_SCOPE::to_double(rs.Y()),
             CGAL_LEDA_SCOPE::to_double(rs.W())),
      Point2(CGAL_LEDA_SCOPE::to_double(rt.X()),
             CGAL_LEDA_SCOPE::to_double(rt.Y()),
             CGAL_LEDA_SCOPE::to_double(rt.W())) );
  }
};
#endif // CGAL_USE_LEDA

} //namespace CGAL

#endif // CGAL_KERNEL_TO_KERNEL_H
