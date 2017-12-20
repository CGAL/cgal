// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_2_H
#define CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>

//--------------------------------------------------------------------

namespace CGAL {

namespace ApolloniusGraph_2 {

template<class K, class MTag>
class Orientation_2
{
public:
  typedef K                        Kernel;
  typedef MTag                     Method_tag;
  typedef typename K::Site_2       Site_2;
  typedef typename K::Point_2      Point_2;
  typedef typename K::Orientation  Orientation;

  typedef Orientation              result_type;
  typedef Site_2                   argument_type;

private:
  typedef Weighted_point_inverter_2<K>   Weighted_point_inverter;
  typedef Inverted_weighted_point_2<K>   Inverted_weighted_point;
  typedef Voronoi_circle_2<K>            Voronoi_circle;
  typedef Bitangent_line_2<K>            Bitangent_line;
  typedef typename Bitangent_line::FT    FT;

private:
  Orientation vv_orientation(const Voronoi_circle& vc, const Point_2& sp1,
			     const Point_2& p1,	const Point_2& p2,
			     const Field_with_sqrt_tag&) const
  {
    FT a = vc.a1() + vc.a2() * CGAL::sqrt(vc.delta());
    FT b = vc.b1() + vc.b2() * CGAL::sqrt(vc.delta());
    FT det1 = a * (p2.y() - p1.y()) - b * (p2.x() - p1.x());
    FT c = vc.c1() + vc.c2() * CGAL::sqrt(vc.delta());
    FT det2 = determinant(p1.x() - sp1.x(), p1.y() - sp1.y(),
				p2.x() - sp1.x(), p2.y() - sp1.y());
    return CGAL::sign(det1 + FT(2) * c * det2);
  }

  Orientation vv_orientation(const Voronoi_circle vc, const Point_2& sp1,
			     const Point_2& p1, const Point_2& p2,
			     const Integral_domain_without_division_tag&) const
  {
    FT dx = p2.x() - p1.x();
    FT dy = p2.y() - p1.y();
    FT det1 = determinant(p1.x() - sp1.x(), p1.y() - sp1.y(),
				p2.x() - sp1.x(), p2.y() - sp1.y());
    FT A = vc.a1() * dy - vc.b1() * dx + FT(2) * vc.c1() * det1;
    FT B = vc.a2() * dy - vc.b2() * dx + FT(2) * vc.c2() * det1;
    return sign_a_plus_b_x_sqrt_c(A, B, vc.delta());
  }

public:
  inline
  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3) const
  {
    return Kernel().orientation_2_object()(s1.point(), s2.point(),
					   s3.point());
  }

  Orientation operator()(const Site_2& s1, const Site_2& s2,
			 const Site_2& s3, const Site_2& p1,
			 const Site_2& p2) const
  {
    // computes the operation of the Voronoi vertex of s1, s2, s3 and
    // the points p1 and p2
    Weighted_point_inverter inverter(s1);
    Inverted_weighted_point u2 = inverter(s2);
    Inverted_weighted_point u3 = inverter(s3);
    Bitangent_line blinv_23(u2, u3);
    Voronoi_circle vc(blinv_23);
    return
      vv_orientation(vc, s1.point(), p1.point(), p2.point(), Method_tag());
  }
};

//--------------------------------------------------------------------

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_ORIENTATION_2_H
