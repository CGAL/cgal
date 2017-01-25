// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_ORIENTED_SIDE_OF_BISECTOR_C2_H
#define CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_ORIENTED_SIDE_OF_BISECTOR_C2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/enum.h>
#include <CGAL/Uncertain.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/Apollonius_graph_2/uncertain/uncertain_functions_on_signs.h>

namespace CGAL {

//--------------------------------------------------------------------

template<class K, class MTag>
class Ag2_uncertain_oriented_side_of_bisector_C2
{
public:
  typedef K                               Kernel;
  typedef MTag                            Method_tag;

  typedef typename K::Point_2             Point_2;
  typedef typename K::Site_2              Site_2;

private:
  typedef typename Kernel::RT             RT;

private:
  Uncertain<Comparison_result>
  compare_distances(const Site_2& p1, const Site_2& p2,
		    const Point_2& p, const Integral_domain_without_division_tag&) const
  {
#ifdef AG2_PROFILE_PREDICATES
    ag2_predicate_profiler::side_of_bisector_counter++;
#endif

    // this function compares the distances of the point(x, y) from the 
    // disks {(x1, y1), w1} and {(x2, y2), w2}
    RT D1 = CGAL::square(p1.x() - p.x()) + CGAL::square(p1.y() - p.y());
    RT D2 = CGAL::square(p2.x() - p.x()) + CGAL::square(p2.y() - p.y());
    RT Dw = p2.weight() - p1.weight();

    Uncertain<Sign> sign_of_Dw = CGAL::sign(Dw);
    if ( is_indeterminate(sign_of_Dw) ) {
      return Uncertain<Comparison_result>::indeterminate();
    }

    Uncertain<Comparison_result> R = CGAL::compare(D1, D2);
    if ( is_indeterminate(R) ) {
      return Uncertain<Comparison_result>::indeterminate();
    }

    if ( sign_of_Dw == ZERO ) {
      return R;
    }
    if ( sign_of_Dw == POSITIVE ) {
      if ( R != SMALLER )  return LARGER;

      return uncertain_sign_a_plus_b_x_sqrt_c(D1 - D2 + CGAL::square(Dw),
					      RT(2) * Dw, D1);
    }

    if ( R != LARGER )  return SMALLER;
    return uncertain_sign_a_plus_b_x_sqrt_c(D1 - D2 - CGAL::square(Dw),
				            RT(2) * Dw, D2);
  }

  Comparison_result
  compare_distances(const Site_2& p1, const Site_2& p2,
		    const Point_2 &p, const Field_with_sqrt_tag&) const
  {
#ifdef AG2_PROFILE_PREDICATES
    ag2_predicate_profiler::side_of_bisector_counter++;
#endif
    // this function compares the distances of the point(x, y) from the 
    // disks {(x1, y1), w1} and {(x2, y2), w2}

    RT D1 = CGAL::square(p1.x() - p.x()) + CGAL::square(p1.y() - p.y());
    RT D2 = CGAL::square(p2.x() - p.x()) + CGAL::square(p2.y() - p.y());

    RT d1 = CGAL::sqrt(D1) - p1.weight();
    RT d2 = CGAL::sqrt(D2) - p2.weight();

    return CGAL::compare(d1, d2);
  }

public:
  typedef Uncertain<Oriented_side>        result_type;
  struct argument_type {};

  inline
  Uncertain<Oriented_side>
  operator()(const Site_2& p1, const Site_2& p2,
	     const Point_2 &p) const
  {
    return - compare_distances(p1, p2, p, Method_tag());
  }

};


//--------------------------------------------------------------------

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_ORIENTED_SIDE_OF_BISECTOR_C2_H
