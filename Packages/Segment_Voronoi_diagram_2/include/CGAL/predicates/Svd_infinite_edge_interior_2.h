// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H
#define CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H

#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>
#include <CGAL/predicates/Svd_are_same_segments_C2.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------

template<class K, class Method_tag>
class Svd_infinite_edge_interior_2
{
public:
  typedef typename K::Site_2           Site_2;
  typedef typename K::RT               RT;
  typedef Svd_are_same_points_C2<K>    Are_same_points_2;
  typedef Svd_are_same_segments_C2<K>  Are_same_segments_2;

  typedef bool                         result_type;
  struct argument_type {};
  typedef Arity_tag<5>                 Arity;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

public:
  bool operator()(const Site_2& q, const Site_2& s, const Site_2& r,
		  const Site_2& t, Sign sgn) const
  {
    if ( t.is_segment() ) {
      return false;
    }

    if ( q.is_segment() ) {
      // in this case r and s must be endpoints of q
      return ( sgn == NEGATIVE );
    }

    if ( s.is_point() && r.is_point() && same_points(s, r) ) {
      // MK::ERROR: write this code using the compare_x_2 and
      //    compare_y_2 predicates instead of computing the inner
      //    product...
      RT dtsx = s.point().x() - t.point().x();
      RT dtsy = s.point().y() - t.point().y();
      RT dtqx = q.point().x() - t.point().x();
      RT minus_dtqy = -q.point().y() + t.point().y();

      Sign sgn1 = sign_of_determinant2x2(dtsx, dtsy, minus_dtqy, dtqx);

      CGAL_assertion( sgn1 != ZERO );

      return (sgn1 == POSITIVE);
    }

    if ( s.is_segment() && r.is_segment() && same_segments(s, r) ) {
      CGAL_assertion( same_points(q, s.source_site()) ||
		      same_points(q, s.target_site()) );
      Site_2 ss;
      if ( same_points(q, s.source_site()) ) {
	ss = s.target_site();
      } else {
	ss = s.source_site();
      }
      // MK::ERROR: write this code using the compare_x_2 and
      //    compare_y_2 predicates instead of computing the inner
      //    product...
      RT dtssx = ss.point().x() - t.point().x();
      RT dtssy = ss.point().y() - t.point().y();
      RT dtqx = q.point().x() - t.point().x();
      RT minus_dtqy = -q.point().y() + t.point().y();

      Sign sgn1 = sign_of_determinant2x2(dtssx, dtssy, minus_dtqy, dtqx);

      CGAL_assertion( sgn1 != ZERO );

      return (sgn1 == POSITIVE);
    }

    return ( sgn == NEGATIVE );
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_INFINITE_EDGE_INTERIOR_2_H
