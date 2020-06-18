// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//-----------------------------------------------------------------------------

template<class K, class Method_tag>
class Infinite_edge_interior_conflict_C2
{
public:
  typedef typename K::Site_2           Site_2;
  typedef typename K::RT               RT;
  typedef typename K::Boolean          Boolean;
  typedef Are_same_points_C2<K>        Are_same_points_2;
  typedef Are_same_segments_C2<K>      Are_same_segments_2;

  typedef Boolean                      result_type;
  struct argument_type {};

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

public:
  Boolean   operator()(const Site_2& q, const Site_2& s, const Site_2& r,
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

      Sign sgn1 = sign_of_determinant(dtsx, dtsy, minus_dtqy, dtqx);

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

      Sign sgn1 = sign_of_determinant(dtssx, dtssy, minus_dtqy, dtqx);

      CGAL_assertion( sgn1 != ZERO );

      return (sgn1 == POSITIVE);
    }

    return ( sgn == NEGATIVE );
  }

};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
