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



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//-----------------------------------------------------------------------------



template<class K, class Method_tag>
class Oriented_side_C2
  : public Basic_predicates_C2<K>
{
private:

  typedef Basic_predicates_C2<K>              Base;
  typedef Voronoi_vertex_C2<K,Method_tag>     Voronoi_vertex_2;

  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Line_2               Line_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;


  using Base::compute_supporting_line;
  using Base::compute_perpendicular;

public:
  typedef typename Base::Oriented_side        Oriented_side;
  typedef Oriented_side                       result_type;
  typedef Site_2                              argument_type;

  // computes the oriented side of the Voronoi vertex of s1, s2, inf
  // wrt the line that passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition(s1.is_point() || s2.is_point());
    CGAL_precondition(s1.is_segment() || s2.is_segment());

    return this->operator()((s1.is_point()? s1: s2), s, p);
  }

  // computes the oriented side of the point q
  // wrt the line that is passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& q,
                           const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( q.is_point() );
    CGAL_precondition( s.is_segment() && p.is_point() );

    Line_2 l = compute_supporting_line( s );
    Line_2 lp = compute_perpendicular(l, p.point());

    return lp.oriented_side(q.point());
  }

  // computes the oriented side of the Voronoi vertex of s1, s2, s3
  // wrt the line that is passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
                           const Site_2& s3,
                           const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    Voronoi_vertex_2 v(s1, s2, s3);
    Line_2 l = compute_supporting_line( s );
    Line_2 lp = compute_perpendicular(l, p.point());

    return v.oriented_side(lp);
  }
};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ORIENTED_SIDE_C2_H
