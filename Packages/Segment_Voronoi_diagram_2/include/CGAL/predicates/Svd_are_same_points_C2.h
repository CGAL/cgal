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





#ifndef CGAL_SVD_ARE_SAME_POINTS_C2_H
#define CGAL_SVD_ARE_SAME_POINTS_C2_H

CGAL_BEGIN_NAMESPACE

template<class K>
class Svd_are_same_points_C2
{
private:
  typedef typename K::Point_2     Point_2;
  typedef typename K::Segment_2   Segment_2;
  typedef typename K::Site_2      Site_2;
  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  typedef typename K::Intersections_tag  ITag;

  Compare_x_2 compare_x_2;
  Compare_y_2 compare_y_2;

  bool are_same(const Point_2& p, const Point_2& q) const
  {
    return
      compare_x_2(p, q) == EQUAL && compare_y_2(p, q) == EQUAL;
  }

  bool are_same(const Segment_2& s, const Segment_2& t) const
  {
    return
      ( are_same(s.source(), t.source()) &&
	are_same(s.target(), t.target()) ) ||
      ( are_same(s.source(), t.target()) &&
	are_same(s.target(), t.source()) );
  }

  bool predicate(const Site_2& p, const Site_2& q, const Tag_false&) const
  {
    return are_same(p.point(), q.point()); 
  }

  bool predicate(const Site_2& p, const Site_2& q, const Tag_true&) const
  {
    if ( !p.is_exact() && !q.is_exact() ) {
      Segment_2 s[2] = { p.supporting_segment(0),
			 p.supporting_segment(1) };
      Segment_2 t[2] = { q.supporting_segment(0),
			 q.supporting_segment(1) };

      if (  ( are_same(s[0], t[0]) && are_same(s[1], t[1]) ) ||
	    ( are_same(s[0], t[1]) && are_same(s[1], t[0]) )  ) {
	return true;
      }
    }

    return predicate(p, q, Tag_false());
  }

public:
  typedef bool           result_type;
  typedef Site_2         argument_type;
  typedef Arity_tag<2>   Arity;

  bool operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    return predicate(p, q, ITag());
  }
};



CGAL_END_NAMESPACE

#endif // CGAL_SVD_ARE_SAME_POINTS_C2_H
