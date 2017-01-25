// Copyright (c) 2015  Università della Svizzera italiana.
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
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_LINF_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_LINF_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

#include <CGAL/Orientation_Linf_2.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

//-----------------------------------------------------------------------------



template<class K>
class Orientation_Linf_C2
  : private Basic_predicates_C2<K>
{
private:
  typedef Basic_predicates_C2<K>              Base;

private:
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;

  typedef typename Base::Comparison_result    Comparison_result;
  typedef typename Base::Oriented_side        Oriented_side;
  typedef typename Base::Sign                 Sign;

  // the orientation Linf predicate for three points 
  typedef Orientation_Linf_2<K>               Orientation_Linf_points_2;
  
  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;

  typedef typename K::Intersections_tag       ITag;

  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  bool have_common_support(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( !p.is_input() && !q.is_input() );

    return
      same_segments(p.supporting_site(0), q.supporting_site(0)) ||
      same_segments(p.supporting_site(0), q.supporting_site(1)) ||
      same_segments(p.supporting_site(1), q.supporting_site(0)) ||
      same_segments(p.supporting_site(1), q.supporting_site(1));
  }

  bool have_common_support(const Site_2& p, const Site_2& q,
			   Site_2& support) const
  {
    CGAL_precondition( !p.is_input() && !q.is_input() );

    if ( same_segments(p.supporting_site(0),
		       q.supporting_site(0)) ||
	 same_segments(p.supporting_site(0),
		       q.supporting_site(1)) ) {
      support = p.supporting_site(0);
      return true;
    } else if ( same_segments(p.supporting_site(1),
			      q.supporting_site(0)) ||
		same_segments(p.supporting_site(1),
			      q.supporting_site(1)) ) {
      support =  p.supporting_site(1);
      return true;
    }
    return false;
  }

  bool is_endpoint(const Site_2& p, const Site_2& s) const
  {
    return same_points(s.source_site(), p) ||
      same_points(s.target_site(), p);
  }

  //-------------------------------------------------------------

  OrientationLinf predicate(const Site_2& p, const Site_2& q,
			const Site_2& r, const Tag_false&) const
  {
    return Orientation_Linf_points_2()(
	    p.point(), q.point(), r.point());
  }

  OrientationLinf predicate(const Site_2& p, const Site_2& q,
			const Site_2& r, const Tag_true&) const
  {
#if 1
    // do geometric filtering
    bool pe = p.is_input();
    bool qe = q.is_input();
    bool re = r.is_input();
    Site_2 support;
    if ( !pe && !qe && !re ) {
      if ( have_common_support(p, q, support) &&
	   have_common_support(support, r) ) {
	return DEGENERATE;
      }
    } else if ( !pe && !qe ) {
      if ( have_common_support(p, q, support) &&
	   is_endpoint(r, support) ) {
	return DEGENERATE;
      }
    } else if ( !pe && !re ) {
      if ( have_common_support(p, r, support) &&
	   is_endpoint(q, support) ) {
	return DEGENERATE;
      }
    } else if ( !qe && !re ) {
      if ( have_common_support(q, r, support) &&
	   is_endpoint(p, support) ) {
	return DEGENERATE;
      }
    } else if ( !pe ) {
      Site_2 s0 = p.supporting_site(0);
      Site_2 s1 = p.supporting_site(1);
      if ( (is_endpoint(q, s0) && is_endpoint(r, s0)) ||
	   (is_endpoint(q, s1) && is_endpoint(r, s1)) ) {
	return DEGENERATE;
      }
    } else if ( !qe ) {
      Site_2 s0 = q.supporting_site(0);
      Site_2 s1 = q.supporting_site(1);
      if ( (is_endpoint(p, s0) && is_endpoint(r, s0)) ||
	   (is_endpoint(p, s1) && is_endpoint(r, s1)) ) {
	return DEGENERATE;
      }
    } else if ( !re ) {
      Site_2 s0 = r.supporting_site(0);
      Site_2 s1 = r.supporting_site(1);
      if ( (is_endpoint(q, s0) && is_endpoint(p, s0)) ||
	   (is_endpoint(q, s1) && is_endpoint(p, s1)) ) {
	return DEGENERATE;
      }
    }
#endif

    return predicate(p, q, r, Tag_false());
  }

public:
  typedef OrientationLinf              result_type;
  typedef Site_2                       argument_type;

  OrientationLinf operator()(const Site_2& p, const Site_2& q,
                             const Site_2& r) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && r.is_point() );

    return predicate(p, q, r, ITag());
  }
};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTATION_LINF_C2_H
