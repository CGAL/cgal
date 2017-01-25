// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_VERTEX_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_VERTEX_CONFLICT_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>

#ifdef CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
#ifndef CGAL_SDG_USE_OLD_INCIRCLE
#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_sqrt_field_C2.h>
#endif // CGAL_SDG_USE_OLD_INCIRCLE
#endif // CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//---------------------------------------------------------------------

template<class K, class Method_tag>
class Vertex_conflict_C2
{
private:
  typedef typename K::Point_2                Point_2;
  typedef typename K::Segment_2              Segment_2;
  typedef typename K::Site_2                 Site_2;
  typedef typename K::FT                     FT;
  typedef typename K::RT                     RT;
  typedef typename K::Orientation            Orientation;
  typedef typename K::Sign                   Sign;

  typedef Voronoi_vertex_C2<K,Method_tag>    Voronoi_vertex_2;

  typedef Are_same_points_C2<K>              Are_same_points_2;
  typedef Are_same_segments_C2<K>            Are_same_segments_2;

  typedef typename K::Intersections_tag      ITag;

private:
  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  bool is_on_common_support(const Site_2& s1, const Site_2& s2,
			    const Point_2& p) const
  {
    CGAL_precondition( !s1.is_input() && !s2.is_input() );

    if (  same_segments(s1.supporting_site(0),
			s2.supporting_site(0)) ||
	  same_segments(s1.supporting_site(0),
			s2.supporting_site(1))  ) {
      Site_2 support = s1.supporting_site(0);
      Site_2 tp = Site_2::construct_site_2(p);

      return (  same_points(support.source_site(), tp) ||
		same_points(support.target_site(), tp)  );
    } else if (  same_segments(s1.supporting_site(1),
			       s2.supporting_site(1)) ||
		 same_segments(s1.supporting_site(1),
			       s2.supporting_site(0))  ) {
      Site_2 support = s1.supporting_site(1);
      Site_2 tp = Site_2::construct_site_2(p);

      return (  same_points(support.source_site(), tp) ||
		same_points(support.target_site(), tp)  );      
    }
    return false;
  }

  bool have_common_support(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( !p.is_input() && !q.is_input() );

    return
      same_segments(p.supporting_site(0), q.supporting_site(0)) ||
      same_segments(p.supporting_site(0), q.supporting_site(1)) ||
      same_segments(p.supporting_site(1), q.supporting_site(1)) ||
      same_segments(p.supporting_site(1), q.supporting_site(0));
  }

  bool have_common_support(const Site_2& s, const Point_2& p1,
			   const Point_2& p2) const
  {
    CGAL_precondition( !s.is_input() );

    Site_2 t = Site_2::construct_site_2(p1, p2);

    return ( same_segments(s.supporting_site(0), t) ||
	     same_segments(s.supporting_site(1), t) );
  }

private:
  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_false&) const
  {
    Point_2 pp = p.point(), qp = q.point(), tp = t.point();

    // MK::ERROR: here I should call a kernel object, not a
    // function...; actually here (and everywhere in this class)
    // use the orientation predicate for sites; it does some
    // geometric filtering...
    Orientation o = orientation(pp, qp, tp);

    if ( o != COLLINEAR ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // MK::ERROR: change the following code to use the compare_x_2
    // and compare_y_2 stuff...
    RT dtpx = pp.x() - tp.x();
    RT dtpy = pp.y() - tp.y();
    RT dtqx = qp.x() - tp.x();
    RT minus_dtqy = -qp.y() + tp.y();
    
    Sign s = sign_of_determinant(dtpx, dtpy, minus_dtqy, dtqx);

    CGAL_assertion( s != ZERO );

    return s;
  }

  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_true&) const
  {
    Orientation o = COLLINEAR; // the initialization was done in
                               // order a compiler warning

    // do some geometric filtering...
    bool p_exact = p.is_input();
    bool q_exact = q.is_input();
    bool t_exact = t.is_input();
    bool filtered = false;
    // the following if-statement does the gometric filtering...
    // maybe it is not so important since this will only be
    // activated if a lot of intersection points appear on the
    // convex hull
    if ( !p_exact || !q_exact || !t_exact ) {
      if ( !p_exact && !q_exact && !t_exact ) {
	if ( have_common_support(p, q) &&
	     have_common_support(q, t) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && !q_exact && t_exact ) {
	if ( is_on_common_support(p, q, t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && q_exact && !t_exact ) {
	if ( is_on_common_support(p, t, q.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && !q_exact && !t_exact ) {
	if ( is_on_common_support(t, q, p.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( !p_exact && q_exact && t_exact ) {
	if ( have_common_support(p, q.point(), t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && !q_exact && t_exact ) {
	if ( have_common_support(q, p.point(), t.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      } else if ( p_exact && q_exact && !t_exact ) {
	if ( have_common_support(t, p.point(), q.point()) ) {
	  o = COLLINEAR;
	  filtered = true;
	}
      }
    }

    Point_2 pp = p.point(), qp = q.point(), tp = t.point();

    if ( !filtered ) {
      // MK::ERROR: here I should call a kernel object, not a
      // function...; actually here (and everywhere in this class)
      // use the orientation predicate for sites; it does some
      // geometric filtering...
      o = orientation(pp, qp, tp);
    }

    if ( o != COLLINEAR ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // MK::ERROR: change the following code to use the compare_x_2
    // and compare_y_2 stuff...
    RT dtpx = pp.x() - tp.x();
    RT dtpy = pp.y() - tp.y();
    RT dtqx = qp.x() - tp.x();
    RT minus_dtqy = -qp.y() + tp.y();
    
    Sign s = sign_of_determinant(dtpx, dtpy, minus_dtqy, dtqx);
    
    CGAL_assertion( s != ZERO );

    return s;
  }


  Sign incircle_p(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_point() );

    if ( p.is_point() && q.is_point() ) {

#if 1
      return incircle_ppp(p, q, t, ITag());

#else
      Orientation o = COLLINEAR; // the initialization was done in
                                 // order a compiler warning

      // do some geometric filtering...
      bool p_exact = p.is_input();
      bool q_exact = q.is_input();
      bool t_exact = t.is_input();
      bool filtered = false;
      // the following if-statement does the gometric filtering...
      // maybe it is not so important since this will only be
      // activated if a lot of intersection points appear on the
      // convex hull
      if ( !p_exact || !q_exact || !t_exact ) {
	if ( !p_exact && !q_exact && !t_exact ) {
	  if ( have_common_support(p, q) &&
	       have_common_support(q, t) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && !q_exact && t_exact ) {
	  if ( is_on_common_support(p, q, t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && q_exact && !t_exact ) {
	  if ( is_on_common_support(p, t, q.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && !q_exact && !t_exact ) {
	  if ( is_on_common_support(t, q, p.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( !p_exact && q_exact && t_exact ) {
	  if ( have_common_support(p, q.point(), t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && !q_exact && t_exact ) {
	  if ( have_common_support(q, p.point(), t.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	} else if ( p_exact && q_exact && !t_exact ) {
	  if ( have_common_support(t, p.point(), q.point()) ) {
	    o = COLLINEAR;
	    filtered = true;
	  }
	}
      }

      Point_2 pp = p.point(), qp = q.point(), tp = t.point();

      if ( !filtered ) {
	// MK::ERROR: here I should call a kernel object, not a
	// function...; actually here (and everywhere in this class)
	// use the orientation predicate for sites; it does some
	// geometric filtering...
	o = orientation(pp, qp, tp);
      }

      if ( o != COLLINEAR ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      // MK::ERROR: change the following code to use the compare_x_2
      // and compare_y_2 stuff...
      RT dtpx = pp.x() - tp.x();
      RT dtpy = pp.y() - tp.y();
      RT dtqx = qp.x() - tp.x();
      RT minus_dtqy = -qp.y() + tp.y();
      
      Sign s = sign_of_determinant(dtpx, dtpy, minus_dtqy, dtqx);

      CGAL_assertion( s != ZERO );

      return s;
#endif
    }

    CGAL_assertion( p.is_point() || q.is_point() );

    Orientation o;
    if ( p.is_point() && q.is_segment() ) {
      Point_2 pq = same_points(p, q.source_site()) ? q.target() : q.source();
      o = orientation(p.point(), pq, t.point());
    } else { // p is a segment and q is a point
      Point_2 pp = same_points(q, p.source_site()) ? p.target() : p.source();
      o = orientation(pp, q.point(), t.point());
    }
    if ( CGAL::is_certain(o == RIGHT_TURN) )
        return CGAL::get_certain( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
    return CGAL::Uncertain<CGAL::Sign>::indeterminate();
  }

  //-----------------------------------------------------------------------


  Sign incircle_pps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    bool is_p_tsrc = same_points(p, t.source_site());
    bool is_p_ttrg = same_points(p, t.target_site());

    bool is_q_tsrc = same_points(q, t.source_site());
    bool is_q_ttrg = same_points(q, t.target_site());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;
    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    if ( is_p_on_t && is_q_on_t ) {
	// if t is the segment joining p and q then t must be a vertex
	// on the convex hull
	return NEGATIVE;
    } else if ( is_p_on_t ) {
      // p is an endpoint of t
      // in this case the p,q,oo vertex is destroyed only if the
      // other endpoint of t is beyond
      Point_2 pt = is_p_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else if ( is_q_on_t ) {
      Point_2 pt = is_q_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // maybe here I should immediately return POSITIVE;
      // since we insert endpoints of segments first, p and q cannot
      // be consecutive points on the convex hull if one of the
      // endpoints of t is to the right of the line pq.
      Point_2 pp = p.point(), qq = q.point();
      Orientation o1 = orientation(pp, qq, t.source());
      Orientation o2 = orientation(pp, qq, t.target());

      if ( o1 == RIGHT_TURN || o2 == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return POSITIVE;
    }
  }


  Sign incircle_sps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_segment() && q.is_point() );

    bool is_q_tsrc = same_points(q, t.source_site());
    bool is_q_ttrg = same_points(q, t.target_site());

    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    if ( is_q_on_t ) {
      Point_2 pp = same_points(q, p.source_site()) ? p.target() : p.source();
      Point_2 pt = is_q_tsrc ? t.target() : t.source();

      Orientation o = orientation(pp, q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      return POSITIVE;
    }
  }


  Sign incircle_pss(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_segment() );

    bool is_p_tsrc = same_points(p, t.source_site());
    bool is_p_ttrg = same_points(p, t.target_site());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;

    if ( is_p_on_t ) {
      Point_2 pq = same_points(p, q.source_site()) ? q.target() : q.source();
      Point_2 pt = is_p_tsrc ? t.target() : t.source();

      Orientation o = orientation(p.point(), pq, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // if p is not an endpoint of t, then either p and q should
      // not be on the convex hull or t does not affect the vertex
      // of p and q.
      return POSITIVE;
    }
  }


  Sign incircle_s(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_segment() );

    if ( p.is_point() && q.is_point() ) {
      return incircle_pps(p, q, t);
    } else if ( p.is_point() && q.is_segment() ) {
      return incircle_pss(p, q, t);
    } else { // p is a segment and q is a point
      return incircle_sps(p, q, t);
    }
  }


public:
  typedef Site_2      argument_type;
  typedef Sign        result_type;


  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& t) const
  {
#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      int np = 0;
      if ( p.is_point() ) ++np;
      if ( q.is_point() ) ++np;
      if ( r.is_point() ) ++np;
      std::string suffix("-failure-log.cin");
      std::string fname;
      if ( np == 3 ) {
	fname = "ppp";
      } else if ( np == 2 ) {
	fname = "pps";
      } else if ( np == 1 ) {
	fname = "pss";
      } else {
	fname = "sss";
      }
      fname += suffix;
      std::ofstream ofs(fname.c_str(), std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl;
      ofs << q << std::endl;
      ofs << r << std::endl;
      ofs << t << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

#ifdef CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
#ifdef CGAL_SDG_USE_OLD_INCIRCLE
    typedef Voronoi_vertex_sqrt_field_new_C2<K>   Alt_Voronoi_vertex_2;
#else
    typedef Voronoi_vertex_sqrt_field_C2<K>       Alt_Voronoi_vertex_2;
#endif

    Voronoi_vertex_2 v(p, q, r);
    Alt_Voronoi_vertex_2 v_alt(p, q, r);

    Sign s = v.incircle(t);
    Sign s_alt = v_alt.incircle(t);

    if ( s != s_alt ) {
      std::cerr << "different results" << std::endl;
      std::cerr << p << std::endl;
      std::cerr << q << std::endl;
      std::cerr << r << std::endl;
      std::cerr << t << std::endl;
      CGAL_assertion( s == s_alt );
      exit(1);
    }

    return s;
#else
    Voronoi_vertex_2 v(p, q, r);

    return v.incircle(t);
#endif // CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
  }


  

  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl;
      ofs << q << std::endl;
      ofs << t << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( same_points(p, q.source_site()) ||
		      same_points(p, q.target_site()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( same_points(p.source_site(), q) ||
		      same_points(p.target_site(), q) );
    }

    if ( t.is_point() ) {
      //      return incircle_p(p, q, t);
      return incircle_p(q, p, t);
    }

    // MK::ERROR: do geometric filtering when orientation is called.
    //    return incircle_s(p, q, t);
    return incircle_s(q, p, t);
  }


};

//---------------------------------------------------------------------

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_VERTEX_CONFLICT_C2_H
