// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_VERTEX_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_VERTEX_CONFLICT_C2_H

#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/nox/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/nox/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/nox/Is_endpoint_of_segment_C2.h>


#ifdef CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
#ifndef CGAL_SDG_USE_OLD_INCIRCLE
#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_sqrt_field_C2.h>
#endif // CGAL_SDG_USE_OLD_INCIRCLE
#endif // CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY

CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

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
  typedef Is_endpoint_of_segment_C2<K>       Is_endpoint_of_segment_2;

  typedef typename K::Intersections_tag      ITag;

private:
  Are_same_points_2         same_points;
  Are_same_segments_2       same_segments;
  Is_endpoint_of_segment_2  is_endpoint_of;

private:
  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_false&) const
  {
    const Point_2& pp = p.point();
    const Point_2& qp = q.point();
    const Point_2& tp = t.point();

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

    const Point_2& pp = p.point();
    const Point_2& qp = q.point();
    const Point_2& tp = t.point();

    // MK::ERROR: here I should call a kernel object, not a
    // function...; actually here (and everywhere in this class)
    // use the orientation predicate for sites; it does some
    // geometric filtering...
    o = orientation(pp, qp, tp);

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
      return incircle_ppp(p, q, t, ITag());
    }

    CGAL_assertion( p.is_point() || q.is_point() );

    Orientation o;
    if ( p.is_point() && q.is_segment() ) {
      const Point_2& pq =
	same_points(p.point(), q.source()) ? q.target() : q.source();
      o = orientation(p.point(), pq, t.point());
    } else { // p is a segment and q is a point
      const Point_2& pp =
	same_points(q.point(), p.source()) ? p.target() : p.source();
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

    bool is_p_tsrc = same_points(p.point(), t.source());
    bool is_p_ttrg = same_points(p.point(), t.target());

    bool is_q_tsrc = same_points(q.point(), t.source());
    bool is_q_ttrg = same_points(q.point(), t.target());

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
      const Point_2& pt = is_p_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else if ( is_q_on_t ) {
      const Point_2& pt = is_q_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // maybe here I should immediately return POSITIVE;
      // since we insert endpoints of segments first, p and q cannot
      // be consecutive points on the convex hull if one of the
      // endpoints of t is to the right of the line pq.
      const Point_2& pp = p.point();
      const Point_2& qq = q.point();
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

    bool is_q_tsrc = same_points(q.point(), t.source());
    bool is_q_ttrg = same_points(q.point(), t.target());

    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    if ( is_q_on_t ) {
      const Point_2& pp =
	same_points(q.point(), p.source()) ? p.target() : p.source();
      const Point_2& pt = is_q_tsrc ? t.target() : t.source();

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

    bool is_p_tsrc = same_points(p.point(), t.source());
    bool is_p_ttrg = same_points(p.point(), t.target());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;

    if ( is_p_on_t ) {
      Point_2 pq = same_points(p.point(), q.source()) ? q.target() : q.source();
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

#if !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG)
    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( is_endpoint_of(p, q) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( is_endpoint_of(q, p) );
    }
#endif

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

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_VERTEX_CONFLICT_C2_H
