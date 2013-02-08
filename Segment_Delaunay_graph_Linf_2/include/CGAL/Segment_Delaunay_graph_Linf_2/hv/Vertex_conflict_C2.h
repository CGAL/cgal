
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_VERTEX_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_VERTEX_CONFLICT_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/Orientation_Linf_2.h>

#ifdef CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
#ifndef CGAL_SDG_USE_OLD_INCIRCLE
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_sqrt_field_C2.h>
#endif // CGAL_SDG_USE_OLD_INCIRCLE
#endif // CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2_hv {

//---------------------------------------------------------------------

template<class K, class Method_tag>
class Vertex_conflict_C2
  : public SegmentDelaunayGraphLinf_2::Basic_predicates_C2<K>
{
private:

  typedef SegmentDelaunayGraphLinf_2::Basic_predicates_C2<K> Base;
  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  using Base::compute_linf_perpendicular;
  using Base::compute_line_from_to;
  using Base::compute_horizontal_side_line;
  using Base::compute_vertical_side_line;
  using Base::compute_perpendicular;
  using Base::compute_cw_perpendicular;
  using Base::compute_horizontal_projection;
  using Base::compute_vertical_projection;
  using Base::intersects_segment_interior_inf_box;
  using Base::intersects_segment_positive_of_wedge;
  using Base::intersects_segment_negative_of_wedge;
  using Base::is_site_h_or_v;

  typedef typename Base::Line_2              Line_2;

  typedef typename K::Point_2                Point_2;
  typedef typename K::Segment_2              Segment_2;
  typedef typename K::Site_2                 Site_2;
  typedef typename K::FT                     FT;
  typedef typename K::RT                     RT;
  //typedef typename K::Orientation            Orientation;
  typedef CGAL::OrientationLinf              OrientationLinf;
  typedef typename K::Sign                   Sign;

  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::
          Voronoi_vertex_C2<K,Method_tag>
          Voronoi_vertex_2;

  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Are_same_points_C2<K>
          Are_same_points_2;
  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Are_same_segments_C2<K>
          Are_same_segments_2;

  typedef typename K::Intersections_tag      ITag;

  // the orientation Linf predicate for three points
  typedef Orientation_Linf_2<K>              Orientation_Linf_points_2;

  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  Compare_x_2 cmpx;
  Compare_y_2 cmpy;
  Orientation_Linf_points_2 orientation_linf;


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
    // philaris: I think we call an object now
    OrientationLinf o = orientation_linf(pp, qp, tp);

    if ( o != DEGENERATE ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // here orientation Linf is degenerate

    // return NEGATIVE if t is strictly inside the bounding box
    // of p and q, otherwise return POSITIVE;
    // ZERO is not possible (because t should not be equal
    // to either p or q)

    Comparison_result cmpxpt = cmpx(pp, tp);
    Comparison_result cmpxtq = cmpx(tp, qp);
    Comparison_result cmpypt = cmpy(pp, tp);
    Comparison_result cmpytq = cmpy(tp, qp);

    //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
    Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

    CGAL_assertion( s != ZERO );

    return s;
  }

  Sign incircle_ppp(const Site_2& p, const Site_2& q,
		    const Site_2& t, const Tag_true&) const
  {
    Orientation o = DEGENERATE; // the initialization was done in
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
      // philaris: I think we call an object now
      o = orientation_linf(pp, qp, tp);
    }

    if ( o != DEGENERATE ) {
      return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
    }

    // here orientation Linf is degenerate

    // return NEGATIVE if t is strictly inside the bounding box
    // of p and q, otherwise return POSITIVE;
    // ZERO is not possible (because t should not be equal
    // to either p or q)

    Comparison_result cmpxpt = cmpx(pp, tp);
    Comparison_result cmpxtq = cmpx(tp, qp);
    Comparison_result cmpypt = cmpy(pp, tp);
    Comparison_result cmpytq = cmpy(tp, qp);

    //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
    Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

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
      Orientation o = DEGENERATE; // the initialization was done in
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
	// philaris: I think we call an object now
	o = orientation_linf(pp, qp, tp);
      }

      if ( o != DEGENERATE ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      // here orientation Linf is degenerate

      // return NEGATIVE if t is strictly inside the bounding box
      // of p and q, otherwise return POSITIVE;
      // ZERO is not possible (because t should not be equal
      // to either p or q)

      Comparison_result cmpxpt = cmpx(pp, tp);
      Comparison_result cmpxtq = cmpx(tp, qp);
      Comparison_result cmpypt = cmpy(pp, tp);
      Comparison_result cmpytq = cmpy(tp, qp);

      //Sign s = -sign_of( cmpxpt * cmpxtq + cmpypt * cmpytq );
      Sign s = CGAL::compare(0, cmpxpt * cmpxtq + cmpypt * cmpytq);

      CGAL_assertion( s != ZERO );

      return s;
#endif
    }

    // here the pair of p, q contains a point and a segment

    CGAL_assertion( p.is_point() || q.is_point() );

    CGAL_SDG_DEBUG( std::cout << "debug incircle_p of p=" << p
                   << " q=" << q << " t=" << t << std::endl; );

    // philaris: here, for Linf there are differences from L2

    if ((q.is_segment() and is_site_h_or_v(q)) or
        (p.is_segment() and is_site_h_or_v(p))   ) {
      // the case of horizontal or vertical segment is the same
      // as in L2
      Orientation o;
      if ( p.is_point() && q.is_segment() ) {
        CGAL_assertion( same_points(p, q.source_site()) or
            same_points(p, q.target_site()) );
        Point_2 pq = same_points(p, q.source_site()) ? q.target() : q.source();
        o = orientation_linf(p.point(), pq, t.point());
      } else { // p is a segment and q is a point
        CGAL_assertion( same_points(q, p.source_site()) or
            same_points(q, p.target_site()) );
        Point_2 pp = same_points(q, p.source_site()) ? p.target() : p.source();
        o = orientation_linf(pp, q.point(), t.point());
      }
      if ( CGAL::is_certain(o == RIGHT_TURN) )  {
        CGAL_SDG_DEBUG(std::cout << "debug incircle_p about to return "
          << (CGAL::get_certain( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE)
                       << std::endl ;);
        return CGAL::get_certain( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
      }
      return CGAL::Uncertain<CGAL::Sign>::indeterminate();
    }
    //Sandeep: Added to avoid compiler warning
    CGAL_assertion(false);
    CGAL_SDG_DEBUG(std::cout
        << "sandeep: debug incircle_p should not reach here "
        << NEGATIVE
        << std::endl ;);
    return NEGATIVE;
  }

  //-----------------------------------------------------------------------

  // philaris: here Linf is very different from L2
  Sign incircle_pps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_pps p=" << p
      << " q=" << q << " t=" << t << std::endl; );

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
    } else {
      // here, neither p nor q is on t
      if (is_site_h_or_v(t)) {
        // here segment is horizontal or vertical;
        // there are similarities with L2 case here

        if ( is_p_on_t ) {
          // p is an endpoint of t
          // in this case the p,q,oo vertex is destroyed only if the
          // other endpoint of t is beyond
          Point_2 pt = is_p_tsrc ? t.target() : t.source();
          Orientation o = CGAL::orientation(p.point(), q.point(), pt);
          CGAL_SDG_DEBUG(std::cout << "debug incircle_pps p is endp of t: "
                         << "p = " << p << ", q = " << q << ", t = " << t
                         << " o = " << o << std::endl; );
          return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
        } else if ( is_q_on_t ) {
          Point_2 pt = is_q_tsrc ? t.target() : t.source();
          Orientation o = CGAL::orientation(p.point(), q.point(), pt);
          CGAL_SDG_DEBUG(std::cout << "debug incircle_pps q is endp of t: "
                         << "p = " << p << ", q = " << q << ", t = " << t
                         << " o = " << o << std::endl; );
          return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
        } else {
          // We immediately return POSITIVE here.
          // This only works because the segment t is hor. or ver. !
          // If t is neither hor. or ver., we have a more complicated
          // decision procedure (see below).
          // since we insert endpoints of segments first, p and q cannot
          // be consecutive points on the convex hull, because of the
          // endpoints of t
          return POSITIVE;
        }
      }
    } // end of case where neither p nor q is on t
    //Sandeep: Added to avoid compiler warning
    CGAL_assertion(false);
    CGAL_SDG_DEBUG(std::cout << "sandeep: debug incircle_pps should not reach here "
                   << NEGATIVE
                   << std::endl ;);
    return NEGATIVE;
  } // end of incircle_pps


  Sign incircle_sps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_segment() && q.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_sps p=" << p
      << " q=" << q << " t=" << t << std::endl; );

    bool is_q_tsrc = same_points(q, t.source_site());
    bool is_q_ttrg = same_points(q, t.target_site());

    bool is_q_on_t = is_q_tsrc or is_q_ttrg;

    bool is_q_psrc = same_points(q, p.source_site());
    bool is_q_ptrg = same_points(q, p.target_site());

    bool is_q_on_p = is_q_psrc or is_q_ptrg;


    if ( is_q_on_t and is_q_on_p) {
      CGAL_SDG_DEBUG(std::cout << "debug incircle_sps "
        << "is_q_on_t and is_q_on_p" << std::endl; );

      Point_2 pp = same_points(q, p.source_site()) ? p.target() : p.source();
      Point_2 pt = is_q_tsrc ? t.target() : t.source();

      Orientation o = orientation_linf(pp, q.point(), pt);

                     CGAL_SDG_DEBUG(std::cout << "debug incircle_sps or( "
        << pp << ", " << q.point() << ", " << pt << " ) = "
        << o << std::endl; );

      if (o == DEGENERATE) {

        Line_2 lpq = compute_line_from_to(pp, q.point());
        Line_2 lperp = compute_linf_perpendicular(lpq, q.point());

        Oriented_side os = oriented_side_of_line(lperp, pt);

        CGAL_assertion( os != ON_ORIENTED_BOUNDARY );

        if (os == ON_NEGATIVE_SIDE) {
          CGAL_SDG_DEBUG(std::cout << "debug incircle_sps degen os=NEG_SIDE "
            << "returns POS" << std::endl; );
          return POSITIVE;
        } else { // here os == ON_POSITIVE_SIDE
          // here, we must make a refinement, based on turn
          Orientation or_l2 = CGAL::orientation(pp, q.point(), pt);

          CGAL_assertion( or_l2 != COLLINEAR );

          CGAL_SDG_DEBUG(std::cout << "debug incircle_sps degen os=POS_SIDE "
            << "returns "
            << ((or_l2 == RIGHT_TURN)? NEGATIVE : POSITIVE)
            << std::endl; );

          return (or_l2 == RIGHT_TURN)? NEGATIVE : POSITIVE;
        }

      }

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // philaris: serious difference for Linf here, related to L2

      if ( is_site_h_or_v(p) ) {
        return POSITIVE;
      }
    } // end of case where q is not on t
    //Sandeep: Added to avoid compiler warning
    CGAL_assertion(false);
    CGAL_SDG_DEBUG(std::cout << "sandeep: debug incircle_sps should not reach here "
                   << NEGATIVE
                   << std::endl ;);
    return NEGATIVE;
  } // end of function incircle_sps


  Sign incircle_pss(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_pss (pqt) = (" <<
      p << ") (" << q << ") (" << t << ")" << std::endl; );

    bool is_p_tsrc = same_points(p, t.source_site());
    bool is_p_ttrg = same_points(p, t.target_site());

    bool is_p_on_t = is_p_tsrc or is_p_ttrg;

    bool is_p_qsrc = same_points(p, q.source_site());
    bool is_p_qtrg = same_points(p, q.target_site());

    bool is_p_on_q = is_p_qsrc or is_p_qtrg;

    if ( is_p_on_t and is_p_on_q ) {

      CGAL_SDG_DEBUG(std::cout << "debug incircle_pss "
        << "is_p_on_t and is_p_on_q" << std::endl; );

      Point_2 pq = same_points(p, q.source_site()) ? q.target() : q.source();
      Point_2 pt = is_p_tsrc ? t.target() : t.source();

      Orientation o = orientation_linf(p.point(), pq, pt);

      if (o == DEGENERATE) {

        Line_2 lpq = compute_line_from_to(pq, p.point());
        Line_2 lperp = compute_linf_perpendicular(lpq, p.point());

        Oriented_side os = oriented_side_of_line(lperp, pt);

        CGAL_assertion( os != ON_ORIENTED_BOUNDARY );

        if (os == ON_NEGATIVE_SIDE) {
          CGAL_SDG_DEBUG(std::cout << "debug incircle_pss degen os=NEG_SIDE "
            << "returns POS" << std::endl; );
          return POSITIVE;
        } else { // here os == ON_POSITIVE_SIDE
          // here, we must make a refinement, based on turn
          Orientation or_l2 = CGAL::orientation(p.point(), pq, pt);

          CGAL_assertion( or_l2 != COLLINEAR );

          CGAL_SDG_DEBUG(std::cout << "debug incircle_pss degen os=POS_SIDE "
            << "returns "
            << ((or_l2 == RIGHT_TURN)? NEGATIVE : POSITIVE)
            << std::endl; );

          return (or_l2 == RIGHT_TURN)? NEGATIVE : POSITIVE;
        }
      }

        CGAL_SDG_DEBUG(std::cout << "debug incircle_pss or( "
            << p.point() << ", " << pq << ", " << pt << " ) = "
            << o << std::endl; );

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // philaris: serious difference for Linf here, related to L2

      if ( is_site_h_or_v(q) ) {
        return POSITIVE;
      }
    } // end of case where p is not on t
    //Sandeep: Added to avoid compiler warning
    CGAL_assertion(false);
    CGAL_SDG_DEBUG(std::cout << "sandeep: debug incircle_pss should not reach here "
                   << NEGATIVE
                   << std::endl ;);
    return NEGATIVE;
  } // end of function incircle_pss


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

    CGAL_SDG_DEBUG(std::cout
        << "debug: Vertex_conflict entering (pqrt)= ("
        << p << ") (" << q << ") (" << r << ") (" << t << ")"
        << std::endl; );


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

    Point_2 temppt = v.point();
    CGAL_SDG_DEBUG(std::cout
        << "debug Voronoi_vertex_2 v=" << temppt << std::endl; );

    Sign retval = v.incircle(t);

              CGAL_SDG_DEBUG(std::cout
                  << "debug: Vertex_conflict (pqrt)= ("
                  << p << ") (" << q << ") ("
                  << r << ") (" << t << ") returns " << retval
                  << std::endl; );

    return retval;
#endif // CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY
  }


  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug: Vertex_conflict (pqt)= ("
              << p << ") (" << q << ") (" << t << ")"
              << std::endl; );

#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << p << std::endl; );
      ofs << q << std::endl; );
      ofs << t << std::endl; );
      ofs << "=======" << std::endl; );
      ofs.close();
    }
#endif

    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    // philaris: the assertion that the point is endpoint of
    // the segment is not true in Linf (in contrast with L2)
    /*
    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( same_points(p, q.source_site()) ||
		      same_points(p, q.target_site()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( same_points(p.source_site(), q) ||
		      same_points(p.target_site(), q) );
    }
    */

    Sign retval;

    if ( t.is_point() ) {
      //      return incircle_p(p, q, t);
      retval = incircle_p(q, p, t);
    } else {
      // MK::ERROR: do geometric filtering when orientation is called.
      //    return incircle_s(p, q, t);
      retval = incircle_s(q, p, t);
    }

              CGAL_SDG_DEBUG(std::cout << "debug: Vertex_conflict (pqt)= ("
              << p << ") (" << q << ") (" << t << ")"
              << " returns " << retval << std::endl; );
    return retval;
  }


};

//---------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2_hv

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_VERTEX_CONFLICT_C2_H
