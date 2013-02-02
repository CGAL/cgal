
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
//#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Orientation_Linf_2.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2_hv {

//-----------------------------------------------------------------------------

template<class K, class Method_tag>
class Infinite_edge_interior_conflict_C2
  : public CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Basic_predicates_C2<K>
{
public:

  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_NS::Basic_predicates_C2<K>
          Base;

  typedef typename K::Site_2           Site_2;
  typedef typename K::Point_2          Point_2;
  typedef typename K::RT               RT;
  typedef typename K::Boolean          Boolean;

  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  Compare_x_2 cmpx;
  Compare_y_2 cmpy;

  typedef Boolean                      result_type;
  struct argument_type {};

private:
  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Are_same_points_C2<K>
          Are_same_points_2;
  typedef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Are_same_segments_C2<K>
          Are_same_segments_2;

  Are_same_points_2    same_points;
  Are_same_segments_2  same_segments;

  typedef Orientation_Linf_2<K>        Orientation_Linf_2_Type;

  Orientation_Linf_2_Type  or_linf;

  using Base::bounded_side_of_bbox;
  using Base::compute_line_from_to;
  using Base::oriented_side_of_line;
  using Base::compute_supporting_line;
  using Base::compute_vertical_projection;
  using Base::compute_horizontal_projection;
  using Base::compute_linf_projection_nonhom;

  typedef typename Base::Line_2        Line_2;

public:
  Boolean   operator()(const Site_2& q, const Site_2& s, const Site_2& r,
		       const Site_2& t, Sign sgn) const
  {

    CGAL_SDG_DEBUG(
        std::cout << "debug infinite-edge-int-cf entering (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << std::endl;);

    if ( t.is_segment() ) {
      if (q.is_point() and s.is_point() and r.is_point()) {
        if (sgn == NEGATIVE) {
          CGAL_SDG_DEBUG(std::cout << "debug return tocheck" << std::endl;);
          if (same_points(q, t.source_site()) or
              same_points(q, t.target_site())   ) {
            // this works because endpoints of a segment are
            // inserted before its interior
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << false << std::endl;);
            return false;
          } else {

            // here q is point
            CGAL_SDG_DEBUG(
                std::cout
                << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
                << q << ' ' << s << ' ' << r << ' ' << t
                << ' ' << sgn << " returns "
                << true << std::endl;);
            return true;
          }
        }
      } // end of case where q, s, r are all points

      CGAL_SDG_DEBUG(std::cout
          << "debug tocheck q,s,r not all points" << std::endl;);

      if (sgn == NEGATIVE) {
        CGAL_SDG_DEBUG(std::cout
            << "debug not all pts return true tocheck" << std::endl;);

        if (q.is_point()) {
          CGAL_assertion(
              not (same_points(q, t.source_site()) or
                   same_points(q, t.target_site())   ) );

          //CGAL_assertion(false);

          // philaris: tocheck more

          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << true << std::endl;);

          return true;

        }

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << true << std::endl;);

        return true;
      }

      // philaris: tocheck
      CGAL_SDG_DEBUG(
          std::cout
          << "debug infinite-edge-int-cf tocheck (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << " returns "
          << false << std::endl;);
      return false;
    } // end of case where t is segment

    // here and below, t is always a point

    if ( q.is_segment() ) {

      // philaris: difference from L2 here;
      // in L2, r and s are endpoints of q
      // in Linf they still have to be points, but
      // they do not have to be endpoints of q
      // (this has to be checked)
      CGAL_assertion(s.is_point() and r.is_point());

      // in this case r and s must be endpoints of q
      CGAL_SDG_DEBUG(
          std::cout
          << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << " returns "
          << ( sgn == NEGATIVE ) << std::endl;);
      return ( sgn == NEGATIVE );
    }

    // here and below both q and t are points

    if ( s.is_point() && r.is_point() && same_points(s, r) ) {
      // MK::ERROR: write this code using the compare_x_2 and
      //    compare_y_2 predicates instead of computing the inner
      //    product...
      // philaris: adaptation to Linf

      CGAL_SDG_DEBUG(std::cout <<
          "debug infcf s=r and is point" << std::endl;);

      Comparison_result cmpxst = cmpx(s.point(), t.point());
      Comparison_result cmpxtq = cmpx(t.point(), q.point());
      Comparison_result cmpyst = cmpy(s.point(), t.point());
      Comparison_result cmpytq = cmpy(t.point(), q.point());

      //Sign sgn1 = -sign_of( cmpxst * cmpxtq + cmpyst * cmpytq );
      Sign sgn1 = CGAL::compare(0, cmpxst * cmpxtq + cmpyst * cmpytq);

      CGAL_assertion( sgn1 != ZERO );

      CGAL_SDG_DEBUG(
          std::cout
          << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << " returns "
          << (sgn1 == POSITIVE) << std::endl;);
      return (sgn1 == POSITIVE);
    }

    // still here, both q and t are points

    if ( s.is_segment() && r.is_segment() && same_segments(s, r) ) {
      CGAL_SDG_DEBUG(std::cout <<
          "debug infcf s=r and is segment" << std::endl;);

      if ( same_points(q, s.source_site()) or
           same_points(q, s.target_site())   )  {
        // here point q is endpoint of segment s
        Site_2 ss;
        if ( same_points(q, s.source_site()) ) {
          ss = s.target_site();
        } else {
          ss = s.source_site();
        }
        // MK::ERROR: write this code using the compare_x_2 and
        //    compare_y_2 predicates instead of computing the inner
        //    product...
        // philaris: adaptation to Linf

        Comparison_result cmpxst = cmpx(ss.point(), t.point());
        Comparison_result cmpxtq = cmpx(t.point(), q.point());
        Comparison_result cmpyst = cmpy(ss.point(), t.point());
        Comparison_result cmpytq = cmpy(t.point(), q.point());

        //Sign sgn1 = -sign_of( cmpxst * cmpxtq + cmpyst * cmpytq );
        Sign sgn1 = CGAL::compare(0, cmpxst * cmpxtq + cmpyst * cmpytq);

        CGAL_assertion( sgn1 != ZERO );

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << (sgn1 == POSITIVE) << std::endl;);

        return (sgn1 == POSITIVE);
      }
    }

    // philaris: here there is significant difference of Linf
    // with relation to L2

    // q is on the Linf convex hull with neighbors (in this hull)
    // s and r

    // here q is a point and s is different from r

    if (sgn == NEGATIVE) {
      CGAL_SDG_DEBUG(std::cout <<
          "debug infinite-edge-int-cf special NEG" << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug infcf special (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << std::endl;);

      if (s.is_point() and r.is_point()) {
        if ((bounded_side_of_bbox(
                q.point(), s.point(), t.point()) == ON_BOUNDED_SIDE) and
            (bounded_side_of_bbox(
                q.point(), r.point(), t.point()) == ON_BOUNDED_SIDE))
        {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << false << std::endl;);
          return false;
          // otherwise, it will return true later
        }
      } else {
        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << true << std::endl;);
        return true;
      }
    } else if (sgn == POSITIVE) {
      // philaris: I assert that s and r are either both points
      //           or both segments
      CGAL_assertion( (s.is_point() and r.is_point())    or
                      (s.is_segment() and r.is_segment())  ) ;
      CGAL_SDG_DEBUG(
          std::cout << "debug infinite-edge-int-cf special POS"
          << std::endl;);
      if (s.is_point() and r.is_point()) {
        // philaris: I assert that sqt and rqt are monotone

        CGAL_assertion(
            or_linf(s.point(), q.point(), t.point()) ==
            DEGENERATE );

        CGAL_assertion(
            or_linf(r.point(), q.point(), t.point()) ==
            DEGENERATE );

        if ((bounded_side_of_bbox(
               t.point(), s.point(), q.point()) ==
                 ON_BOUNDED_SIDE) and
            (bounded_side_of_bbox(
               t.point(), r.point(), q.point()) ==
                 ON_BOUNDED_SIDE))
        {
          CGAL_SDG_DEBUG(
              std::cout
              << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
              << q << ' ' << s << ' ' << r << ' ' << t
              << ' ' << sgn << " returns "
              << true << std::endl;);
          return true;
          // otherwise it will return false later
        }
      } else {
        // here s and r are both segments

        // philaris: I assert that q is an endpoint of both s and r

        CGAL_assertion(
            same_points(q, s.source_site()) or
            same_points(q, s.target_site())
            );

        CGAL_assertion(
            same_points(q, r.source_site()) or
            same_points(q, r.target_site())
            );

        // Since s and r are different, the infinite vertices
        // (q,s,inf) and (q,inf,r) of the Voronoi diagram
        // have the following property:
        // if they are connected, their line has +-45 deg slope
        // and this line goes through point q.
        // This line is Linf-perpendicular to both s and r.
        // Since q is not in conflict with any of (q,s,inf),
        // (q,inf,r), it is enough to check if t is in the
        // infinite box with corner q containing neither
        // of these: (q,s,inf), (q,inf,r), any point of s, r.

        Point_2 otherpnt =
          (same_points(q, s.source_site())) ?
          s.segment().target():
          s.segment().source();

        CGAL_assertion(
            or_linf(t.point(), q.point(), otherpnt) ==
            DEGENERATE );

        Bounded_side bside =
          bounded_side_of_bbox(t.point(), otherpnt, q.point());

        CGAL_assertion(bside != ON_BOUNDARY);

        CGAL_SDG_DEBUG(
            std::cout
            << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
            << q << ' ' << s << ' ' << r << ' ' << t
            << ' ' << sgn << " returns "
            << (bside == ON_BOUNDED_SIDE) << std::endl;);
        return (bside == ON_BOUNDED_SIDE);

      }
    }
    // here it might be sgn == ZERO

    CGAL_SDG_DEBUG(
        std::cout
        << "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
        << q << ' ' << s << ' ' << r << ' ' << t
        << ' ' << sgn << " returns "
        << ( sgn == NEGATIVE ) << std::endl;);
    return ( sgn == NEGATIVE );
  }

};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2_hv

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_HV_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
