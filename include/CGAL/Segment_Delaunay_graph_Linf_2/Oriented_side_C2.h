
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

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

  typedef Are_same_points_C2<K>    Are_same_points_2;
  Are_same_points_2                same_points;

  typedef Are_same_segments_C2<K>  Are_same_segments_2;
  Are_same_segments_2              same_segments;

  using Base::compute_supporting_line;
  using Base::compute_linf_perpendicular;
  using Base::oriented_side_of_line;
  using Base::compute_vertical_projection;
  using Base::compute_horizontal_projection;

public:
  typedef typename Base::Oriented_side        Oriented_side;
  typedef Oriented_side                       result_type;
  typedef Site_2                              argument_type;

  // philaris: not used in Linf
  // computes the oriented side of the point q
  // wrt the line that is passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& q, 
			   const Site_2& s, const Site_2& p) const
  {
    // philaris: q might also be a segment in Linf
    CGAL_precondition( q.is_point() );
    CGAL_precondition( s.is_segment() && p.is_point() );

    Line_2 l = compute_supporting_line( s );
    Line_2 lp = compute_linf_perpendicular(l, p.point());

    Oriented_side retval = lp.oriented_side(q.point());

    std::cout << "debug: Oriented_side_C2 (qsp)= ("
              << q << ") (" << s << ") (" << p << ") "
              << "returns " << retval 
              << std::endl;

    return retval;
  }

  // computes the oriented side of the Voronoi vertex of s1, s2, s3
  // wrt the line that passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
			   const Site_2& s3,
			   const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    Voronoi_vertex_2 v(s1, s2, s3);
    Line_2 l = compute_supporting_line( s );
    Line_2 lp = compute_linf_perpendicular(l, p.point());

    Oriented_side retval = v.oriented_side(lp);

    std::cout << "debug: Oriented_side_C2 (s1,s2,s3,s,p)= ("
              << s1 << ") (" << s2 << ") (" << s3 << ") ("
              << s << ") (" << p << ") "
              << "returns " << retval 
              << std::endl;

    return retval;
  }

  // computes the oriented side of the Voronoi vertex of s1, s2, inf
  // wrt the line that passes through the point p and its direction
  // is the direction of the supporting line of s, rotated by 90
  // degrees counterclockwise.
  Oriented_side operator()(const Site_2& s1, const Site_2& s2,
			   const Site_2& s, const Site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    Line_2 lseg = compute_supporting_line( s );
    Line_2 lp = compute_linf_perpendicular(lseg, p.point());

    // Voronoi_vertex_2 v(s1, s2, inf);
    // compute linf projection of v(s1, s2, inf) on segment s;

    Point_2 proj_of_infv;

    bool is_s1_segment = s1.is_segment();
    bool is_s2_segment = s2.is_segment();

    CGAL_assertion(
        (is_s1_segment and same_segments(s, s1)) or 
        (is_s2_segment and same_segments(s, s2))   );

    if (is_s1_segment and is_s2_segment) {
      // the two segments must have a common endpoint,
      // which is the linf projection

      CGAL_assertion(
          s1.segment().is_horizontal() or
          s1.segment().is_vertical()   or
          s2.segment().is_horizontal() or
          s2.segment().is_vertical()      );

      if (same_points(s1.source_site(), s2.source_site()) or
          same_points(s1.source_site(), s2.target_site())   ) {
        proj_of_infv = s1.source_site().point();
      } else {
        CGAL_assertion(
          same_points(s1.target_site(), s2.source_site()) or
          same_points(s1.target_site(), s2.target_site())   );
        proj_of_infv = s1.target_site().point();
      }

    } else {
      // here, there is a point and a segment in {s1, s2}

      if ( 
           (is_s1_segment and 
            ( same_points(s2, s1.source_site()) or 
              same_points(s2, s1.target_site())   ) )
           or 
           (is_s2_segment and 
            ( same_points(s1, s2.source_site()) or 
              same_points(s1, s2.target_site())   ) )
         )
      { // here the point in {s1,s2} 
        // is endpoint of the segment in {s1,s2}
        if (is_s1_segment) {
          proj_of_infv = s2.point();
        } else {
          CGAL_assertion(is_s2_segment);
          proj_of_infv = s1.point();
        }
      } // end of case: point is endpoint of segment
      else {
        // here, the point is not endpoint of the segment

        CGAL_assertion(not (s.segment().is_horizontal() or
                            s.segment().is_vertical()     ) );

        bool has_lseg_neg_slope =
          CGAL::sign(lseg.a()) == CGAL::sign(lseg.b());

        if (has_lseg_neg_slope) {
          if (is_s1_segment) {
            proj_of_infv = 
              compute_horizontal_projection(lseg, s2.point());
          } else {
            proj_of_infv = 
              compute_vertical_projection(lseg, s1.point());
          }
        } else {
          // here, segment has positive slope
          if (is_s1_segment) {
            proj_of_infv = 
              compute_vertical_projection(lseg, s2.point());
          } else {
            proj_of_infv = 
              compute_horizontal_projection(lseg, s1.point());
          }
        } // end of case: seg has positive slope
      } // end of case: point is not endpoint of segment

    } // end of case: a point and a segment in {s1, s2}


    Oriented_side retval = 
      oriented_side_of_line(lp, proj_of_infv);

    std::cout << "debug: Oriented_side_C2 (s1,s2,s,p)= ("
              << s1 << ") (" << s2 << ") (" 
              << s << ") (" << p << ") "
              << "returns " << retval 
              << std::endl;

    return retval;
  }
};


//-----------------------------------------------------------------------------

} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_ORIENTED_SIDE_C2_H
