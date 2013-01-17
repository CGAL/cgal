
#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Are_same_segments_C2.h>
#include <CGAL/Orientation_Linf_2.h>

namespace CGAL {
  
  namespace SegmentDelaunayGraphLinf_2 {
    
    //-----------------------------------------------------------------------------
    
    template<class K, class Method_tag>
    class Infinite_edge_interior_conflict_C2
    : public Basic_predicates_C2<K>
    {
    public:
    
    typedef Basic_predicates_C2<K>              Base;
    
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
    typedef Are_same_points_C2<K>        Are_same_points_2;
    typedef Are_same_segments_C2<K>      Are_same_segments_2;
    
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

      CGAL_SDG_DEBUG(std::cout <<
          "debug infinite-edge-int-cf entering (q,s,r,t,sgn)= "
          << q << ' ' << s << ' ' << r << ' ' << t
          << ' ' << sgn << std::endl;);

      if ( t.is_segment() ) {
        //Sandeep: This is because end points of a segment are inserted first
        CGAL_SDG_DEBUG(std::cout <<
              "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                     << q << ' ' << s << ' ' << r << ' ' << t
                     << ' ' << sgn << " returns " << false <<
                     std::endl;);
        return false;
      }


      // here and below, t is always a point

      if ( q.is_segment() ) {

        // philaris: difference from L2 here;
        // in L2, r and s are endpoints of q
        // in Linf they still have to be points, but
        // they do not have to be endpoints of q
        // (this has to be checked)
        // Sandeep: the following are the examples
        //
        //         r.
        // _____________________ s
        //         q
        //     s.
        //  r____________________
        //           q
        //
        CGAL_assertion(s.is_point() and r.is_point());

        CGAL_assertion( q.segment().is_horizontal() or q.segment().is_vertical() );

        bool is_s_endp_of_q =
        same_points(s, q.source_site()) or
        same_points(s, q.target_site());
        bool is_r_endp_of_q =
        same_points(r, q.source_site()) or
        same_points(r, q.target_site());

        if (is_s_endp_of_q and is_r_endp_of_q) {
          CGAL_SDG_DEBUG(std::cout << "sandeep: debug infcf s,r are end points of q" << std::endl;);
          CGAL_SDG_DEBUG(std::cout <<
              "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
                     << q << ' ' << s << ' ' << r << ' ' << t
                     << ' ' << sgn << " returns " << (sgn==POSITIVE) <<
                     std::endl;);
          return (sgn == POSITIVE);
        }

      }
      if ( s.is_point() && r.is_point() && same_points(s, r) ) {
        // MK::ERROR: write this code using the compare_x_2 and
        //    compare_y_2 predicates instead of computing the inner
        //    product...
        // philaris: adaptation to Linf
        
        CGAL_SDG_DEBUG(std::cout << "debug infcf s=r and is point" << std::endl;);
        
        Comparison_result cmpxst = cmpx(s.point(), t.point());
        Comparison_result cmpxtq = cmpx(t.point(), q.point());
        Comparison_result cmpyst = cmpy(s.point(), t.point());
        Comparison_result cmpytq = cmpy(t.point(), q.point());
        
        //Sign sgn1 = -sign_of( cmpxst * cmpxtq + cmpyst * cmpytq );
        Sign sgn1 = CGAL::compare(0, cmpxst * cmpxtq + cmpyst * cmpytq);

        CGAL_assertion( sgn1 != ZERO );
        CGAL_SDG_DEBUG(std::cout << "sandeep: debug infcf about to return " << (sgn1 == POSITIVE) << std::endl;);
        CGAL_SDG_DEBUG(std::cout <<
            "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
               << q << ' ' << s << ' ' << r << ' ' << t
               << ' ' << sgn << " returns " << (sgn1==POSITIVE)
               << std::endl;);
        return (sgn1 == POSITIVE);
      }
      CGAL_SDG_DEBUG(std::cout <<
          "sandeep: debug infcf about to return "
          << (sgn == NEGATIVE) << std::endl;);
      CGAL_SDG_DEBUG(std::cout <<
          "debug infinite-edge-int-cf with (q,s,r,t,sgn)= "
             << q << ' ' << s << ' ' << r << ' ' << t
             << ' ' << sgn << " returns " << (sgn==NEGATIVE)
             << std::endl;);
      return ( sgn == NEGATIVE );
      }

    };
    
    
    //-----------------------------------------------------------------------------
    
  } //namespace SegmentDelaunayGraphLinf_2
  
} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_INFINITE_EDGE_INTERIOR_CONFLICT_C2_H
