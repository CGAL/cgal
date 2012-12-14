#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H


#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Compare_x_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Compare_y_2.h>
#include <CGAL/Side_of_oriented_square_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {


template<class K>
class Voronoi_vertex_ring_C2
  : public Basic_predicates_C2<K>
{
public:
  typedef Basic_predicates_C2<K> Base;

  typedef enum {PPP = 0, PPS, PSS, SSS} vertex_t;
  struct PPP_Type {};
  struct PPS_Type {};
  struct PSS_Type {};
  struct SSS_Type {};

  typedef typename Base::Point_2             Point_2;
  typedef typename Base::Segment_2           Segment_2;
  typedef typename Base::Line_2              Line_2;
  typedef typename Base::Site_2              Site_2;
  typedef typename Base::FT                  FT;
  typedef typename Base::RT                  RT;

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Sign                Sign;

  typedef typename Base::Polychainline_2     Polychainline_2;

  using Base::compute_supporting_line;
  using Base::compute_linf_projection_hom;
  using Base::oriented_side_of_line;
  using Base::opposite_line;
  using Base::to_ft;
  using Base::compute_line_from_to;
  using Base::compute_horizontal_projection;
  using Base::compute_vertical_projection;
  using Base::compute_linf_perpendicular;

private:
  typedef Are_same_points_C2<K>          Are_same_points_2;
  typedef Are_same_segments_C2<K>        Are_same_segments_2;
  typedef Side_of_oriented_square_2<K>   Side_of_oriented_square_2_Type;
  typedef Bisector_Linf<K>               Bisector_Linf_Type;
  typedef Compare_x_2<K>                 Compare_x_2_Sites_Type;
  typedef Compare_y_2<K>                 Compare_y_2_Sites_Type;


  typedef typename K::Intersections_tag ITag;

  Are_same_points_2                same_points;
  Are_same_segments_2              same_segments;
  Side_of_oriented_square_2_Type   side_of_oriented_square;
  Bisector_Linf_Type               bisector_linf;
  Compare_x_2_Sites_Type           scmpx;
  Compare_y_2_Sites_Type           scmpy;

private:
  //--------------------------------------------------------------------------

  void
  compute_ppp(const Site_2& sp, const Site_2& sq, const Site_2& sr)
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug vring entering compute_ppp" << std::endl;);

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

    v_type = PPP;

    RT x_min, x_max, y_min, y_max;
    RT x_center, y_center;
    RT half(0.5);
    RT two(2);

    bool is_set_x_center(false);
    bool is_set_y_center(false);
    bool is_set_x_max(false);
    bool is_set_y_max(false);
    bool is_set_x_min(false);
    bool is_set_y_min(false);

    Comparison_result cmpxqp = CGAL::compare(q.x(), p.x());

    if (cmpxqp == SMALLER) { // q.x() < p.x()
      x_min = q.x();
      x_max = p.x();
    } else if (cmpxqp == LARGER) { // q.x() > p.x()
      x_min = p.x();
      x_max = q.x();
    } else { // q.x() = p.x()
      x_min = p.x();
      x_max = p.x();
      y_center = half * (p.y() + q.y());
      is_set_y_center = true;

      CGAL_SDG_DEBUG(std::cout << "debug set y_center=" << 
        y_center << std::endl;);

      Comparison_result cmpxrothers = CGAL::compare(r.x(), p.x());
      if (cmpxrothers == SMALLER) {
        CGAL_SDG_DEBUG(std::cout << "debug r is left of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  and (cmpyrq == LARGER)) or
            ((cmpyrp == SMALLER) and (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two*y_center - r.y();
            is_set_y_min = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_min=" << 
              y_min << std::endl;);
          } else {
            y_max = two*y_center - r.y();
            is_set_y_max = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_max=" << 
              y_max << std::endl;);
          }
        } 
      } else if (cmpxrothers == LARGER) {
        CGAL_SDG_DEBUG(std::cout << "debug r is right of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  and (cmpyrq == LARGER)) or
            ((cmpyrp == SMALLER) and (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two*y_center - r.y();
            is_set_y_min = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_min=" << 
              y_min << std::endl;);
          } else {
            y_max = two*y_center - r.y();
            is_set_y_max = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_max=" << 
              y_max << std::endl;);
          }
        } 
      } else {
        // not possible
      }
    }

    Comparison_result cmpyqp = CGAL::compare(q.y(), p.y());

    if (cmpyqp == SMALLER) { // q.y() < p.y()
      if (not is_set_y_min) {
        y_min = q.y();
      }
      if (not is_set_y_max) {
        y_max = p.y();
      }
    } else if (cmpyqp == LARGER) { // q.y() > p.y()
      if (not is_set_y_min) {
        y_min = p.y();
      }
      if (not is_set_y_max) {
        y_max = q.y();
      }
    } else { //  q.y() = p.y()
      if (not is_set_y_min) {
        y_min = p.y();
      }
      if (not is_set_y_max) {
        y_max = p.y();
      }
      x_center = half * (p.x() + q.x());
      is_set_x_center = true;

      Comparison_result cmpyrothers = CGAL::compare(r.y(), p.y());
      if (cmpyrothers == SMALLER) {
        Comparison_result cmpxrp = CGAL::compare(r.x(), p.x());
        Comparison_result cmpxrq = CGAL::compare(r.x(), q.x());
        if (((cmpxrp == LARGER)  and (cmpxrq == LARGER)) or
            ((cmpxrp == SMALLER) and (cmpxrq == SMALLER))
           ) {
          // do fix
          if (cmpxrp == LARGER) {
            x_min = two*x_center - r.x();
            is_set_x_min = true;
          } else {
            x_max = two*x_center - r.x();
            is_set_x_max = true;
          }
        } 
      } else if (cmpyrothers == LARGER) {
        Comparison_result cmpxrp = CGAL::compare(r.x(), p.x());
        Comparison_result cmpxrq = CGAL::compare(r.x(), q.x());
        if (((cmpxrp == LARGER)  and (cmpxrq == LARGER)) or
            ((cmpxrp == SMALLER) and (cmpxrq == SMALLER))
           ) {
          // do fix
          if (cmpxrp == LARGER) {
            x_min = two*x_center - r.x();
            is_set_x_min = true;
          } else {
            x_max = two*x_center - r.x();
            is_set_x_max = true;
          }
        } 
      } else {
        // not possible
      }

    }

    Comparison_result cmpxrmin = CGAL::compare(r.x(), x_min);
    Comparison_result cmpxrmax = CGAL::compare(r.x(), x_max);
    if (cmpxrmin == SMALLER) { 
	// here r.x() < x_min <= x_max
        if (not is_set_x_min) {
          x_min = r.x();
        }
    } else if (cmpxrmin == LARGER) {
      // here r.x() > x_min
      if (cmpxrmax == LARGER) {
        // here x_min <= x_max < r.x()
        if (not is_set_x_max) {
          x_max = r.x();
        }
      } else if (cmpxrmax == SMALLER) {
        // x_min < r.x() < x_max
        // do nothing
      } else { // r.x() = x_max
        // r.x() = p.x() or r.x() = q.x()
        if (CGAL::compare(r.x(), p.x()) == EQUAL) {
          y_center = half * (p.y() + r.y());
          //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
          Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
          if ((cmpyqp == LARGER) and (cmpyqr == LARGER)) {
            y_min = two*y_center - q.y();
            is_set_y_min = true;
          }
          if ((cmpyqp == SMALLER) and (cmpyqr == SMALLER)) {
            y_max = two*y_center - q.y();
            is_set_y_max = true;
          }
        } else {
          y_center = half * (q.y() + r.y());
          Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
          Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
          if ((cmpypq == LARGER) and (cmpypr == LARGER)) {
            y_min = two*y_center - p.y();
            is_set_y_min = true;
          }
          if ((cmpypq == SMALLER) and (cmpypr == SMALLER)) {
            y_max = two*y_center - p.y();
            is_set_y_max = true;
          }
        }
        is_set_y_center = true;
      }
    } else {
      // here r.x() = x_min
      // r.x() = p.x() or r.x() = q.x()
      if (CGAL::compare(r.x(), p.x()) == EQUAL) {
        CGAL_SDG_DEBUG(std::cout << "debug r.x = p.x" << std::endl;);
        // r.x() = p.x()
        y_center = half * (p.y() + r.y());
        //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
        Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
        if ((cmpyqp == LARGER) and (cmpyqr == LARGER)) {
          CGAL_SDG_DEBUG(std::cout << "debug q is above p, r" << std::endl;);
          y_min = two*y_center - q.y();
          is_set_y_min = true;
        }
        if ((cmpyqp == SMALLER) and (cmpyqr == SMALLER)) {
          CGAL_SDG_DEBUG(std::cout << "debug q is below p, r" << std::endl;);
          y_max = two*y_center - q.y();
          is_set_y_max = true;
        }
      } else { 
        // r.x() = q.x()
        CGAL_SDG_DEBUG(std::cout << "debug r.x = q.x" << std::endl;);
        y_center = half * (q.y() + r.y());
        Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
        Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
        if ((cmpypq == LARGER) and (cmpypr == LARGER)) {
          CGAL_SDG_DEBUG(std::cout << "debug p is above q, r" << std::endl;);
          y_min = two*y_center - p.y();
          is_set_y_min = true;
        }
        if ((cmpypq == SMALLER) and (cmpypr == SMALLER)) {
          CGAL_SDG_DEBUG(std::cout << "debug p is below q, r" << std::endl;);
          y_max = two*y_center - p.y();
          is_set_y_max = true;
        }
      }
      is_set_y_center = true;
    }

    Comparison_result cmpyrmin = CGAL::compare(r.y(), y_min);
    Comparison_result cmpyrmax = CGAL::compare(r.y(), y_max);
    if (cmpyrmin == SMALLER) { 
      // here r.y() < y_min <= y_max
      if (not is_set_y_min) { 
        y_min = r.y();
      }
    } else if (cmpyrmin == LARGER) {
      // here r.y() > y_min
      if (cmpyrmax == LARGER) {
        // here y_min <= y_max < r.y()
        if (not is_set_y_max) { 
          y_max = r.y();
        }
      } else if (cmpyrmax == SMALLER) {
        // y_min < r.y() < y_max
        // do nothing
      } else { // r.y() = y_max
        // r.y() = p.y() or r.y() = q.y()
        if (CGAL::compare(r.y(), p.y()) == EQUAL) {
          x_center = half * (p.x() + r.x());
          //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
          Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
          if ((cmpxqp == LARGER) and (cmpxqr == LARGER)) {
            x_min = two*x_center - q.x();
            is_set_x_min = true;
          }
          if ((cmpxqp == SMALLER) and (cmpxqr == SMALLER)) {
            x_max = two*x_center - q.x();
            is_set_x_max = true;
          }
        } else {
          x_center = half * (q.x() + r.x());
          Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
          Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
          if ((cmpxpq == LARGER) and (cmpxpr == LARGER)) {
            x_min = two*x_center - p.x();
            is_set_x_min = true;
          }
          if ((cmpxpq == SMALLER) and (cmpxpr == SMALLER)) {
            x_max = two*x_center - p.x();
            is_set_x_max = true;
          }
        }
        is_set_x_center = true;
      }
    } else {
      // here r.y() = y_min
      // r.y() = p.y() or r.y() = q.y()
      if (CGAL::compare(r.y(), p.y()) == EQUAL) {
        x_center = half * (p.x() + r.x());
        //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
        Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
        if ((cmpxqp == LARGER) and (cmpxqr == LARGER)) {
          x_min = two*x_center - q.x();
          is_set_x_min = true;
        }
        if ((cmpxqp == SMALLER) and (cmpxqr == SMALLER)) {
          x_max = two*x_center - q.x();
          is_set_x_max = true;
        }
      } else {
        x_center = half * (q.x() + r.x());
        Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
        Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
        if ((cmpxpq == LARGER) and (cmpxpr == LARGER)) {
          x_min = two*x_center - p.x();
          is_set_x_min = true;
        }
        if ((cmpxpq == SMALLER) and (cmpxpr == SMALLER)) {
          x_max = two*x_center - p.x();
          is_set_x_max = true;
        }
      }
      is_set_x_center = true;
    }

    Comparison_result cmpsides = 
	CGAL::compare(x_max - x_min, y_max - y_min);

    // if bounding box is non-square and points are not 
    // on corners of it, then grow it to become square
    switch(cmpsides) {
      case SMALLER:
        CGAL_SDG_DEBUG(std::cout << "rectangle has to be made fatter" << std::endl;);
        // make rectangle fatter
        if (is_set_x_center) {
          CGAL_SDG_DEBUG(std::cout << "x_center already set" << std::endl;);
          // grow in both sides
          break;
        }
        // grow only if any point is inside vertical sides
	if (((CGAL::compare(p.x(), x_min) == EQUAL)   and
	     (CGAL::compare(p.y(), y_max) == SMALLER) and
	     (CGAL::compare(p.y(), y_min) == LARGER)     ) or 
	    ((CGAL::compare(q.x(), x_min) == EQUAL)   and
	     (CGAL::compare(q.y(), y_max) == SMALLER) and
	     (CGAL::compare(q.y(), y_min) == LARGER)     ) or 
	    ((CGAL::compare(r.x(), x_min) == EQUAL)   and
	     (CGAL::compare(r.y(), y_max) == SMALLER) and
	     (CGAL::compare(r.y(), y_min) == LARGER)     )   )
        { // grow rectangle to the right
          CGAL_SDG_DEBUG(std::cout << "debug vring grow right" << std::endl;);
          x_max = x_min + y_max - y_min;
        } else 
        { // grow rectangle to the left
          CGAL_SDG_DEBUG(std::cout << "debug vring grow left" << std::endl;);
          x_min = x_max - y_max + y_min;
        }
        break;
      case LARGER:
        CGAL_SDG_DEBUG(std::cout << "rectangle has to be made taller" << std::endl;);
        // make rectangle taller
        if (is_set_y_center) {
          // grow in both sides
          CGAL_SDG_DEBUG(std::cout << "y_center already set" << std::endl;);
          break;
        }
        // grow only if any point is inside horizontal sides
	if (((CGAL::compare(p.y(), y_min) == EQUAL)   and
	     (CGAL::compare(p.x(), x_max) == SMALLER) and
	     (CGAL::compare(p.x(), x_min) == LARGER)     ) or 
	    ((CGAL::compare(q.y(), y_min) == EQUAL)   and
	     (CGAL::compare(q.x(), x_max) == SMALLER) and
	     (CGAL::compare(q.x(), x_min) == LARGER)     ) or 
	    ((CGAL::compare(r.y(), y_min) == EQUAL)   and
	     (CGAL::compare(r.x(), x_max) == SMALLER) and
	     (CGAL::compare(r.x(), x_min) == LARGER)     )   )
        { // grow rectangle upwards
          CGAL_SDG_DEBUG(std::cout << "debug vring grow upwards" << std::endl;);
          y_max = y_min + x_max - x_min;
        } else 
        { // grow rectangle downwards
          CGAL_SDG_DEBUG(std::cout << "debug vring grow downwards" << std::endl;);
          y_min = y_max - x_max + x_min;
        }
        break;
      case EQUAL:
        // do nothing
        break;
    }

    ux_ = x_min + x_max;
    uy_ = y_min + y_max;
    uz_ = RT(2) ;
  }

  //--------------------------------------------------------------------------

  void
  compute_pss(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_point() && q.is_segment() &&
		       r.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug: compute_pss entering p=" << p 
       << " q=" << q << " r=" << r << std::endl;);

    v_type = PSS;

    bool pq =
      same_points(p, q.source_site()) || same_points(p, q.target_site());
    bool pr =
      same_points(p, r.source_site()) || same_points(p, r.target_site());

    Point_2 pp = p.point();

    if ( pq && pr ) {
      // philaris: result should be point p

      ux_ = pp.x();
      uy_ = pp.y();
      uz_ = RT(1);
      return;
    }
    //sandeep: We will find the vertex by comparisons and simple arithmatic
    
    Point_2 qsrc = q.source_site().point();
    Point_2 rsrc = r.source_site().point();
  
    Segment_2 qseg = q.segment();
    Segment_2 rseg = r.segment();
    //case1:when q and r are both horizontal
    if (qseg.is_horizontal() and rseg.is_horizontal()) {
      //if (CGAL::compare(q.source_site().y(),r.source_site().y()) == LARGER) 
      ux_ = RT(2) * pp.x() - qsrc.y() + rsrc.y();
      uy_ = qsrc.y() + rsrc.y();
      uz_ = RT(2);
      return;
    }
    //case2:when q and r are both vertical
    if (qseg.is_vertical() and rseg.is_vertical()) {
      ux_ = qsrc.x() + rsrc.x();
      uy_ = RT(2) * pp.y() + qsrc.x() - rsrc.x();
      uz_ = RT(2);
      return;
    }
    //case3:when q is horizontal and r is vertical
    if (qseg.is_horizontal() and rseg.is_vertical()) {
      if (CGAL::compare( abs(pp.x() - rsrc.x()),
                         abs(pp.y() - qsrc.y()) )
                         == SMALLER ) {
        //diff x is smaller than diff y
        ux_ = RT(2) * rsrc.x() + qsrc.y() - pp.y();
        uy_ = qsrc.y() + pp.y();
        uz_ = RT(2);
      } else {//diff x is larger or equal to diff y
        ux_ = pp.x() + rsrc.x();
        uy_ = RT(2) * qsrc.y() + rsrc.x() - pp.x();
        uz_ = RT(2);
      }
      return;
    }
    //case4:when q is vertical and r is horizontal
    if (CGAL::compare( abs(pp.x() - qsrc.x()),
                       abs(pp.y() - rsrc.y()) )
                       == SMALLER ) {
      //diff x is smaller than diff y
      ux_ = RT(2) * qsrc.x() - rsrc.y() + pp.y();
      uy_ = rsrc.y() + pp.y();
      uz_ = RT(2);
    } else {//diff x is larger or equal to diff y
      ux_ = pp.x() + qsrc.x();
      uy_ = RT(2) * rsrc.y() - qsrc.x() + pp.x();
      uz_ = RT(2);
    }

    return;
  }



  //--------------------------------------------------------------------------


  void
  compute_pps(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_point() && q.is_point() &&
		                   r.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug: compute_pps entering p=" << p 
      << " q=" << q << " r=" << r << std::endl;);

    v_type = PPS;

    bool p_endp_r = is_endpoint_of(p, r);
    bool q_endp_r = is_endpoint_of(q, r);
    //incase p or q is end point of r no square is possible in ccw sense
    CGAL_assertion(not (p_endp_r or q_endp_r));

    bool samexpq = (scmpx(p, q) == EQUAL);
    bool sameypq = (scmpy(p, q) == EQUAL);

    bool samecoordpq = samexpq or sameypq ;
   
    Point_2 pp = p.point(), qp = q.point();
    Segment_2 rs = r.segment();
    //Sandeep: We divide the computation in two cases
    //case 1: r is horizontal
    //Case 2: r is vertical
    if (rs.is_horizontal()) {
      if( CGAL::compare( abs(pp.x() - qp.x()),
                         abs(pp.y() - rs.source().y()) )
                         == LARGER &&
          CGAL::compare( abs(pp.x() - qp.x()),
                         abs(qp.y() - rs.source().y()) )
                         == LARGER ) {
      // This is the case when diffx is greater than diffy
      
        ux_ = pp.x() + qp.x();
        uy_ = RT(2) * rs.source().y() + pp.x() - qp.x();
        uz_ = RT(2);
      } else {
        // diffy is greater
        if (sameypq) {
          ux_ = pp.x() + qp.x();
          uy_ = pp.y() + rs.source().y();
          uz_ = RT(2);
        } else if ( CGAL::compare( abs(pp.y() - rs.source().y()),
                                   abs(qp.y() - rs.source().y()) )
                                   == LARGER) {
          // p is more distant than q from r 
          ux_ = RT(2) * qp.x() + pp.y() - rs.source().y();
          uy_ = pp.y() + rs.source().y();
          uz_ = RT(2);
          
        } else {
          // q is more distant than p from r
          ux_ = RT(2) * pp.x() - qp.y() + rs.source().y();
          uy_ = qp.y() + rs.source().y();
          uz_ = RT(2);
          
        }
      }// end of diffy greater
    
    } else {//rs is vertical
      if( CGAL::compare( abs(pp.y() - qp.y()),
                         abs(pp.x() - rs.source().x()) )
                         == LARGER &&
          CGAL::compare( abs(pp.y() - qp.y()),
                         abs(qp.x() - rs.source().x()) )
                         == LARGER ) {
        // This is the case when diffy is greater than diffx
        
        ux_ = RT(2) * rs.source().x() + pp.y() - qp.y();
        uy_ = pp.y() + qp.y();
        uz_ = RT(2);
      } else {
        // diffx is greater
        if (samexpq) {
          ux_ = pp.x() + rs.source().x();
          uy_ = pp.y() + qp.y();
          uz_ = RT(2);
        } else if ( CGAL::compare( abs(pp.x() - rs.source().x()),
                                   abs(qp.x() - rs.source().x()) )
                                   == LARGER) {
          // p is more distant than q from r
          ux_ = pp.x() + rs.source().x();
          uy_ = RT(2) * qp.y() + pp.x() - rs.source().x();
          uz_ = RT(2);
          
        } else {
          // q is more distant than p from r
          ux_ = qp.x() + rs.source().x();
          uy_ = RT(2) * pp.y() - qp.x() + rs.source().x();
          uz_ = RT(2);
          
        }
      }// end of diffx greater
    }//end of vertical rs
    return;
  }

  //--------------------------------------------------------------------------

  bool check_if_exact(const Site_2& , unsigned int ,
		      const Tag_false&) const
  {
    return true;
  }

  bool check_if_exact(const Site_2& s, unsigned int i,
		      const Tag_true&) const
  {
    return s.is_input(i);
  }

  // determines of the segment s is on the positive halfspace as
  // defined by the supporting line of the segment supp; the line l
  // is supposed to be the supporting line of the segment supp and we
  // pass it so that we do not have to recompute it
  bool
  is_on_positive_halfspace(const Site_2& supp,
			   const Site_2& s, const Line_2& l) const
  {
    CGAL_precondition( supp.is_segment() && s.is_segment() );

    if ( same_segments(supp.supporting_site(),
		       s.supporting_site()) ) {
      return false;
    }

    if ( same_points(supp.source_site(), s.source_site()) ||
	 same_points(supp.target_site(), s.source_site()) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( same_points(supp.source_site(), s.target_site()) ||
	 same_points(supp.target_site(), s.target_site()) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    ITag itag;

    if ( !check_if_exact(s, 0, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(0)) ) {
      return oriented_side_of_line(l, s.target()) == ON_POSITIVE_SIDE;
    }

    if ( !check_if_exact(s, 1, itag) &&
	 same_segments(supp.supporting_site(),
		       s.crossing_site(1)) ) {
      return oriented_side_of_line(l, s.source()) == ON_POSITIVE_SIDE;
    }

    return Base::is_on_positive_halfspace(l, s.segment());
  }

  //--------------------------------------------------------------------------

  void
  compute_sss(const Site_2& p, const Site_2& q, const Site_2& r)
  {
    CGAL_precondition( p.is_segment() && q.is_segment() &&
		       r.is_segment() );

    v_type = SSS;
    
    int numh=0, numv=0;
    Segment_2 ps = p.segment(), qs = q.segment(), rs = r.segment();
    
    (ps.is_horizontal()) ? numh++ : numv++;
    (qs.is_horizontal()) ? numh++ : numv++;
    (rs.is_horizontal()) ? numh++ : numv++;
  
    CGAL_assertion(numh == 2 or numv == 2);
    
    if (numv == 2) { // two vertical and one horizontal segment
      if (ps.is_horizontal()) {
        ux_ = qs.source().x() + rs.source().x();
        uy_ = RT(2) * ps.source().y() + qs.source().x() - rs.source().x();
        uz_ = RT(2);
      }
      if (qs.is_horizontal()) {
        ux_ = ps.source().x() + rs.source().x();
        uy_ = RT(2) * qs.source().y() + rs.source().x() - ps.source().x();
        uz_ = RT(2);
      }
      if (rs.is_horizontal()) {
        ux_ = qs.source().x() + ps.source().x();
        uy_ = RT(2) * rs.source().y() + ps.source().x() - qs.source().x();
        uz_ = RT(2);
      }
     
    } else { // two horizontal and one vertical segment
      if (ps.is_vertical()) {
        uy_ = qs.source().y() + rs.source().y();
        ux_ = RT(2) * ps.source().x() + rs.source().y() - qs.source().y();
        uz_ = RT(2);
      }
      if (qs.is_vertical()) {
        uy_ = ps.source().y() + rs.source().y();
        ux_ = RT(2) * qs.source().x() + ps.source().y() - rs.source().y();
        uz_ = RT(2);
      }
      if (rs.is_vertical()) {
        uy_ = qs.source().y() + ps.source().y();
        ux_ = RT(2) * rs.source().x() + qs.source().y() - ps.source().y();
        uz_ = RT(2);
      }

    }//end of two horizontal and one vertical segment case
  
    return;
  }

  //--------------------------------------------------------------------------

  void
  compute_vertex(const Site_2& s1, const Site_2& s2, const Site_2& s3)
  {
    if ( s1.is_point() && s2.is_point() && s3.is_point() ) {
      compute_ppp(s1, s2, s3);

    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      compute_vertex(s2, s3, s1);
      pps_idx = 1;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
      pps_idx = 2;

    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      compute_pps(s1, s2, s3);
      pps_idx = 0;

    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      compute_pss(s1, s2, s3);
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      compute_vertex(s2, s3, s1);
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      compute_vertex(s3, s1, s2);
    } else {
      compute_sss(s1, s2, s3);
    }

  }

  //--------------------------------------------------------------------------

  bool is_endpoint_of(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    return ( same_points(p, s.source_site()) ||
	     same_points(p, s.target_site()) );
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           the orientation test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line_2& l, PPP_Type) const
  {
    Sign s_uz = CGAL::sign(uz_);
    Sign s_l =
      CGAL::sign(l.a() * ux_ + l.b() * uy_ + l.c() * uz_);

    return s_uz * s_l;
  }

  
  //--------------------------------------------------------------------------

  Orientation
  orientation(const Line_2& l, PPS_Type) const
  {
    Sign s_uz = CGAL::sign(uz_);
    Sign s_l = CGAL::sign(l.a() * ux_ + l.b() * uy_ + l.c() * uz_);

    return s_uz * s_l;
  }


  //--------------------------------------------------------------------------

  // the cases PSS and SSS are identical
  template<class Type>
  Orientation
  orientation(const Line_2& l, Type) const
  {
    Sign s_uz = CGAL::sign(uz_);
    Sign s_l = CGAL::sign(l.a() * ux_ + l.b() * uy_ + l.c() * uz_);

    return s_uz * s_l;
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                              the incircle test
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a point
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Sign check_easy_degeneracies(const Site_2& t, PPS_Type,
			       bool& use_result) const
  {
    CGAL_precondition( t.is_point() );

    use_result = false;
    if (  ( p_.is_point() && same_points(p_, t) ) ||
	  ( q_.is_point() && same_points(q_, t) ) ||
	  ( r_.is_point() && same_points(r_, t) )  ) {
      use_result = true;
      return ZERO;
    }

    // philaris: following might have to be changed
    // philaris: I remove the following line to be on the safe side
    /*
    if (  ( p_.is_segment() && is_endpoint_of(t, p_) ) || 
	  ( q_.is_segment() && is_endpoint_of(t, q_) ) ||
	  ( r_.is_segment() && is_endpoint_of(t, r_) )  ) {
      use_result = true;
      return POSITIVE;
    }
    */

    return ZERO;
  }

  inline
  Sign check_easy_degeneracies(const Site_2& t, PSS_Type,
			       bool& use_result) const
  {
    CGAL_precondition( t.is_point() );
    return check_easy_degeneracies(t, PPS_Type(), use_result);
  }

  inline
  Sign check_easy_degeneracies(const Site_2& t, SSS_Type,
			       bool& use_result) const
  {
    CGAL_precondition( t.is_point() );
    use_result = false;
    // ADD THE CASES WHERE t IS AN ENDPOINT OF ONE OF THE SEGMENTS
    return ZERO;
  }

  //--------------------------------------------------------------------------

  template<class Type>
  inline  
  Sign incircle_p(const Site_2& st, Type type) const
  {
    CGAL_precondition( st.is_point() );

    bool use_result(false);
    Sign s = check_easy_degeneracies(st, type, use_result);
    if ( use_result ) { return s; }

    return incircle_p_no_easy(st, type);
  }

  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& st, PPP_Type) const
  {
    CGAL_precondition( st.is_point() );

    Point_2 t = st.point();

    Point_2 pref = p_.point();
  
    RT dup =
       CGAL::max(CGAL::abs(ux_ - pref.x() * uz_),
                 CGAL::abs(uy_ - pref.y() * uz_));
  
    RT dut =
       CGAL::max(CGAL::abs(ux_ - t.x() * uz_),
                 CGAL::abs(uy_ - t.y() * uz_));
  
    return CGAL::sign(dut - dup);
  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PPS_Type ) const
  {
    CGAL_precondition( st.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_p_no_easy PPS p=" 
      << p_ << " q=" << q_  << " r=" << r_ << " t=" << st 
      << std::endl;);

    Point_2 t = st.point();
  
    Point_2 pref = p_ref().point();
  
    RT dup =
       CGAL::max(CGAL::abs(ux_ - pref.x() * uz_),
                 CGAL::abs(uy_ - pref.y() * uz_));
  
    RT dut =
       CGAL::max(CGAL::abs(ux_ - t.x() * uz_),
                 CGAL::abs(uy_ - t.y() * uz_));
  
   return CGAL::sign(dut - dup);
  

  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PSS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    Point_2 t = st.point();

    Point_2 pref = p_ref().point();

    RT dup =
       CGAL::max(CGAL::abs(ux_ - pref.x() * uz_),
                 CGAL::abs(uy_ - pref.y() * uz_));

    RT dut =
       CGAL::max(CGAL::abs(ux_ - t.x() * uz_),
                 CGAL::abs(uy_ - t.y() * uz_));

    return CGAL::sign(dut - dup);

  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, SSS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    // Sandeep: This is fine, we do not need corner points
    Point_2 t = st.point();
    Segment_2 pseg = p_.segment();
    
    RT dup;
    if (pseg.is_horizontal()) {
      dup = CGAL::abs(uy_ - pseg.source().y() * uz_);
    } else {//p_ is vertical
      dup = CGAL::abs(ux_ - pseg.source().x() * uz_);
    }

    RT dut = 
      CGAL::max(CGAL::abs(ux_ - t.x() * uz_),
                CGAL::abs(uy_ - t.y() * uz_));

    return CGAL::sign(dut - dup);
  }



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& t) const 
  {
    CGAL_SDG_DEBUG(std::cout << "debug: entering vring incircle_p p="
      << p_ << " q=" << q_ << " r=" << r_ << " t=" << t 
      << std::endl;);

    if ( is_degenerate_Voronoi_circle() ) {
      return POSITIVE;
    }

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(t, PPP_Type());
      break;
    case PPS:
      s = incircle_p(t, PPS_Type());
      break;
    case PSS:
      s = incircle_p(t, PSS_Type());
      break;
    case SSS:
      s = incircle_p(t, SSS_Type());
      break;
    }

    return s;
  }

  Sign incircle_p_no_easy(const Site_2& t) const 
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(t, PPP_Type());
      break;
    case PPS:
      s = incircle_p_no_easy(t, PPS_Type());
      break;
    case PSS:
      s = incircle_p_no_easy(t, PSS_Type());
      break;
    case SSS:
      s = incircle_p_no_easy(t, SSS_Type());
      break;
    }

    return s;
  }

  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //  the incircle test when the fourth site is a segment
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------

  Oriented_side
  oriented_side_l2(const Line_2& l, const Point_2& p, PPP_Type) const
  {
    Sign s_uz = CGAL::sign(uz_);

    RT px = uz_ * p.x() - ux_;
    RT py = uz_ * p.y() - uy_;

    Sign s1 = CGAL::sign(l.b() * px - l.a() * py);

    return s_uz * s1;
  }

  Oriented_side
  oriented_side_l2(const Line_2& l, const Point_2& p, PPS_Type) const
  {
    RT dx = ux_ - uz_ * p.x();
    RT dy = uy_ - uz_ * p.y();

    return CGAL::sign(uz_) * CGAL::sign(dy * l.a() - dx * l.b());
  }

  // the cases PSS and SSS are identical
  template<class Type>
  Oriented_side
  oriented_side_l2(const Line_2& l, const Point_2& p, Type) const
  {
    RT px = p.x();
    RT py = p.y();

    RT dx = ux_ - px * uz_;
    RT dy = uy_ - py * uz_;

    RT a = l.a();
    RT b = l.b();

    return CGAL::sign(uz_) * CGAL::sign(a * dy - b * dx);
  }


  // philaris: 
  template<class Type>
  Oriented_side
  oriented_side_linf(const Line_2& l, const Point_2& p, Type) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug oriented_side_linf " << std::endl;);

    Point_2 vv (ux_, uy_, uz_);

    Line_2 l1 = compute_linf_perpendicular(l, vv);

    return oriented_side_of_line(l1, p);
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  template<class Type>
  Sign incircle_s(const Site_2& t, Type type) const
  {
    CGAL_precondition( t.is_segment() );

    if ( v_type == PPP || v_type == PPS ) {
      if (  p_.is_point() && q_.is_point() &&
	    is_endpoint_of(p_, t) && is_endpoint_of(q_, t)  ) {
	return NEGATIVE;
      }

      if (  p_.is_point() && r_.is_point() &&
	    is_endpoint_of(p_, t) && is_endpoint_of(r_, t)  ){
	return NEGATIVE;
      }

      if (  q_.is_point() && r_.is_point() &&
	    is_endpoint_of(q_, t) && is_endpoint_of(r_, t)  ){
	return NEGATIVE;
      }
    }

    if ( v_type == PSS ) {
      if ( p_.is_segment() &&
	   same_segments(p_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
      if ( q_.is_segment() &&
	   same_segments(q_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
      if ( r_.is_segment() &&
	   same_segments(r_.supporting_site(),
			 t.supporting_site()) ) {
	return POSITIVE;
      }
    }

    Sign retval = incircle_s_no_easy(t, type);

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s: about to return retval of"
      << " incircle_s_no_easy = " << retval << std::endl;);

    return retval;
  }

  template<class Type>
  Sign incircle_s_no_easy(const Site_2& t, Type type) const
  {

    CGAL_SDG_DEBUG(std::cout << "debug fn incircle_s_no_easy pqrt= (" << p_ << ") ("
      << q_ << ") (" << r_ << ") (" << t << ")" << std::endl;);

    bool is_p_point = p_.is_point();
    bool is_q_point = q_.is_point();
    bool is_r_point = r_.is_point();

    unsigned int numpts_in_pqr =
      ((is_p_point)? 1 : 0) +
      ((is_q_point)? 1 : 0) +
      ((is_r_point)? 1 : 0)  ;

    bool is_p_tsrc(false);
    bool has_p_endp_tsrc(false);
    if ( is_p_point ) {
      if ( same_points(p_, t.source_site()) ) {
        is_p_tsrc = true;
      }
    } else { // p is segment
      if (same_points(p_.source_site(), t.source_site()) or
          same_points(p_.target_site(), t.source_site())   ) {
        has_p_endp_tsrc = true;
      }
    }

    bool is_q_tsrc(false);
    bool has_q_endp_tsrc(false);
    if ( is_q_point ) {
      if ( same_points(q_, t.source_site()) ) {
        is_q_tsrc = true;
      }
    } else { // q is segment
      if (same_points(q_.source_site(), t.source_site()) or
          same_points(q_.target_site(), t.source_site())   ) {
        has_q_endp_tsrc = true;
      }
    }

    bool is_r_tsrc(false);
    bool has_r_endp_tsrc(false);
    if ( is_r_point ) {
      if ( same_points(r_, t.source_site()) ) {
        is_r_tsrc = true;
      }
    } else { // r is segment
      if (same_points(r_.source_site(), t.source_site()) or
          same_points(r_.target_site(), t.source_site())   ) {
        has_r_endp_tsrc = true;
      }
    }

    unsigned int num_common_endp_tsrc =
      ((has_p_endp_tsrc)? 1 : 0) +
      ((has_q_endp_tsrc)? 1 : 0) +
      ((has_r_endp_tsrc)? 1 : 0)  ;

    CGAL_SDG_DEBUG(std::cout << "debug num_common_endp_tsrc="
      << num_common_endp_tsrc << std::endl;);

    unsigned int numendpts_of_t = 0;

    Sign d1, d2;
    if ( is_p_tsrc or is_q_tsrc or is_r_tsrc ) {
      d1 = ZERO;
      ++numendpts_of_t;
    } else if ( num_common_endp_tsrc >= 2 ) {
      d1 = ZERO;
    } else {
      d1 = incircle_p(t.source_site());
    }
    if ( d1 == NEGATIVE ) { return NEGATIVE; }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy d1=" << d1 << std::endl;);

    bool is_p_ttrg(false);
    bool has_p_endp_ttrg(false);
    if ( is_p_point ) {
      if ( same_points(p_, t.target_site()) ) {
        is_p_ttrg = true;
      }
    } else { // p is segment
      if (same_points(p_.source_site(), t.target_site()) or
          same_points(p_.target_site(), t.target_site())   ) {
        has_p_endp_ttrg = true;
      }
    }

    bool is_q_ttrg(false);
    bool has_q_endp_ttrg(false);
    if ( is_q_point ) {
      if ( same_points(q_, t.target_site()) ) {
        is_q_ttrg = true;
      }
    } else { // q is segment
      if (same_points(q_.source_site(), t.target_site()) or
          same_points(q_.target_site(), t.target_site())   ) {
        has_q_endp_ttrg = true;
      }
    }

    bool is_r_ttrg(false);
    bool has_r_endp_ttrg(false);
    if ( is_r_point ) {
      if ( same_points(r_, t.target_site()) ) {
        is_r_ttrg = true;
      }
    } else { // r is segment
      if (same_points(r_.source_site(), t.target_site()) or
          same_points(r_.target_site(), t.target_site())   ) {
        has_r_endp_ttrg = true;
      }
    }

    unsigned int num_common_endp_ttrg =
      ((has_p_endp_ttrg)? 1 : 0) +
      ((has_q_endp_ttrg)? 1 : 0) +
      ((has_r_endp_ttrg)? 1 : 0)  ;

    CGAL_SDG_DEBUG(std::cout << "debug num_common_endp_ttrg="
      << num_common_endp_ttrg << std::endl;);

    if ( is_p_ttrg or is_q_ttrg or is_r_ttrg ) {
      d2 = ZERO;
      ++numendpts_of_t;
    } else if ( num_common_endp_ttrg >= 2 ) {
      d2 = ZERO;
    } else {
      d2 = incircle_p(t.target_site());
    }
    if ( d2 == NEGATIVE ) { return NEGATIVE; }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy d2=" << d2 << std::endl;);

    CGAL_assertion(numendpts_of_t < 2);

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy numendpts_of_t= "
      << numendpts_of_t << std::endl;);

   // Line_2 l = compute_supporting_line(t.supporting_site());
   // Sign sl = incircle(l, type);
  
    Segment_2 tseg = t.segment();
  
    RT dup;
    if (p_.is_point()) {
      Point_2 pp = p_.point();
      dup = CGAL::max(CGAL::abs(ux_ - pp.x() * uz_),
                      CGAL::abs(uy_ - pp.y() * uz_));
    } else {//p_ is segment
      Segment_2 pseg = p_.segment();
      if (pseg.is_horizontal()) {
        dup = CGAL::abs(uy_ - pseg.source().y() * uz_);
      } else {//p_ is vertical
        dup = CGAL::abs(ux_ - pseg.source().x() * uz_);
      }
    }
  
    if (tseg.is_horizontal()) {
      Point_2 txmax = ( CGAL::compare(tseg.source().x(), tseg.target().x())
                       == SMALLER ) ? tseg.target()
                                    : tseg.source();
      Point_2 txmin = ( CGAL::compare(tseg.source().x(), tseg.target().x())
                       == SMALLER ) ? tseg.source()
                                    : tseg.target();
      if (CGAL::compare(uy_ - dup , uz_ * tseg.source().y()) == SMALLER &&
          CGAL::compare(uz_ * tseg.source().y() , uy_ + dup) == SMALLER &&
          CGAL::compare(ux_ - dup , uz_ * txmax.x()) == SMALLER ) {
        return NEGATIVE;
      } else if (CGAL::compare(uy_ - dup , uz_ * tseg.source().y()) == LARGER or
                 CGAL::compare(uz_ * tseg.source().y() , uy_ + dup) == LARGER or
                 CGAL::compare(ux_ - dup , uz_ * txmax.x()) == LARGER or
                 CGAL::compare(ux_ + dup , uz_ * txmin.x()) == SMALLER) {
        return POSITIVE;
      } else {
        return ZERO;
      }
                     
    } else {// tseg is vertical
      Point_2 tymax = ( CGAL::compare(tseg.source().y(), tseg.target().y())
                       == SMALLER ) ? tseg.target()
                                    : tseg.source();
      Point_2 tymin = ( CGAL::compare(tseg.source().y(), tseg.target().y())
                       == SMALLER ) ? tseg.source()
                                    : tseg.target();
      if (CGAL::compare(ux_ - dup , uz_ * tseg.source().x()) == SMALLER &&
          CGAL::compare(uz_ * tseg.source().x() , ux_ + dup) == SMALLER &&
          CGAL::compare(uy_ - dup , uz_ * tymax.y()) == SMALLER ) {
        return NEGATIVE;
      } else if (CGAL::compare(ux_ - dup , uz_ * tseg.source().x()) == LARGER or
                 CGAL::compare(uz_ * tseg.source().x() , ux_ + dup) == LARGER or
                 CGAL::compare(uy_ - dup , uz_ * tymax.y()) == LARGER or
                 CGAL::compare(uy_ + dup , uz_ * tymin.y()) == SMALLER) {
        return POSITIVE;
      } else {
        return ZERO;
      }
      
    }//end of tseg vertical
    
  }

  inline
  bool
  compute_helper(const Site_2& p, const Site_2& q, const Site_2& r,
      const Site_2& t, 
      Site_2& sqpnt, Site_2& other_of_t, Site_2& other_of_seg)
  const
  {
    CGAL_assertion(t.is_segment());

    bool is_p_point = p.is_point();
    bool is_q_point = q.is_point();
    bool is_r_point = r.is_point();

    unsigned int numpts = 
      ((is_p_point)? 1 : 0) + 
      ((is_q_point)? 1 : 0) +
      ((is_r_point)? 1 : 0)  ;

    CGAL_SDG_DEBUG(std::cout << "debug compute_helper #pts=" << numpts << std::endl;);

    if (numpts == 3) {
      return false;
    }

    // here and on, there are 1 or 2 points in {p,q,r}


    bool is_p_tsrc(false);
    bool is_p_ttrg(false);
    bool is_p_endp_of_t(false);

    if (is_p_point) {
      is_p_tsrc = same_points(p, t.source_site());
      is_p_ttrg = same_points(p, t.target_site());
      is_p_endp_of_t = is_p_tsrc or is_p_ttrg;

      if (is_p_endp_of_t) {
        sqpnt = p;
      }
    }

    bool is_q_tsrc(false);
    bool is_q_ttrg(false);
    bool is_q_endp_of_t(false);

    if (is_q_point) {
      is_q_tsrc = same_points(q, t.source_site());
      is_q_ttrg = same_points(q, t.target_site());
      is_q_endp_of_t = is_q_tsrc or is_q_ttrg;
      if (is_q_endp_of_t) {
        sqpnt = q;
      }
    }

    bool is_r_tsrc(false);
    bool is_r_ttrg(false);
    bool is_r_endp_of_t(false);

    if (is_r_point) {
      is_r_tsrc = same_points(r, t.source_site());
      is_r_ttrg = same_points(r, t.target_site());
      is_r_endp_of_t = is_r_tsrc or is_r_ttrg;
      if (is_r_endp_of_t) {
        sqpnt = r;
      }
    }

    unsigned int numendpts_of_t = 
      ((is_p_endp_of_t)? 1 : 0) + 
      ((is_q_endp_of_t)? 1 : 0) +
      ((is_r_endp_of_t)? 1 : 0)  ; 

    CGAL_SDG_DEBUG(std::cout << "debug compute_helper #endpts_of_t=" << 
      numendpts_of_t << std::endl;);

    if (numendpts_of_t == 0) {

      bool is_psrc_tsrc(false),
           is_ptrg_tsrc(false),
           is_psrc_ttrg(false),
           is_ptrg_ttrg(false),
           have_common_p_tsrc(false),
           have_common_p_ttrg(false),
           have_common_p_t(false);

      if (not is_p_point) {
        CGAL_assertion( not same_segments(p, t) );
        is_psrc_tsrc = same_points(p.source_site(), t.source_site());
        is_ptrg_tsrc = same_points(p.target_site(), t.source_site());
        is_psrc_ttrg = same_points(p.source_site(), t.target_site());
        is_ptrg_ttrg = same_points(p.target_site(), t.target_site());
        have_common_p_tsrc = is_psrc_tsrc or is_ptrg_tsrc;
        have_common_p_ttrg = is_psrc_ttrg or is_ptrg_ttrg;
        have_common_p_t = have_common_p_tsrc or have_common_p_ttrg;
      }

      bool is_qsrc_tsrc(false),
           is_qtrg_tsrc(false),
           is_qsrc_ttrg(false),
           is_qtrg_ttrg(false),
           have_common_q_tsrc(false),
           have_common_q_ttrg(false),
           have_common_q_t(false);

      if (not is_q_point) {
        CGAL_assertion( not same_segments(q, t) );
        is_qsrc_tsrc = same_points(q.source_site(), t.source_site());
        is_qtrg_tsrc = same_points(q.target_site(), t.source_site());
        is_qsrc_ttrg = same_points(q.source_site(), t.target_site());
        is_qtrg_ttrg = same_points(q.target_site(), t.target_site());
        have_common_q_tsrc = is_qsrc_tsrc or is_qtrg_tsrc;
        have_common_q_ttrg = is_qsrc_ttrg or is_qtrg_ttrg;
        have_common_q_t = have_common_q_tsrc or have_common_q_ttrg;
      }

      bool is_rsrc_tsrc(false),
           is_rtrg_tsrc(false),
           is_rsrc_ttrg(false),
           is_rtrg_ttrg(false),
           have_common_r_tsrc(false),
           have_common_r_ttrg(false),
           have_common_r_t(false);

      if (not is_r_point) {
        CGAL_assertion( not same_segments(r, t) );
        is_rsrc_tsrc = same_points(r.source_site(), t.source_site());
        is_rtrg_tsrc = same_points(r.target_site(), t.source_site());
        is_rsrc_ttrg = same_points(r.source_site(), t.target_site());
        is_rtrg_ttrg = same_points(r.target_site(), t.target_site());
        have_common_r_tsrc = is_rsrc_tsrc or is_rtrg_tsrc;
        have_common_r_ttrg = is_rsrc_ttrg or is_rtrg_ttrg;
        have_common_r_t = have_common_r_tsrc or have_common_r_ttrg;
      }

      unsigned int numcommon = 
      ((have_common_p_t)? 1 : 0) + 
      ((have_common_q_t)? 1 : 0) +
      ((have_common_r_t)? 1 : 0)  ; 
 
      CGAL_SDG_DEBUG(std::cout << "debug compute_helper #numcommon=" << 
        numcommon << std::endl;);

      CGAL_assertion(numcommon < 3);

      if (numcommon < 2) {
        return false;
      }

      // here, numcommon == 2

      unsigned int numcommon_tsrc = 
      ((have_common_p_tsrc)? 1 : 0) + 
      ((have_common_q_tsrc)? 1 : 0) +
      ((have_common_r_tsrc)? 1 : 0)  ; 

      unsigned int numcommon_ttrg = 
      ((have_common_p_ttrg)? 1 : 0) + 
      ((have_common_q_ttrg)? 1 : 0) +
      ((have_common_r_ttrg)? 1 : 0)  ; 

      CGAL_assertion( numcommon_tsrc + numcommon_ttrg == 2 );

      if ( numcommon_tsrc == numcommon_ttrg )  { // both equal 1
        return false;
      }

      // here either numcommon_tsrc==2 or numcommon_ttrg==2

      if (numcommon_tsrc > 0) {
        // here, numcommon_tsrc == 2
        sqpnt = t.source_site();
        other_of_t = t.target_site();
      } else {
        // here, numcommon_ttrg == 2
        sqpnt = t.target_site();
        other_of_t = t.source_site();
      }

      if (have_common_p_t and have_common_q_t) {
        compute_helper_two_seg(p, q, sqpnt, other_of_seg); 
      } else if (have_common_q_t and have_common_r_t) {
        compute_helper_two_seg(q, r, sqpnt, other_of_seg); 
      } else if (have_common_r_t and have_common_p_t) {
        compute_helper_two_seg(r, p, sqpnt, other_of_seg); 
      } else {
        CGAL_assertion(false);
      }

      return true;
    }

    // philaris: tocheck 
    CGAL_assertion( numendpts_of_t == 1 );

    if (is_p_tsrc or is_q_tsrc or is_r_tsrc) {
      other_of_t = t.target_site();
    } else {
      other_of_t = t.source_site();
    }

    if (is_p_endp_of_t) {
      if (q.is_segment()) {
        bool is_p_qsrc = same_points(p, q.source_site());
        bool is_p_qtrg = same_points(p, q.target_site());
        if (is_p_qsrc or is_p_qtrg) {
          other_of_seg = is_p_qsrc ? q.target_site() : q.source_site();
          return true;
        }
      }
      if (r.is_segment()) {
        bool is_p_rsrc = same_points(p, r.source_site());
        bool is_p_rtrg = same_points(p, r.target_site());
        if (is_p_rsrc or is_p_rtrg) {
          other_of_seg = is_p_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

    } // end of case: is_p_endp_of_t

    if (is_q_endp_of_t) {
      if (r.is_segment()) {
        bool is_q_rsrc = same_points(q, r.source_site());
        bool is_q_rtrg = same_points(q, r.target_site());
        if (is_q_rsrc or is_q_rtrg) {
          other_of_seg = is_q_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

      if (p.is_segment()) {
        bool is_q_psrc = same_points(q, p.source_site());
        bool is_q_ptrg = same_points(q, p.target_site());
        if (is_q_psrc or is_q_ptrg) {
          other_of_seg = is_q_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

    } // end of case: is_q_endp_of_t

    if (is_r_endp_of_t) {
      if (p.is_segment()) {
        bool is_r_psrc = same_points(r, p.source_site());
        bool is_r_ptrg = same_points(r, p.target_site());
        if (is_r_psrc or is_r_ptrg) {
          other_of_seg = is_r_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

      if (q.is_segment()) {
        bool is_r_qsrc = same_points(r, q.source_site());
        bool is_r_qtrg = same_points(r, q.target_site());
        if (is_r_qsrc or is_r_qtrg) {
          other_of_seg = is_r_qsrc ? q.target_site() : q.source_site();
          return true;
        }
      }

    } // end of case: is_r_endp_of_t

    return false;

  }

  inline
  void
  compute_helper_two_seg(
      const Site_2& a, const Site_2& b, 
      const Site_2& common_site, Site_2& other_of_seg)
  const
  {
    CGAL_assertion(a.is_segment());
    CGAL_assertion(b.is_segment());

    CGAL_SDG_DEBUG(std::cout << "debug compute_helper_two_seg entering with "
      << a << " and " << b << " having common " 
      << common_site << std::endl;);

    if (a.segment().is_horizontal() or a.segment().is_vertical()) {
      if ( same_points(common_site, b.source_site()) ) {
        other_of_seg = b.target_site();
      } else {
        other_of_seg = b.source_site();
      }
    } else {
      CGAL_assertion( 
          b.segment().is_horizontal() or b.segment().is_vertical() );

      if ( same_points(common_site, a.source_site()) ) {
        other_of_seg = a.target_site();
      } else {
        other_of_seg = a.source_site();
      }

    }
  } // end of compute_helper_two_seg


  //--------------------------------------------------------------------------

  

  Sign incircle_s(const Site_2& t) const 
  {
    CGAL_precondition( t.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s (pqrt) = "
      << "(" << p_ << ") (" << q_ << ") (" << r_ << ") "
      << "(" << t << ")" << std::endl;);

    if ( is_degenerate_Voronoi_circle() ) {
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle
      if (  !same_points( p_ref(), t.source_site() ) &&
	    !same_points( p_ref(), t.target_site() )  ) {
	return POSITIVE;
      }

      CGAL_assertion( v_type == PSS );

      if ( p_.is_segment() &&
	   same_segments(p_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      if ( q_.is_segment() &&
	   same_segments(q_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      if ( r_.is_segment() &&
	   same_segments(r_.supporting_site(),
			 t.supporting_site()) ) {
	return ZERO;
      }

      Site_2 pr;
      Site_2 sp, sq;
      if ( p_.is_point() ) {
	CGAL_assertion( q_.is_segment() && r_.is_segment() );
	pr = p_;
	sp = q_;
	sq = r_;
      } else if ( q_.is_point() ) {
	CGAL_assertion( r_.is_segment() && p_.is_segment() );
	pr = q_;
	sp = r_;
	sq = p_;
      } else {
	CGAL_assertion( p_.is_segment() && q_.is_segment() );
	pr = r_;
	sp = p_;
	sq = q_;
      }

      Point_2 pq = sq.source(), pp = sp.source(), pt = t.source();

      if ( same_points(sp.source_site(), pr) ) { pp = sp.target(); }
      if ( same_points(sq.source_site(), pr) ) { pq = sq.target(); }
      if ( same_points( t.source_site(), pr) ) { pt =  t.target(); }

      Point_2 pr_ = pr.point();

      if ( CGAL::orientation(pr_, pp, pt) == LEFT_TURN &&
	   CGAL::orientation(pr_, pq, pt) == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return ZERO;
    } // if ( is_degenerate_Voronoi_circle() )

    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_s(t, PPP_Type());
      break;
    case PPS:
      s = incircle_s(t, PPS_Type());
      break;
    case PSS:
      s = incircle_s(t, PSS_Type());
      break;
    case SSS:
      s = incircle_s(t, SSS_Type());
      break;
    }

    return s;
  }

  Sign incircle_s_no_easy(const Site_2& t) const
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_s_no_easy(t, PPP_Type());
      break;
    case PPS:
      s = incircle_s_no_easy(t, PPS_Type());
      break;
    case PSS:
      s = incircle_s_no_easy(t, PSS_Type());
      break;
    case SSS:
      s = incircle_s_no_easy(t, SSS_Type());
      break;
    }

    return s;
  }

  //--------------------------------------------------------------------------
  //  subpredicates for the incircle test
  //--------------------------------------------------------------------------


public:
  bool is_degenerate_Voronoi_circle() const
  {
    if ( v_type != PSS ) { return false; }

    if ( p_.is_point() ) {
      return ( is_endpoint_of(p_, q_) && is_endpoint_of(p_, r_) );
    } else if ( q_.is_point() ) {
      return ( is_endpoint_of(q_, p_) && is_endpoint_of(q_, r_) );
    } else {
      CGAL_assertion( r_.is_point() );
      return ( is_endpoint_of(r_, p_) && is_endpoint_of(r_, q_) );
    }
  }


  //--------------------------------------------------------------------------

private:

  //--------------------------------------------------------------------------
  //  the reference point (valid if v_type != SSS)
  //--------------------------------------------------------------------------

  Site_2 p_ref() const
  {
    CGAL_precondition ( v_type != SSS );


    if ( v_type == PPS ) {
      //CGAL_SDG_DEBUG(std::cout << "debug p_ref pps_idx=" << pps_idx << std::endl;);

      if ( pps_idx == 0 ) { 
        CGAL_assertion( p_.is_point());
        return p_; 
      }

      if ( pps_idx == 1 ) { 
        CGAL_assertion( q_.is_point());
        return q_; 
      }

      //CGAL_SDG_DEBUG(std::cout << "debug p_ref about to return r=" << r_ << std::endl;);

      CGAL_assertion( r_.is_point());
      return r_;
    }

    if ( p_.is_point() ) {
      return p_;
    } else if ( q_.is_point() ) {
      return q_;
    } else {
      CGAL_assertion( r_.is_point() );
      return r_;
    }
  }


public:
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //                           access methods
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  inline FT x(Integral_domain_without_division_tag) const {
    return CGAL::to_double(hx()) / CGAL::to_double(hw());
  }
  inline FT y(Integral_domain_without_division_tag) const {
    return CGAL::to_double(hy()) / CGAL::to_double(hw());
  }

  inline FT x(Field_tag) const { return hx() / hw(); }
  inline FT y(Field_tag) const { return hy() / hw(); }

  inline FT x() const {
      typedef Algebraic_structure_traits<FT> AST;
      return x(typename AST::Algebraic_category());
  }

  inline FT y() const {
      typedef Algebraic_structure_traits<FT> AST;
      return y(typename AST::Algebraic_category());
  }

  FT hx() const {
    // philaris: changed to one type
    //if ( v_type == PPP ) { return ux_ppp; }
    //if ( v_type == PPS ) { return to_ft(ux_pps); }
    //return to_ft(ux);
    return ux_;
  }

  FT hy() const {
    // philaris: changed to one type
    //if ( v_type == PPP ) { return uy_ppp; }
    //if ( v_type == PPS ) { return to_ft(uy_pps); }
    //return to_ft(uy);
    return uy_;
  }

  FT hw() const {
    // philaris: changed to one type
    //if ( v_type == PPP ) { return uz_ppp; }
    //if ( v_type == PPS ) { return to_ft(uz_pps); }
    //return to_ft(uz);
    return uz_;
  }

  FT radius() const {
    switch (v_type) {
    case PPP:    case PPS:    case PSS:
      {
	Point_2 pref = p_ref().point();
	//FT absdx = CGAL::abs(x() - pref.x());
	//FT absdy = CGAL::abs(y() - pref.y());
        return CGAL::max( CGAL::abs(x() - pref.x()), 
		          CGAL::abs(y() - pref.y()) );
      }
      break;
    case SSS:
      {
	Line_2 l = compute_supporting_line(p_.supporting_site());
	Homogeneous_point_2 q = compute_linf_projection_hom(l, point());

	FT dx = CGAL::abs(x() - q.x());
	FT dy = CGAL::abs(y() - q.y());
	return CGAL::max(dx, dy);
      }
      break;
    default:
      return FT(0);
    }
  }

  Point_2 point() const {
    if ( is_degenerate_Voronoi_circle() ) {
      return degenerate_point();
    }

    return Point_2(x(), y());
  }


  Point_2 degenerate_point() const
  {
    CGAL_precondition( is_degenerate_Voronoi_circle() );
    return p_ref().point();
  }

  // philaris: the circle is in fact an Iso_rectangle_2
  typename K::Iso_rectangle_2 circle() const
  {
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    Point_2 pleftbot (point().x()-radius(), point().y()-radius());
    Point_2 prghttop (point().x()+radius(), point().y()+radius());
    return Iso_rectangle_2(pleftbot, prghttop);
  }

  vertex_t type() const { return v_type; }

public:
  Voronoi_vertex_ring_C2(const Site_2& p,
			 const Site_2& q,
			 const Site_2& r)
    : p_(p), q_(q), r_(r)
  {
    compute_vertex(p, q, r);
  }

  //--------------------------------------------------------------------------

  Sign incircle(const Site_2& t) const 
  {
    Sign s;

    CGAL_SDG_DEBUG(std::cout << "debug ring incircle t=" << t << std::endl;);

    if ( t.is_point() ) {
      s = incircle_p(t);
    } else {
      CGAL_SDG_DEBUG(std::cout << "debug about to run incircle_s with t=" 
        << t << std::endl;);
      s = incircle_s(t);
    }

    return s;
  }

  Sign incircle_no_easy(const Site_2& t) const 
  {
    Sign s;

    if ( t.is_point() ) {
      s = incircle_p_no_easy(t);
    } else {
      s = incircle_s_no_easy(t);
    }

    return s;
  }

  //--------------------------------------------------------------------------


  Orientation orientation(const Line_2& l) const 
  {
    Orientation o(COLLINEAR);
    switch ( v_type ) {
    case PPP:
      o = orientation(l, PPP_Type());
      break;
    case PPS:
      o = orientation(l, PPS_Type());
      break;
    case PSS:
      o = orientation(l, PSS_Type());
      break;
    case SSS:
      o = orientation(l, SSS_Type());
      break;
    }

    return o;
  }

  Oriented_side oriented_side(const Line_2& l) const
  {
    Orientation o = orientation(l);

    if ( o == COLLINEAR ) { return ON_ORIENTED_BOUNDARY; }
    return ( o == LEFT_TURN ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }


  // L_inf refinement

  inline
  Comparison_result
  linf_refinement( Homogeneous_point_2& lrefhp ) const
  {
    Point_2 vv ( ux_, uy_, uz_ );

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    FT difxvl = vv.x() - lrefhp.x();
    FT difyvl = vv.y() - lrefhp.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);

    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    if (p_.is_point()) {
      Point_2 pp = p_.point();
      FT difxvp = vv.x() - pp.x();
      FT difyvp = vv.y() - pp.y();
      FT absdifxvp = CGAL::abs(difxvp);
      FT absdifyvp = CGAL::abs(difyvp);
      Comparison_result cmppabsxy = CGAL::compare(absdifxvp, absdifyvp);
      if (not ( (cmplabsxy == SMALLER) and (cmppabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvp) == EQUAL) {
          compare_p = CGAL::compare(absdifyvl, absdifyvp);
        }
      }
      if (not ( (cmplabsxy == LARGER ) and (cmppabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvp) == EQUAL) {
          CGAL_assertion(compare_p == EQUAL);
          compare_p = CGAL::compare(absdifxvl, absdifxvp);
        }
      }
    }

    if (q_.is_point()) {
      Point_2 qq = q_.point();
      FT difxvq = vv.x() - qq.x();
      FT difyvq = vv.y() - qq.y();
      FT absdifxvq = CGAL::abs(difxvq);
      FT absdifyvq = CGAL::abs(difyvq);
      Comparison_result cmpqabsxy = CGAL::compare(absdifxvq, absdifyvq);
      if (not ( (cmplabsxy == SMALLER) and (cmpqabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvq) == EQUAL) {
          compare_q = CGAL::compare(absdifyvl, absdifyvq);
        }
      }
      if (not ( (cmplabsxy == LARGER ) and (cmpqabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvq) == EQUAL) {
          CGAL_assertion(compare_q == EQUAL);
          compare_q = CGAL::compare(absdifxvl, absdifxvq);
        }
      }
    }

    if (r_.is_point()) {
      Point_2 rr = r_.point();
      FT difxvr = vv.x() - rr.x();
      FT difyvr = vv.y() - rr.y();
      FT absdifxvr = CGAL::abs(difxvr);
      FT absdifyvr = CGAL::abs(difyvr);
      Comparison_result cmprabsxy = CGAL::compare(absdifxvr, absdifyvr);
      if (not ( (cmplabsxy == SMALLER) and (cmprabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvr) == EQUAL) {
          compare_r = CGAL::compare(absdifyvl, absdifyvr);
        }
      }
      if (not ( (cmplabsxy == LARGER ) and (cmprabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvr) == EQUAL) {
          CGAL_assertion(compare_r == EQUAL);
          compare_r = CGAL::compare(absdifxvl, absdifxvr);
        }
      }
    }

    CGAL_SDG_DEBUG(std::cout << "debug compare p q r = " 
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) or 
        (compare_q == SMALLER) or
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    /*
    if ((compare_p == LARGER) or 
        (compare_q == LARGER) or
        (compare_r == LARGER)   ) {
      // tocheck
      return LARGER;
    }
    */
    return EQUAL;

  }

  

  //--------------------------------------------------------------------------

private:
  const Site_2& p_, q_, r_;

  vertex_t v_type;

  // index that indicates the refence point for the case PPS
  short pps_idx;

  // philaris: different types are not needed any more
  // the case ppp
  //RT ux_ppp, uy_ppp, uz_ppp;

  // the case pps
  //Sqrt_1 ux_pps, uy_pps, uz_pps;

  // the case pss and sss
  //Sqrt_3 ux, uy, uz;

  RT ux_, uy_, uz_;
};


} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H
