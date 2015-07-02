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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H


#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_x_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_y_2.h>
#include <CGAL/Side_of_bounded_square_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
#include <CGAL/tuple.h>


namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

#ifndef CGAL_CFG_NO_CPP0X_TUPLE
#if (defined(_MSC_VER) && (_MSC_VER < 1700))
#define sdg_tuple_maker cpp11::make_tuple
#else
#define sdg_tuple_maker std::forward_as_tuple
#endif
#else
#define sdg_tuple_maker cpp11::make_tuple
#endif

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
  typedef typename Base::Direction_2         Direction_2;
  typedef typename Base::Site_2              Site_2;
  typedef typename Base::FT                  FT;
  typedef typename Base::RT                  RT;

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Bounded_side        Bounded_side;
  typedef typename Base::Sign                Sign;

  typedef typename Base::Polychainline_2     Polychainline_2;

  typedef typename Base::Bearing Bearing;

  using Base::compute_supporting_line;
  using Base::compute_linf_projection_hom;
  using Base::compute_linf_projection_nonhom;
  using Base::oriented_side_of_line;
  using Base::opposite_line;
  using Base::compute_line_from_to;
  using Base::compute_horizontal_projection;
  using Base::compute_vertical_projection;
  using Base::compute_linf_perpendicular;
  using Base::is_site_horizontal;
  using Base::is_site_vertical;
  using Base::is_site_h_or_v;
  using Base::is_line_h_or_v;
  using Base::test_star;
  using Base::compute_neg_45_line_at;
  using Base::compute_pos_45_line_at;
  using Base::compute_hor_line_at;
  using Base::compute_ver_line_at;
  using Base::has_positive_slope;
  using Base::are_in_same_open_halfspace_of;
  using Base::horseg_y_coord;
  using Base::verseg_x_coord;
  using Base::hvseg_coord;
  using Base::compute_intersection_of_lines;
  using Base::is_orth_dist_smaller_than_pt_dist;
  using Base::touch_same_side;
  using Base::coord_at;
  using Base::orient_lines_linf;
  using Base::compute_line_dir;
  using Base::bisector_linf_line;
  using Base::is_endpoint_of;
  using Base::orient_line_endp;
  using Base::orient_line_nonendp;
  using Base::bearing;
  using Base::bearing_diff;
  using Base::center_from_corner_and_pt;
  using Base::points_inside_touching_sides_v;
  using Base::center_from_opposite_corners;
  using Base::center_from_same_side_corners;
  using Base::is_on_hv_seg_line;

private:
  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>
          Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K>
          Are_same_segments_2;
  typedef Side_of_bounded_square_2<K>    Side_of_bounded_square_2_Type;
  typedef Bisector_Linf<K>               Bisector_Linf_Type;

  typedef SegmentDelaunayGraph_2::Compare_x_2<K> Compare_x_2_Sites_Type;
  typedef SegmentDelaunayGraph_2::Compare_y_2<K> Compare_y_2_Sites_Type;

  Are_same_points_2                same_points;
  Are_same_segments_2              same_segments;
  Side_of_bounded_square_2_Type    side_of_bounded_square;
  Bisector_Linf_Type               bisector_linf;
  Compare_x_2_Sites_Type           scmpx;
  Compare_y_2_Sites_Type           scmpy;

private:
  //--------------------------------------------------------------------------

  inline void
  compute_v_if_not_computed() const {
    if (! is_v_computed) {
      compute_vertex(p_, q_, r_);
      is_v_computed = true;
    }
  }

  void
  compute_ppp(const Site_2& sp, const Site_2& sq, const Site_2& sr)
  const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
		       sr.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug vring entering compute_ppp"
        << std::endl;);

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

    compute_ppp(p, q, r);
  }

  void
  compute_ppp(const Point_2& p, const Point_2& q, const Point_2& r)
  const
  {
    RT x_min, x_max, y_min, y_max;
    RT two_x_center, two_y_center;
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
      two_y_center = p.y() + q.y();
      is_set_y_center = true;

      CGAL_SDG_DEBUG(std::cout << "debug set " <<
        " py, qy =" << p.y() << ' ' << q.y() <<
        " two_y_center=" << two_y_center << std::endl;);

      Comparison_result cmpxrothers = CGAL::compare(r.x(), p.x());
      if (cmpxrothers == SMALLER) {
        CGAL_SDG_DEBUG(std::cout << "debug r is left of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  && (cmpyrq == LARGER)) ||
            ((cmpyrp == SMALLER) && (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two_y_center - r.y();
            is_set_y_min = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_min=" <<
              y_min << std::endl;);
          } else {
            y_max = two_y_center - r.y();
            is_set_y_max = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_max=" <<
              y_max << std::endl;);
          }
        }
      } else if (cmpxrothers == LARGER) {
        CGAL_SDG_DEBUG(std::cout << "debug r is right of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  && (cmpyrq == LARGER)) ||
            ((cmpyrp == SMALLER) && (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two_y_center - r.y();
            is_set_y_min = true;
            CGAL_SDG_DEBUG(std::cout << "debug set y_min=" <<
              y_min << std::endl;);
          } else {
            y_max = two_y_center - r.y();
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
      if (! is_set_y_min) {
        y_min = q.y();
      }
      if (! is_set_y_max) {
        y_max = p.y();
      }
    } else if (cmpyqp == LARGER) { // q.y() > p.y()
      if (! is_set_y_min) {
        y_min = p.y();
      }
      if (! is_set_y_max) {
        y_max = q.y();
      }
    } else { //  q.y() = p.y()
      if (! is_set_y_min) {
        y_min = p.y();
      }
      if (! is_set_y_max) {
        y_max = p.y();
      }
      two_x_center = p.x() + q.x();
      is_set_x_center = true;

      Comparison_result cmpyrothers = CGAL::compare(r.y(), p.y());
      if (cmpyrothers == SMALLER) {
        Comparison_result cmpxrp = CGAL::compare(r.x(), p.x());
        Comparison_result cmpxrq = CGAL::compare(r.x(), q.x());
        if (((cmpxrp == LARGER)  && (cmpxrq == LARGER)) ||
            ((cmpxrp == SMALLER) && (cmpxrq == SMALLER))
           ) {
          // do fix
          if (cmpxrp == LARGER) {
            x_min = two_x_center - r.x();
            is_set_x_min = true;
          } else {
            x_max = two_x_center - r.x();
            is_set_x_max = true;
          }
        }
      } else if (cmpyrothers == LARGER) {
        Comparison_result cmpxrp = CGAL::compare(r.x(), p.x());
        Comparison_result cmpxrq = CGAL::compare(r.x(), q.x());
        if (((cmpxrp == LARGER)  && (cmpxrq == LARGER)) ||
            ((cmpxrp == SMALLER) && (cmpxrq == SMALLER))
           ) {
          // do fix
          if (cmpxrp == LARGER) {
            x_min = two_x_center - r.x();
            is_set_x_min = true;
          } else {
            x_max = two_x_center - r.x();
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
        if (! is_set_x_min) {
          x_min = r.x();
        }
    } else if (cmpxrmin == LARGER) {
      // here r.x() > x_min
      if (cmpxrmax == LARGER) {
        // here x_min <= x_max < r.x()
        if (! is_set_x_max) {
          x_max = r.x();
        }
      } else if (cmpxrmax == SMALLER) {
        // x_min < r.x() < x_max
        // do nothing
      } else { // r.x() = x_max
        // r.x() = p.x() || r.x() = q.x()
        if (CGAL::compare(r.x(), p.x()) == EQUAL) {
          two_y_center = p.y() + r.y();
          //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
          Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
          if ((cmpyqp == LARGER) && (cmpyqr == LARGER)) {
            y_min = two_y_center - q.y();
            is_set_y_min = true;
          }
          if ((cmpyqp == SMALLER) && (cmpyqr == SMALLER)) {
            y_max = two_y_center - q.y();
            is_set_y_max = true;
          }
        } else {
          two_y_center = q.y() + r.y();
          Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
          Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
          if ((cmpypq == LARGER) && (cmpypr == LARGER)) {
            y_min = two_y_center - p.y();
            is_set_y_min = true;
          }
          if ((cmpypq == SMALLER) && (cmpypr == SMALLER)) {
            y_max = two_y_center - p.y();
            is_set_y_max = true;
          }
        }
        is_set_y_center = true;
      }
    } else {
      // here r.x() = x_min
      // r.x() = p.x() || r.x() = q.x()
      if (CGAL::compare(r.x(), p.x()) == EQUAL) {
        CGAL_SDG_DEBUG(std::cout << "debug r.x = p.x" << std::endl;);
        // r.x() = p.x()
        two_y_center = p.y() + r.y();
        //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
        Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
        if ((cmpyqp == LARGER) && (cmpyqr == LARGER)) {
          CGAL_SDG_DEBUG(std::cout << "debug q is above p, r" << std::endl;);
          y_min = two_y_center - q.y();
          is_set_y_min = true;
        }
        if ((cmpyqp == SMALLER) && (cmpyqr == SMALLER)) {
          CGAL_SDG_DEBUG(std::cout << "debug q is below p, r" << std::endl;);
          y_max = two_y_center - q.y();
          is_set_y_max = true;
        }
      } else {
        // r.x() = q.x()
        CGAL_SDG_DEBUG(std::cout << "debug r.x = q.x" << std::endl;);
        two_y_center = q.y() + r.y();
        Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
        Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
        if ((cmpypq == LARGER) && (cmpypr == LARGER)) {
          CGAL_SDG_DEBUG(std::cout << "debug p is above q, r" << std::endl;);
          y_min = two_y_center - p.y();
          is_set_y_min = true;
        }
        if ((cmpypq == SMALLER) && (cmpypr == SMALLER)) {
          CGAL_SDG_DEBUG(std::cout << "debug p is below q, r" << std::endl;);
          y_max = two_y_center - p.y();
          is_set_y_max = true;
        }
      }
      is_set_y_center = true;
    }

    Comparison_result cmpyrmin = CGAL::compare(r.y(), y_min);
    Comparison_result cmpyrmax = CGAL::compare(r.y(), y_max);
    if (cmpyrmin == SMALLER) {
      // here r.y() < y_min <= y_max
      if (! is_set_y_min) {
        y_min = r.y();
      }
    } else if (cmpyrmin == LARGER) {
      // here r.y() > y_min
      if (cmpyrmax == LARGER) {
        // here y_min <= y_max < r.y()
        if (! is_set_y_max) {
          y_max = r.y();
        }
      } else if (cmpyrmax == SMALLER) {
        // y_min < r.y() < y_max
        // do nothing
      } else { // r.y() = y_max
        // r.y() = p.y() || r.y() = q.y()
        if (CGAL::compare(r.y(), p.y()) == EQUAL) {
          two_x_center = p.x() + r.x();
          //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
          Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
          if ((cmpxqp == LARGER) && (cmpxqr == LARGER)) {
            x_min = two_x_center - q.x();
            is_set_x_min = true;
          }
          if ((cmpxqp == SMALLER) && (cmpxqr == SMALLER)) {
            x_max = two_x_center - q.x();
            is_set_x_max = true;
          }
        } else {
          two_x_center = q.x() + r.x();
          Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
          Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
          if ((cmpxpq == LARGER) && (cmpxpr == LARGER)) {
            x_min = two_x_center - p.x();
            is_set_x_min = true;
          }
          if ((cmpxpq == SMALLER) && (cmpxpr == SMALLER)) {
            x_max = two_x_center - p.x();
            is_set_x_max = true;
          }
        }
        is_set_x_center = true;
      }
    } else {
      // here r.y() = y_min
      // r.y() = p.y() || r.y() = q.y()
      if (CGAL::compare(r.y(), p.y()) == EQUAL) {
        two_x_center = p.x() + r.x();
        //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
        Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
        if ((cmpxqp == LARGER) && (cmpxqr == LARGER)) {
          x_min = two_x_center - q.x();
          is_set_x_min = true;
        }
        if ((cmpxqp == SMALLER) && (cmpxqr == SMALLER)) {
          x_max = two_x_center - q.x();
          is_set_x_max = true;
        }
      } else {
        two_x_center = q.x() + r.x();
        Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
        Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
        if ((cmpxpq == LARGER) && (cmpxpr == LARGER)) {
          x_min = two_x_center - p.x();
          is_set_x_min = true;
        }
        if ((cmpxpq == SMALLER) && (cmpxpr == SMALLER)) {
          x_max = two_x_center - p.x();
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
        CGAL_SDG_DEBUG(std::cout
            << "debug vring rectangle has to be made fatter" << std::endl;);
        // make rectangle fatter
        if (is_set_x_center) {
          CGAL_SDG_DEBUG(std::cout
              << "debug vring x_center already set" << std::endl;);
          // grow in both sides
          break;
        }
        // grow only if any point is inside vertical sides
	if (((CGAL::compare(p.x(), x_min) == EQUAL)   &&
	     (CGAL::compare(p.y(), y_max) == SMALLER) &&
	     (CGAL::compare(p.y(), y_min) == LARGER)     ) ||
	    ((CGAL::compare(q.x(), x_min) == EQUAL)   &&
	     (CGAL::compare(q.y(), y_max) == SMALLER) &&
	     (CGAL::compare(q.y(), y_min) == LARGER)     ) ||
	    ((CGAL::compare(r.x(), x_min) == EQUAL)   &&
	     (CGAL::compare(r.y(), y_max) == SMALLER) &&
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
        CGAL_SDG_DEBUG(std::cout
            << "debug vring rectangle has to be made taller" << std::endl;);
        // make rectangle taller
        if (is_set_y_center) {
          // grow in both sides
          CGAL_SDG_DEBUG(std::cout << "y_center already set" << std::endl;);
          break;
        }
        // grow only if any point is inside horizontal sides
	if (((CGAL::compare(p.y(), y_min) == EQUAL)   &&
	     (CGAL::compare(p.x(), x_max) == SMALLER) &&
	     (CGAL::compare(p.x(), x_min) == LARGER)     ) ||
	    ((CGAL::compare(q.y(), y_min) == EQUAL)   &&
	     (CGAL::compare(q.x(), x_max) == SMALLER) &&
	     (CGAL::compare(q.x(), x_min) == LARGER)     ) ||
	    ((CGAL::compare(r.y(), y_min) == EQUAL)   &&
	     (CGAL::compare(r.x(), x_max) == SMALLER) &&
	     (CGAL::compare(r.x(), x_min) == LARGER)     )   )
        { // grow rectangle upwards
          CGAL_SDG_DEBUG(std::cout
              << "debug vring grow upwards" << std::endl;);
          y_max = y_min + x_max - x_min;
        } else
        { // grow rectangle downwards
          CGAL_SDG_DEBUG(std::cout
              << "debug vring grow downwards" << std::endl;);
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
  const
  {
    CGAL_precondition( p.is_point() && q.is_segment() &&
		       r.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug: compute_pss entering p=" << p
       << " q=" << q << " r=" << r << std::endl;);

    const bool pq =
      same_points(p, q.source_site()) || same_points(p, q.target_site());
    const bool pr =
      same_points(p, r.source_site()) || same_points(p, r.target_site());

    if ( pq && pr ) {
      // philaris: result should be point p
      const Point_2 pp = p.point();
      ux_ = pp.hx();
      uy_ = pp.hy();
      uz_ = pp.hw();
      return;
    }
    const bool is_q_hor = is_site_horizontal(q);
    const bool is_q_ver = is_site_vertical(q);
    const bool is_r_hor = is_site_horizontal(r);
    const bool is_r_ver = is_site_vertical(r);
    const bool is_q_hv = is_q_hor || is_q_ver;
    const bool is_r_hv = is_r_hor || is_r_ver;
    if (is_q_hv && is_r_hv) {
      compute_pss_both_hv(p, q, r, is_q_hor, is_r_hor, pq, pr);
    } else {
      if (pq || pr) {
        compute_pss_endp(p, q, r,
            is_q_hv, is_q_hor, pq, is_r_hv, is_r_hor, pr);
      } else {
        compute_pss_nonendp(p, q, r,
            is_q_hv, is_q_hor, is_r_hv, is_r_hor);
      }
    }
  }

  // PSS case when not both segments are axis-parallel and p is
  // not an endpoint of any of the segments q and r
  inline void
  compute_pss_nonendp(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hv, const bool is_q_hor,
      const bool is_r_hv, const bool is_r_hor)
  const
  {
    CGAL_USE(is_q_hv);
    CGAL_USE(is_q_hor);
    CGAL_USE(is_r_hv);
    CGAL_USE(is_r_hor);
    const Line_2 lq = orient_line_nonendp(p, q);
    const Line_2 lr = orient_line_nonendp(p, r);
    const unsigned int bq = bearing(lq);
    const unsigned int br = bearing(lr);
    const unsigned int bdiff = bearing_diff(bq, br);
    CGAL_assertion( bdiff != 0 );
    CGAL_assertion( bdiff != 7 );

    if (bdiff == 1) {
      compute_pss_corner_and_pt(p, q, r, lq, lr, bq, br);
    } else if (bdiff == 2) {
      compute_pss_nonhv_consecutive(p, q, r, lq, lr, bq, br);
    } else if ((bdiff == 3) || (bdiff == 4)) {
      compute_pss_ortho_wedge(p, q, r, lq, lr, bq, br);
    } else if (bdiff == 5) {
      compute_pss_side_p_known(p, q, r, lq, lr, bq, br);
    } else if (bdiff == 6) {
      compute_pss_lines_side(p, lq, lr, (br+1)%8);
    } else {
      CGAL_assertion( false );
    }
    CGAL_assertion_code( is_v_computed = true );
    CGAL_assertion( oriented_side_of_line(lq, this->point()) != ZERO );
    CGAL_assertion( oriented_side_of_line(lr, this->point()) != ZERO );
  }

  inline void
  compute_pss_lines_side(const Site_2& p,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bside)
  const
  {
    CGAL_precondition(bside % 2 == 1);
    const bool side_ver = (bside % 4 == 1);
    const FT pcoord = (side_ver) ? p.point().x() : p.point().y();
    const FT qcoord = coord_at(lq, pcoord, side_ver);
    const FT rcoord = coord_at(lr, pcoord, side_ver);
    const FT sidelen = CGAL::abs(qcoord-rcoord);
    const int sgn = (bside < 4) ? -1 : +1;
    const RT two(2);
    ux_ = side_ver ? (two*pcoord + sgn*sidelen) : (qcoord+rcoord);
    uy_ = side_ver ? (qcoord+rcoord) : (two*pcoord + sgn*sidelen);
    uz_ = two;
  }

  inline void
  compute_pss_side_p_known(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br)
  const
  {
    CGAL_USE(q);
    CGAL_USE(r);
    CGAL_USE(bq);
    const Bearing bside = (br + ((br % 2 == 0) ? 1 : 2)) % 8;
    const bool l_compute_y = (bside % 4 == 1) ? true : false;
    const FT pcoord = l_compute_y ? p.point().x() : p.point().y();
    const FT qcoord = coord_at(lq, pcoord, l_compute_y);
    const FT rcoord = coord_at(lr, pcoord, l_compute_y);
    const Point_2 qcorner =
      l_compute_y ? Point_2(pcoord, qcoord) : Point_2(qcoord, pcoord);
    const Point_2 rcorner =
      l_compute_y ? Point_2(pcoord, rcoord) : Point_2(rcoord, pcoord);
    const Point_2 v =
      center_from_same_side_corners(rcorner, qcorner, bside);
    ux_ = v.hx();
    uy_ = v.hy();
    uz_ = v.hw();
  }

  inline void
  compute_pss_ortho_wedge(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br)
  const
  {
    CGAL_USE(q);
    CGAL_USE(r);
    const FT xp = p.point().x();
    const FT yp = p.point().y();
    const bool lq_compute_y = ((bq / 2) % 2 == 0) ? false : true;
    const FT & lq_from_p = lq_compute_y ? xp : yp;
    const FT & lr_from_p = lq_compute_y ? yp : xp;
    const FT qcoord = coord_at(lq, lq_from_p, lq_compute_y);
    const FT rcoord = coord_at(lr, lr_from_p, ! lq_compute_y);
    const FT qdist = (bq < 4) ? qcoord - lr_from_p :
                                lr_from_p - qcoord;
    CGAL_assertion(CGAL::sign(qdist) == POSITIVE);
    const FT rdist = (bq <= 1) || (bq >= 6) ? rcoord - lq_from_p :
                                              lq_from_p - rcoord;
    CGAL_assertion(CGAL::sign(rdist) == POSITIVE);
    const Comparison_result cmpqr = CGAL::compare(qdist, rdist);
    const bool q_closer = (cmpqr == SMALLER);
    const Point_2 corner =
      q_closer ?
      (lq_compute_y ? Point_2(xp, qcoord) : Point_2(qcoord, yp)) :
      (lq_compute_y ? Point_2(rcoord, yp) : Point_2(xp, rcoord)) ;
    const Bearing bnonhv = (bq % 2 == 1) ? br : bq;
    CGAL_assertion(bnonhv % 2 == 0);
    const Line_2 lcorner = (bnonhv % 4 == 0)?
        compute_neg_45_line_at(corner) :
        compute_pos_45_line_at(corner) ;
    const Line_2 & lother = q_closer ? lr : lq;
    RT hx, hy, hw;
    compute_intersection_of_lines(lother, lcorner, hx, hy, hw);
    const Point_2 v =
      center_from_opposite_corners(Point_2(hx, hy, hw), corner);
    ux_ = v.hx();
    uy_ = v.hy();
    uz_ = v.hw();
  }

  inline void
  compute_pss_nonhv_consecutive(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br)
  const
  {
    const Bearing bqr = (bq+1)%8;
    return (bqr % 4) == 1 ?
      compute_pss_x_consecutive(p, q, r, lq, lr, bq, br, bqr) :
      compute_pss_y_consecutive(p, q, r, lq, lr, bq, br, bqr) ;
  }

  inline void
  compute_pss_x_consecutive(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br,
      const Bearing bqr)
  const
  {
    CGAL_precondition((bqr == 1) || (bqr == 5));
    CGAL_USE(q);
    CGAL_USE(r);
    CGAL_USE(bq);
    CGAL_USE(br);
    const FT xp = p.point().x();
    const FT x =
      (lr.b()*(lq.b()*xp + lq.c()) - lq.b()*lr.c()) /
      (lr.b()*(lq.b() -lq.a()) + lq.b()*lr.a()) ;
    const FT yq = (lq.a()*x + lq.c())/(-lq.b());
    const FT yr = (lr.a()*x + lr.c())/(-lr.b());

    const FT yp = p.point().y();
    if (CGAL::compare(yp, yq) == ((bqr == 1) ? SMALLER : LARGER)) {
      // p close to q
      const FT xs = coord_at(lq, yp, false);
      const FT ys = coord_at(lr, xs, true);
      ux_ = RT(2)*xs + (yp - ys);
      uy_ = yp + ys;
    } else if (CGAL::compare(yp, yr) == ((bqr == 1) ? LARGER : SMALLER)) {
      // p close to r
      const FT xs = coord_at(lr, yp, false);
      const FT ys = coord_at(lq, xs, true);
      ux_ = RT(2)*xs + (ys - yp);
      uy_ = yp + ys;
    } else {
      // p on opposite side of two lines (or on its corners)
      ux_ = xp + x;
      uy_ = yq + yr;
    }
    uz_ = RT(2);
  }

  inline void
  compute_pss_y_consecutive(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br,
      const Bearing bqr)
  const
  {
    CGAL_precondition((bqr == 3) || (bqr == 7));
    CGAL_USE(q);
    CGAL_USE(r);
    CGAL_USE(bq);
    CGAL_USE(br);
    const FT yp = p.point().y();
    const FT y =
      (lr.a()*(lq.a()*yp - lq.c()) + lq.a()*lr.c()) /
      (lr.a()*(lq.a() + lq.b()) - lq.a()*lr.b()) ;
    const FT xq = (lq.b()*y + lq.c())/(-lq.a());
    const FT xr = (lr.b()*y + lr.c())/(-lr.a());

    const FT xp = p.point().x();
    if (CGAL::compare(xp, xq) == ((bqr == 3) ? LARGER : SMALLER)) {
      // p close to q
      const FT ys = coord_at(lq, xp, true);
      const FT xs = coord_at(lr, ys, false);
      ux_ = xp + xs;
      uy_ = RT(2)*ys + (xs - xp);
    } else if (CGAL::compare(xp, xr) == ((bqr == 3) ? SMALLER : LARGER)) {
      // p close to r
      const FT ys = coord_at(lr, xp, true);
      const FT xs = coord_at(lq, ys, false);
      ux_ = xp + xs;
      uy_ = RT(2)*ys + (xp - xs);
    } else {
      // p on opposite side of two lines (or on its corners)
      ux_ = xq + xr;
      uy_ = yp + y;
    }
    uz_ = RT(2);
  }

  inline void
  compute_pss_corner_and_pt(const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br)
  const
  {
    const Bearing cb = (bq % 2 == 0) ? bq : br;
    Point_2 v;
    const bool is_qsrc_r = is_endpoint_of(q.source_site(), r);
    if (is_qsrc_r) {
      v = center_from_corner_and_pt(q.source(), cb, p.point());
    } else {
      const bool is_qtrg_r = is_endpoint_of(q.target_site(), r);
      if (is_qtrg_r) {
        v = center_from_corner_and_pt(q.target(), cb, p.point());
      } else {
        RT cx, cy, cw;
        compute_intersection_of_lines(lq, lr, cx, cy, cw);
        v = center_from_corner_and_pt(Point_2(cx, cy, cw), cb, p.point());
      }
    }
    ux_ = v.hx();
    uy_ = v.hy();
    uz_ = v.hw();
  }


  // PSS case when not both segments are axis-parallel and p is
  // an endpoint of one of the segments
  inline void
  compute_pss_endp(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hv, const bool is_q_hor, const bool pq,
      const bool is_r_hv, const bool is_r_hor, const bool pr)
  const
  {
    CGAL_precondition(pq || pr);
    CGAL_USE(pr);
    const Line_2 lendp = orient_line_endp(p, (pq ? q : r), pq);
    const Line_2 lnon = orient_line_nonendp(p, (pq ? r : q));
    CGAL_SDG_DEBUG(std::cout << "debug compute_pss_endp lendp="
        << lendp.a() << ' ' << lendp.b() << ' ' << lendp.c() << ' '
        << std::endl ; );
    CGAL_SDG_DEBUG(std::cout << "debug compute_pss_endp lnon="
        << lnon.a() << ' ' << lnon.b() << ' ' << lnon.c() << ' '
        << std::endl ; );
    const Line_2 llbis = bisector_linf_line(
        (pq ? q : r), (pq ? r : q), lendp, lnon);
    Line_2 lperp;
    const bool is_hv = pq ? is_q_hv : is_r_hv;
    if (is_hv) {
      const bool is_hor = pq ? is_q_hor : is_r_hor;
      lperp = is_hor ? compute_ver_line_at(p.point()) :
                       compute_hor_line_at(p.point()) ;
    } else {
      lperp = has_positive_slope(pq ? q : r) ?
        compute_neg_45_line_at(p.point()) :
        compute_pos_45_line_at(p.point()) ;
    }
    CGAL_SDG_DEBUG(std::cout << "debug compute_pss_endp llbis="
        << llbis.a() << ' ' << llbis.b() << ' ' << llbis.c() << ' '
        << std::endl ; );
    CGAL_SDG_DEBUG(std::cout << "debug compute_pss_endp lperp="
        << lperp.a() << ' ' << lperp.b() << ' ' << lperp.c() << ' '
        << std::endl ; );
    compute_intersection_of_lines(llbis, lperp, ux_, uy_, uz_);
    CGAL_SDG_DEBUG(std::cout << "debug compute_pss_endp vertex="
        << ux_ << ' ' << uy_ << ' ' << uz_ << ' '
        << Point_2(ux_, uy_, uz_) << std::endl ; );
    CGAL_assertion_code( is_v_computed = true );
    CGAL_assertion( oriented_side_of_line(lendp, this->point()) != ZERO );
    CGAL_assertion( oriented_side_of_line(lnon, this->point()) != ZERO );
  }

  // both segments are axis-parallel
  inline void
  compute_pss_both_hv(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hor, const bool is_r_hor,
      const bool pq, const bool pr)
  const
  {
    CGAL_precondition(! (pq && pr));
    if (is_q_hor == is_r_hor) {
      // parallel segments
      const RT q_coord = hvseg_coord(q, is_q_hor);
      const RT r_coord = hvseg_coord(r, is_r_hor);
      RT & upar = is_q_hor ? ux_ : uy_;
      RT & uort = is_q_hor ? uy_ : ux_;
      upar = RT(2)*(is_q_hor ? p.point().x() : p.point().y())
        + (( pq || pr ) ? RT(0) :
                          RT(is_q_hor ? +1 : -1)*(r_coord - q_coord));
      uort = q_coord + r_coord;
      uz_ = RT(2);
    } else {
      return compute_pss_both_hv_nonpar(
          p, q, r, is_q_hor, is_r_hor, pq, pr);
    }
  }

  // one segment is horizontal and the other is vertical
  inline void
  compute_pss_both_hv_nonpar(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hor, const bool is_r_hor,
      const bool pq, const bool pr)
  const
  {
    CGAL_precondition(is_q_hor != is_r_hor);
    if (pq || pr) {
      const RT q_coord = hvseg_coord(q, is_q_hor);
      const RT r_coord = hvseg_coord(r, is_r_hor);
      const bool is_touched_hor = pq ? is_q_hor : is_r_hor;
      const RT coord_c = is_touched_hor ? p.point().x() : p.point().y();
      const RT radius = CGAL::abs(coord_c - (pq ? r_coord : q_coord));
      RT & upar = is_touched_hor ? ux_ : uy_;
      RT & uort = is_touched_hor ? uy_ : ux_;
      const Site_2 & sother =
        pq ? (same_points(p, q.source_site()) ?
              q.target_site() : q.source_site()) :
             (same_points(p, r.source_site()) ?
              r.target_site() : r.source_site());
      const bool test = is_touched_hor ?
        (scmpx(p, sother) == LARGER) : (scmpy(p, sother) == SMALLER);
      const RT sgn = RT( (pq ? +1: -1)* (test ? -1 : +1) );
      upar = coord_c;
      uort = (pq ? q_coord : r_coord) + sgn*radius;
      uz_ = RT(1);
      CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss vv="
          << Point_2(ux_, uy_, uz_) << " radius=" << radius << std::endl;);
      CGAL_assertion_code( const Point_2 pother = sother.point() );
      CGAL_assertion(pq ?
          CGAL::left_turn(Point_2(ux_,uy_,uz_), p.point(), pother) :
          CGAL::left_turn(pother, p.point(), Point_2(ux_,uy_,uz_)) );
      return;
    } else {
      return compute_pss_both_hv_nonpar_nonendp(
               p, q, r, is_q_hor, is_r_hor, pq, pr);
    }
  }

  // one segment is horizontal and the other is vertical and
  // the point p is not an endpoint of the segments
  inline void
  compute_pss_both_hv_nonpar_nonendp(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hor, const bool is_r_hor,
      const bool pq, const bool pr)
  const
  {
    CGAL_precondition(! (pq || pr));
    CGAL_USE(pq);
    CGAL_USE(pr);
    const RT q_coord = hvseg_coord(q, is_q_hor);
    const RT r_coord = hvseg_coord(r, is_r_hor);
    const RT p_coord_q = is_q_hor ? p.point().y() : p.point().x();
    const RT p_coord_r = is_r_hor ? p.point().y() : p.point().x();
    const RT sdistq = p_coord_q - q_coord;
    const RT sdistr = p_coord_r - r_coord;
    const RT distq = CGAL::abs(sdistq);
    const RT distr = CGAL::abs(sdistr);
    const RT & dx = is_r_hor ? distq : distr;
    const RT & dy = is_r_hor ? distr : distq;
    const Comparison_result cmp = CGAL::compare(dx, dy);
    if (cmp == LARGER) {
      ux_ = is_q_hor ? (r_coord + p_coord_r) : (q_coord + p_coord_q);
      uy_ = is_q_hor ? (RT(2)*q_coord + CGAL::sign(sdistq)*dx) :
                       (RT(2)*r_coord + CGAL::sign(sdistr)*dx) ;
    } else if (cmp == SMALLER) {
      uy_ = is_r_hor ? (r_coord + p_coord_r) : (q_coord + p_coord_q);
      ux_ = is_r_hor ? (RT(2)*q_coord + CGAL::sign(sdistq)*dy) :
                       (RT(2)*r_coord + CGAL::sign(sdistr)*dy) ;
    } else {
      ux_ = is_q_hor ? (r_coord + p_coord_r) : (q_coord + p_coord_q);
      uy_ = is_q_hor ? (q_coord + p_coord_q) : (r_coord + p_coord_r);
    }
    uz_ = RT(2);
  }

  inline void
  compute_pss_bisectors(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    const bool pq =
      same_points(p, q.source_site()) || same_points(p, q.target_site());
    Polychainline_2 goodbisector;
    if (pq) {
      goodbisector = bisector_linf(p, q);
      CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss bpq p=" << p
          << " q=" << q << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss bpq ="
          << goodbisector << std::endl;);
    } else {
      goodbisector = bisector_linf(r, p);
      CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss brp r=" << r
          << " p=" << p << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss brp ="
          << goodbisector << std::endl;);
    }

    Polychainline_2 bqr = bisector_linf(q, r);
    CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss bqr q=" << q
        << " r=" << r << std::endl;);
    CGAL_SDG_DEBUG(std::cout << "debug: vring compute_pss bqr ="
        << bqr << std::endl;);

    Point_2 vv = goodbisector.first_intersection_point_with(bqr);
    CGAL_SDG_DEBUG(std::cout
        << "debug: vring compute_pss vv=" << vv << std::endl;);

    ux_ = vv.hx();
    uy_ = vv.hy();
    uz_ = vv.hw();

  }



  //--------------------------------------------------------------------------

  inline void
  compute_pps_endp_hv(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r, const bool is_r_horizontal)
  const
  {
    CGAL_USE(r);
    const Site_2 & A = p_endp_r ? p : q;
    const Site_2 & B = p_endp_r ? q : p;
    const RT Apar = is_r_horizontal ? A.point().x() : A.point().y();
    const RT Aort = is_r_horizontal ? A.point().y() : A.point().x();
    const RT Bpar = is_r_horizontal ? B.point().x() : B.point().y();
    const RT Bort = is_r_horizontal ? B.point().y() : B.point().x();
    const RT dpar = Apar - Bpar;
    const RT dort = Aort - Bort;
    const RT absdpar = CGAL::abs(dpar);

    RT & upar = is_r_horizontal ? ux_ : uy_;
    RT & uort = is_r_horizontal ? uy_ : ux_;

    if (2*absdpar < CGAL::abs(dort)) {
      upar = RT(2)*Apar;
      uort = RT(2)*Aort - dort;
      uz_ = RT(2);
    } else {
      upar = Apar;
      uort = Aort - CGAL::sign(dort)*absdpar;
      uz_ = RT(1);
    }
  }

  inline void
  compute_pps_nonendp_hv_samecoord(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_r_horizontal)
  const
  {
    const RT ppar = is_r_horizontal ? p.point().x() : p.point().y();
    const RT port = is_r_horizontal ? p.point().y() : p.point().x();
    const RT qort = is_r_horizontal ? q.point().y() : q.point().x();
    RT & upar = is_r_horizontal ? ux_ : uy_;
    RT & uort = is_r_horizontal ? uy_ : ux_;
    const RT segort = (is_r_horizontal)?
      horseg_y_coord(r) : verseg_x_coord(r);
    const RT sumort = port + qort;
    uort = sumort;
    const int vhsign = is_r_horizontal ? +1 : -1;
    const int distsign = CGAL::abs(segort-qort) < CGAL::abs(segort-port) ?
      +1: -1;
    upar = RT(2)*ppar - RT(vhsign*distsign)*(RT(2)*segort-sumort);
    uz_ = RT(2);
  }

  /* compute pps vertex when the points p, q are not endpoints of
   * segment r, and when segment r is axis-parallel
   */
  inline void
  compute_pps_nonendp_hv(const Site_2& p, const Site_2& q, const Site_2& r,
                         const bool is_r_horizontal)
  const
  {
    if ((is_r_horizontal     && (scmpx(p, q) == EQUAL)) ||
        ((! is_r_horizontal) && (scmpy(p, q) == EQUAL))   ) {
      return compute_pps_nonendp_hv_samecoord(p, q, r, is_r_horizontal);
    }

    // here, the segment is axis-parallel and the two points:
    // either: are not sharing a coordinate
    // or: they share a coordinate AND
    //     the line connecting the two points is parallel to the segment
    const Point_2 pp = p.point();
    const Point_2 qq = q.point();
    Line_2 l = compute_supporting_line(r);
    if (oriented_side_of_line(l, pp) == NEGATIVE) {
      l = opposite_line(l);
    }
    CGAL_assertion(oriented_side_of_line(l, pp) == POSITIVE);
    CGAL_assertion(oriented_side_of_line(l, qq) == POSITIVE);

    const Comparison_result perpcomp =
      is_r_horizontal ? scmpy(p, q) : scmpx(p, q);

    const RT coordr = hvseg_coord(r, is_r_horizontal);

    RT & upar = is_r_horizontal ? ux_ : uy_;
    RT & uort = is_r_horizontal ? uy_ : ux_;

    if (perpcomp == EQUAL) {
      const RT pqdist = is_r_horizontal ?
        CGAL::abs(pp.x()-qq.x()) : CGAL::abs(pp.y()-qq.y());
      const RT signrdist = (is_r_horizontal ? pp.y() : pp.x()) - coordr;
      Comparison_result comp = CGAL::compare(pqdist, CGAL::abs(signrdist));
      upar = is_r_horizontal ? pp.x() + qq.x() : pp.y() + qq.y();
      if (comp == LARGER) {
        uort = RT(2)*coordr + CGAL::sign(signrdist)*pqdist;
      } else {
        uort = coordr + (is_r_horizontal ? pp.y() : pp.x());
      }
      uz_ = RT(2);
      return;
    }
    // here, perpcomp is not EQUAL
    const Sign signla = CGAL::sign(l.a());
    const Sign signlb = CGAL::sign(l.b());
    const Sign & testsign = is_r_horizontal ? signlb : signla;
    CGAL_assertion(testsign != ZERO);
    const Point_2 & farp =
      testsign == POSITIVE ?
        (perpcomp == SMALLER ? qq : pp) :
        (perpcomp == SMALLER ? pp : qq) ;
    CGAL_assertion(Base::compare_linf_distances_to_line(
          l, farp == pp ? qq : pp , farp) == SMALLER);
    const RT pqdist = (CGAL::max)(
      CGAL::abs(pp.x()-qq.x()), CGAL::abs(pp.y()-qq.y()));

    const RT sdistf = (is_r_horizontal ? farp.y() : farp.x()) - coordr;
    CGAL_assertion(CGAL::sign(sdistf) == testsign);

    if (CGAL::compare(CGAL::abs(sdistf), pqdist) == LARGER) {
      const bool is_p_farthest = farp == pp;
      const Point_2 & closep = (is_p_farthest)? qq : pp;
      upar = RT(2)* (is_r_horizontal ? closep.x() : closep.y()) +
        (is_r_horizontal? -1 : +1) * (is_p_farthest? -1 : +1) * sdistf;
      uort = RT(2)*coordr + sdistf;
    } else {
      upar = is_r_horizontal ? pp.x() + qq.x() : pp.y() + qq.y();
      uort = RT(2)*coordr + CGAL::sign(sdistf)*pqdist;
    }
    uz_ = RT(2);
    return;
  }

  inline void
  compute_pps_endp_slope(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r, const bool pos_slope)
  const
  {
    CGAL_USE(r);
    const Site_2 & A = p_endp_r ? p : q;
    const Site_2 & B = p_endp_r ? q : p;
    const RT Ax = A.point().x();
    const RT Ay = A.point().y();
    const RT Bx = B.point().x();
    const RT By = B.point().y();
    const RT dx = Ax - Bx;
    const RT dy = Ay - By;
    const RT absdx = CGAL::abs(dx);
    const RT absdy = CGAL::abs(dy);

    if (absdx > absdy) {
      ux_ = RT(2)*Ax - dx;
      uy_ = RT(2)*Ay - RT(pos_slope? -1 : +1)*dx;
    } else {
      ux_ = RT(2)*Ax - RT(pos_slope? -1 : +1)*dy;
      uy_ = RT(2)*Ay - dy;
    }
    uz_ = RT(2);
  }

  inline void
  compute_pps_endp(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r)
  const
  {
    const bool is_r_horizontal = is_site_horizontal(r);
    if (is_r_horizontal || is_site_vertical(r)) {
      return compute_pps_endp_hv(p, q, r, p_endp_r, is_r_horizontal);
    } else {
      const bool pos_slope = has_positive_slope(r);
      return compute_pps_endp_slope(p, q, r, p_endp_r, pos_slope);
    }
  }

  /* compute pps vertex when the points p, q are not endpoints of
   * segment r, when segment r is not axis-parallel, and when the
   * two points p and q do not share any coordinate
   */
  inline void
  compute_pps_nonendp_nonhv_nonsamec
  (const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    Line_2 l = compute_supporting_line(r);
    if (oriented_side_of_line(l, p.point()) == NEGATIVE) {
      l = opposite_line(l);
    }
    CGAL_assertion(oriented_side_of_line(l, p.point()) == POSITIVE);
    CGAL_assertion(oriented_side_of_line(l, q.point()) == POSITIVE);

    const bool pos_slope = has_positive_slope(r);
    const Comparison_result first_comp =
      (pos_slope) ? scmpy(p, q) : scmpx(p, q);
    const Comparison_result second_comp =
      (pos_slope) ? scmpx(p, q) : scmpy(p, q);

    const Sign signla = CGAL::sign(l.a());
    const Sign signlb = CGAL::sign(l.b());
    const Comparison_result first_value =
      (signlb == POSITIVE)? SMALLER : LARGER;

    const Comparison_result second_value =
      (signla == NEGATIVE)? SMALLER : LARGER;

    if (first_comp == first_value) {
      const RT pcoord = pos_slope ? p.point().x() : p.point().y();
      const RT lineval = coord_at(l, pcoord, pos_slope);
      const Point_2 corner = pos_slope?
        Point_2(pcoord, lineval) : Point_2(lineval, pcoord);
      const RT sidelen = (CGAL::max)(CGAL::abs(corner.x() - q.point().x()),
                                     CGAL::abs(corner.y() - q.point().y()));
      ux_ = RT(2)*corner.x() + signla*sidelen;
      uy_ = RT(2)*corner.y() + signlb*sidelen;
      uz_ = RT(2);
      return;
    }
    if (second_comp == second_value) {
      const RT qcoord = pos_slope ? q.point().y() : q.point().x();
      const RT lineval = coord_at(l, qcoord, ! pos_slope);
      const Point_2 corner = pos_slope?
        Point_2(lineval, qcoord) : Point_2(qcoord, lineval);
      const RT sidelen = (CGAL::max)(CGAL::abs(corner.x() - p.point().x()),
                                     CGAL::abs(corner.y() - p.point().y()));
      ux_ = RT(2)*corner.x() + signla*sidelen;
      uy_ = RT(2)*corner.y() + signlb*sidelen;
      uz_ = RT(2);
      return;
    }

    CGAL_assertion((first_comp  == -first_value ) &&
                   (second_comp == -second_value)    );

    const RT px = p.point().x();
    const RT py = p.point().y();
    const RT qx = q.point().x();
    const RT qy = q.point().y();
    const RT pqdist = (CGAL::max)(CGAL::abs(px - qx), CGAL::abs(py - qy));

    CGAL_SDG_DEBUG(std::cout
        << "debug: vring pqdist=" << pqdist << std::endl;);

    const RT & pcoord = pos_slope ? px : py;
    const RT plineval = coord_at(l, pcoord, pos_slope);
    const RT & pothercoord = pos_slope ? py : px;
    const RT plen = CGAL::abs(plineval -  pothercoord);
    CGAL_SDG_DEBUG(std::cout
        << "debug: vring plen=" << plen << std::endl;);
    if (CGAL::compare(pqdist, plen) != SMALLER) {
      // here, appropriate projection of p on supporting line of segment r
      // is shorter than Linf p, q distance
      const Point_2 corner = pos_slope?
        Point_2(pcoord, plineval) : Point_2(plineval, pcoord);
      ux_ = RT(2)*corner.x() + signla*pqdist;
      uy_ = RT(2)*corner.y() + signlb*pqdist;
      uz_ = RT(2);
      return;
    }

    const RT & qcoord = pos_slope ? qy : qx;
    const RT qlineval = coord_at(l, qcoord, ! pos_slope);
    const RT & qothercoord = pos_slope ? qx : qy;
    const RT qlen = CGAL::abs(qlineval -  qothercoord);
    CGAL_SDG_DEBUG(std::cout
        << "debug: vring qlen=" << qlen << std::endl;);
    if (CGAL::compare(pqdist, qlen) != SMALLER) {
      // here, appropriate projection of q on supporting line of segment r
      // is shorter than Linf p, q distance
      const Point_2 corner = pos_slope?
        Point_2(qlineval, qcoord) : Point_2(qcoord, qlineval);
      ux_ = RT(2)*corner.x() + signla*pqdist;
      uy_ = RT(2)*corner.y() + signlb*pqdist;
      uz_ = RT(2);
      return;
    }

    CGAL_assertion((pqdist < plen) && (pqdist < qlen));

    // here, compute corner opposite of corner on line of segment r
    const Point_2 opposite_corner = pos_slope ?
      Point_2(qx, py) : Point_2(px, qy);
    CGAL_SDG_DEBUG(std::cout << "debug: vring opposite_corner="
        << opposite_corner << std::endl;);

    const Point_2 corner =
      compute_linf_projection_nonhom(l, opposite_corner);

    ux_ = corner.x() + opposite_corner.x();
    uy_ = corner.y() + opposite_corner.y();
    uz_ = RT(2);
  }

  inline void
  compute_pps_nonendp_nonhv(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    const bool samexpq = scmpx(p, q) == EQUAL;
    const bool sameypq = (samexpq)? false : make_certain(scmpy(p, q) == EQUAL);
    if (! (samexpq || sameypq)) {
      return compute_pps_nonendp_nonhv_nonsamec(p, q, r);
    } else {
      // samexpq || sameypq
      CGAL_assertion(samexpq != sameypq);
      Line_2 l = compute_supporting_line(r);
      const FT common_coord = (samexpq) ? p.point().x() : p.point().y();
      const FT sumdiffpq = (samexpq) ?
        p.point().y() + q.point().y() :
        p.point().x() + q.point().x();
      const bool pos_slope = has_positive_slope(r);
      FT vsamecoord;
      if (touch_same_side(p, q, l, samexpq, pos_slope)) {
        vsamecoord = common_coord +
          (pos_slope? +1: -1)*
            (coord_at(l, common_coord, samexpq) - (sumdiffpq/FT(2)));
      } else {
        const FT closest_coord =
          (samexpq)? ((pos_slope)? q.point().y() : p.point().y()):
                     ((pos_slope)? p.point().x() : q.point().x());
        if (is_orth_dist_smaller_than_pt_dist(
              closest_coord, l, p, q, samexpq)) {
          vsamecoord =
            coord_at(l, closest_coord, sameypq) +
            (((samexpq) ? (q.point().y() - p.point().y()) :
                          (p.point().x() - q.point().x())  ) / FT(2)) ;
        } else {
          const Line_2 lc (
              (samexpq) ? 1 : (2 * ((pos_slope) ? +1 : -1)),
              (samexpq) ? (2 * ((pos_slope) ? +1 : -1)) : 1 ,
              ((pos_slope)? -1 : +1 ) * sumdiffpq - common_coord);
          RT hx, hy, hz;
          compute_intersection_of_lines(l, lc, hx, hy, hz);
          vsamecoord = ((samexpq ? hx/hz : hy/hz) + common_coord)/ FT(2);
        }
      }
      const FT vdiffcoord = sumdiffpq/FT(2);
      ux_ = (samexpq) ? vsamecoord : vdiffcoord;
      uy_ = (samexpq) ? vdiffcoord : vsamecoord;
      uz_ = RT(1);
    }
  }

  inline void
  compute_pps_nonendp(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    const bool is_r_horizontal = is_site_horizontal(r);
    if (is_r_horizontal || is_site_vertical(r)) {
      return compute_pps_nonendp_hv(p, q, r, is_r_horizontal);
    } else {
      return compute_pps_nonendp_nonhv(p, q, r);
    }
  }

  void
  compute_pps(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    CGAL_precondition( p.is_point() && q.is_point() &&
		       r.is_segment() );

    CGAL_SDG_DEBUG(std::cout
        << "debug: vring compute_pps entering p=" << p
        << " q=" << q << " r=" << r << std::endl;);

    bool p_endp_r = is_endpoint_of(p, r);
    bool q_endp_r = is_endpoint_of(q, r);

    if (p_endp_r || q_endp_r) {
      return compute_pps_endp(p, q, r, p_endp_r);
    }
    CGAL_assertion(are_in_same_open_halfspace_of(p, q, r));
    return compute_pps_nonendp(p, q, r);
  }


  void
  compute_pps_bisectors(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    CGAL_precondition( p.is_point() && q.is_point() &&
		       r.is_segment() );

    CGAL_SDG_DEBUG(std::cout
        << "debug: vring compute_pps_bisectors entering p=" << p
        << " q=" << q << " r=" << r << std::endl;);

    bool p_endp_r = is_endpoint_of(p, r);
    bool q_endp_r = is_endpoint_of(q, r);
    Polychainline_2 bpq = bisector_linf(p, q);
    CGAL_SDG_DEBUG(std::cout << "debug: bpq p="
        << p << " q=" << q << std::endl;);
    CGAL_SDG_DEBUG(std::cout << "debug: bpq =" << bpq << std::endl;);

    CGAL_assertion(! (p_endp_r && q_endp_r));

    bool samexpq = (scmpx(p, q) == EQUAL);
    bool sameypq = (scmpy(p, q) == EQUAL);

    bool samecoordpq = samexpq || sameypq ;

    Polychainline_2 goodbisector;
    if (p_endp_r) {
      goodbisector = bisector_linf(r, p);
      CGAL_SDG_DEBUG(std::cout << "debug: vring brp r="
          << r << " p=" << p << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: vring brp res="
          << goodbisector << std::endl;);
    } else if (q_endp_r) {
      goodbisector = bisector_linf(q, r);
      CGAL_SDG_DEBUG(std::cout << "debug: vring bqr q="
          << q << " r=" << r << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: vring bqr res="
          << goodbisector << std::endl;);
    } else if (samecoordpq) {
      CGAL_SDG_DEBUG(std::cout << "debug vring PPS samecoordpq"
          << std::endl;);

      // check which of points p, q is closer to segment r

      bool use_bqr;

      Line_2 l = compute_supporting_line(r.supporting_site());

      if (((CGAL::sign(l.a()) == ZERO) && sameypq) ||
          ((CGAL::sign(l.b()) == ZERO) && samexpq)   )  {
        // here l is horizontal or vertical and parallel to pq;
        // bqr or brp are equally good
        use_bqr = true;
      } else {
        // here l and segment are neither hor. nor ver.
        Point_2 proj;

        // philaris: tofix with RT
        FT projft, pft, qft;
        if (samexpq) {
          // compute vertical projection
          proj = compute_vertical_projection(l, p.point());
          projft = proj.y();
          pft = p.point().y();
          qft = q.point().y();

        } else {
          CGAL_assertion(sameypq);
          // compute horizontal projection
          proj = compute_horizontal_projection(l, p.point());
          projft = proj.x();
          pft = p.point().x();
          qft = q.point().x();
        }
        Comparison_result cpq, cqproj;
        cpq    = CGAL::compare(pft, qft);
        cqproj = CGAL::compare(qft, projft);
        if (cpq == cqproj) {
          use_bqr = true;
        } else {
          use_bqr = false;
        }
      } // end of case of neither hor nor ver segment

      if (use_bqr) {
        goodbisector = bisector_linf(q, r);
        CGAL_SDG_DEBUG(std::cout << "debug: vring bqr q="
            << q << " r=" << r << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug: vring bqr res="
            << goodbisector << std::endl;);
      } else {
        goodbisector = bisector_linf(r, p);
        CGAL_SDG_DEBUG(std::cout << "debug: vring brp r="
            << r << " p=" << p << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug: vring brp res="
            << goodbisector << std::endl;);
      }
    } else {
      goodbisector = bisector_linf(q, r);
      CGAL_SDG_DEBUG(std::cout << "debug: vring bqr q="
          << q << " r=" << r << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: vring bqr res="
          << goodbisector << std::endl;);
    }


    Point_2 vv = bpq.first_intersection_point_with(goodbisector);

    CGAL_SDG_DEBUG(std::cout << "debug: PPS returns with vv="
        << vv << std::endl;);

    ux_ = vv.hx();
    uy_ = vv.hy();
    uz_ = vv.hw();
  }

  //--------------------------------------------------------------------------

  void
  compute_sss(const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() &&
		       r.is_segment() );
    const bool is_psrc_q = is_endpoint_of(p.source_site(), q);
    const bool is_psrc_r = is_endpoint_of(p.source_site(), r);
    const bool is_ptrg_q = is_endpoint_of(p.target_site(), q);
    const bool is_ptrg_r = is_endpoint_of(p.target_site(), r);

    if (is_psrc_q && is_psrc_r) {
      ux_ = p.source().hx();
      uy_ = p.source().hy();
      uz_ = p.source().hw();
    } else if (is_ptrg_q && is_ptrg_r) {
      ux_ = p.target().hx();
      uy_ = p.target().hy();
      uz_ = p.target().hw();
    } else {
      // here, not all segments have a common point

      const bool is_p_hor = is_site_horizontal(p);
      const bool is_q_hor = is_site_horizontal(q);
      const bool is_r_hor = is_site_horizontal(r);

      const bool is_p_hv = is_p_hor || is_site_vertical(p);
      const bool is_q_hv = is_q_hor || is_site_vertical(q);
      const bool is_r_hv = is_r_hor || is_site_vertical(r);

      if (is_p_hv && is_q_hv && is_r_hv) {
        return compute_sss_hv(p, q, r, is_p_hor, is_q_hor, is_r_hor);
      }

      Line_2 lines[3];
      orient_lines_linf(p, q, r, lines);

      compute_sss_bisectors(p, q, r, lines);
      CGAL_assertion_code( is_v_computed = true );
      CGAL_assertion( oriented_side_of_line(lines[0], this->point()) != ZERO );
      CGAL_assertion( oriented_side_of_line(lines[1], this->point()) != ZERO );
      CGAL_assertion( oriented_side_of_line(lines[2], this->point()) != ZERO );
    }
  }

  // SSS: all sites are axis-parallel
  inline void
  compute_sss_hv(const Site_2 & p, const Site_2 & q, const Site_2 & r,
      const bool is_p_hor, const bool is_q_hor, const bool is_r_hor)
  const
  {
    CGAL_precondition(! (is_p_hor && is_q_hor && is_r_hor));
    CGAL_precondition(is_p_hor || is_q_hor || is_r_hor);
    const unsigned int num_hor =
      (is_p_hor ? 1 : 0) + (is_q_hor ? 1 : 0) + (is_r_hor ? 1 : 0);
    CGAL_assertion((num_hor == 1) || (num_hor == 2));
    const bool are_common_hor = num_hor == 2;
    const bool is_odd_hor = ! are_common_hor;

    const Site_2 & odd = (is_odd_hor) ?
      (is_p_hor ? p : (is_q_hor ? q : r)) :
      (is_p_hor ? (is_q_hor ? r : q) : p);
    CGAL_assertion( (! (num_hor == 1)) || is_site_horizontal(odd) );
    CGAL_assertion( (! (num_hor == 2)) || is_site_vertical(odd) );
    const Site_2 & prev = (is_odd_hor) ?
      (is_p_hor ? r : (is_q_hor ? p : q)) :
      (is_p_hor ? (is_q_hor ? q : p) : r);
    CGAL_assertion( (! (num_hor == 1)) || is_site_vertical(prev) );
    CGAL_assertion( (! (num_hor == 2)) || is_site_horizontal(prev) );
    const Site_2 & next = (is_odd_hor) ?
      (is_p_hor ? q : (is_q_hor ? r : p)) :
      (is_p_hor ? (is_q_hor ? p : r) : q);
    CGAL_assertion( (! (num_hor == 1)) || is_site_vertical(next) );
    CGAL_assertion( (! (num_hor == 2)) || is_site_horizontal(next) );

    const RT prevc = hvseg_coord(prev, are_common_hor);
    const RT nextc = hvseg_coord(next, are_common_hor);
    RT & umid = is_odd_hor ? ux_ : uy_;
    RT & udis = is_odd_hor ? uy_ : ux_;
    umid = prevc + nextc;
    udis = RT(2)*hvseg_coord(odd, is_odd_hor) +
           RT(are_common_hor ? +1 : -1) * (prevc - nextc);
    uz_ = RT(2);
  }

  inline void
  compute_sss_bisectors(const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2 lines[])
  const
  {
    CGAL_SDG_DEBUG(std::cout << "debug vring compute_sss_bisectors"
        << " p=" << p << " q=" << q  << " r=" << r << std::endl;);
    Line_2 bpq = bisector_linf_line(p, q, lines[0], lines[1]);
    Line_2 bqr = bisector_linf_line(q, r, lines[1], lines[2]);
    compute_intersection_of_lines(bpq, bqr, ux_, uy_, uz_);
  }

  inline void
  compute_sss_bisectors_old(
      const Site_2& p, const Site_2& q, const Site_2& r)
  const
  {
    CGAL_SDG_DEBUG(std::cout << "debug vring compute_sss_bisectors_old"
        << " p=" << p << " q=" << q  << " r=" << r << std::endl;);
    Polychainline_2 bpq = bisector_linf(p, q);
    Polychainline_2 bqr = bisector_linf(q, r);
    Point_2 vv = bpq.first_intersection_point_with(bqr);
    ux_ = vv.hx();
    uy_ = vv.hy();
    uz_ = vv.hw();
  }

  //--------------------------------------------------------------------------

  void
  analyze_vertex(const Site_2& s1, const Site_2& s2, const Site_2& s3)
  {
    if ( s1.is_point() && s2.is_point() && s3.is_point() ) {
      v_type = PPP;
    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      pps_idx = 1;
      v_type = PPS;
    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      pps_idx = 2;
      v_type = PPS;
    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      pps_idx = 0;
      v_type = PPS;
    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      v_type = PSS;
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      v_type = PSS;
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      v_type = PSS;
    } else {
      v_type = SSS;
    }
  }

  void
  compute_vertex(const Site_2& s1, const Site_2& s2, const Site_2& s3)
  const
  {
    CGAL_assertion( ! is_v_computed );
    CGAL_SDG_DEBUG(std::cout << "debug vring compute_vertex "
        << s1 << ' ' << s2 << ' ' << s3 << std::endl;);
    if ( v_type == PPP ) {
      compute_ppp(s1, s2, s3);

    } else if ( s1.is_segment() && s2.is_point() && s3.is_point() ) {
      compute_pps(s2, s3, s1);
    } else if ( s1.is_point() && s2.is_segment() && s3.is_point() ) {
      compute_pps(s3, s1, s2);
    } else if ( s1.is_point() && s2.is_point() && s3.is_segment() ) {
      compute_pps(s1, s2, s3);

    } else if ( s1.is_point() && s2.is_segment() && s3.is_segment() ) {
      compute_pss(s1, s2, s3);
    } else if ( s1.is_segment() && s2.is_point() && s3.is_segment() ) {
      compute_pss(s2, s3, s1);
    } else if ( s1.is_segment() && s2.is_segment() && s3.is_point() ) {
      compute_pss(s3, s1, s2);
    } else {
      compute_sss(s1, s2, s3);
    }
    is_v_computed = true;
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

    if (
      ( p_.is_segment() && is_endpoint_of(t, p_) ) ||
      ( q_.is_segment() && is_endpoint_of(t, q_) ) ||
      ( r_.is_segment() && is_endpoint_of(t, r_) )  )
    {
      use_result = true;
      return POSITIVE;
    }

    if (
      ( p_.is_segment() && is_on_hv_seg_line(t, p_) ) ||
      ( q_.is_segment() && is_on_hv_seg_line(t, q_) ) ||
      ( r_.is_segment() && is_on_hv_seg_line(t, r_) )  )
    {
      use_result = true;
      return POSITIVE;
    }

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
    if (is_endpoint_of(t, p_) || is_endpoint_of(t, q_) ||
        is_endpoint_of(t, r_) ) {
      use_result = true;
      return POSITIVE;
    }
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

    Bounded_side bs =
      side_of_bounded_square(p_.point(), q_.point(), r_.point(), t);

    switch(bs) {
      case ON_UNBOUNDED_SIDE:
        CGAL_SDG_DEBUG(std::cout
            << "debug incircle_p returns POSITIVE" << std::endl;);
        return POSITIVE;
      case ON_BOUNDED_SIDE:
        CGAL_SDG_DEBUG(std::cout
            << "debug incircle_p returns NEGATIVE" << std::endl;);
        return NEGATIVE;
      default:
        CGAL_SDG_DEBUG(std::cout
            << "debug incircle_p returns ZERO" << std::endl;);
        return ZERO;
    }
  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PPS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    compute_v_if_not_computed();

    CGAL_SDG_DEBUG(std::cout << "debug vring incircle_p_no_easy PPS p="
      << p_ << " q=" << q_  << " r=" << r_ << " t=" << st
      << std::endl;);

    Point_2 t = st.point();

    Point_2 pointref = p_ref().point();

    CGAL_SDG_DEBUG(std::cout
        << "debug vring incircle_p_no_easy PPS pointref="
        << pointref << std::endl;);

    RT vx = ux_ - pointref.x() * uz_;
    RT vy = uy_ - pointref.y() * uz_;

    RT Rs =
      (CGAL::max) ( CGAL::abs(vx), CGAL::abs(vy) );

    RT scalediffdvtx = ux_ - t.x() * uz_;
    RT scalediffdvty = uy_ - t.y() * uz_;

    RT Rs1 =
      (CGAL::max)(
        CGAL::abs(scalediffdvtx),
        CGAL::abs(scalediffdvty) );


    Sign crude = CGAL::sign(Rs1 - Rs);

    if (crude != ZERO) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout
          << "debug vring refining in incircle_p_no_easy PPS pqr=("
          << p_ << ", " << q_ << ", " << r_ << "), "
          << "t=" << t
          << std::endl;);
      // here crude == ZERO, so
      // we might have to refine

      // tocheck

      const cpp11::tuple<
        const Site_2 &, const Site_2 &, const Site_2 &> sites =
         r_.is_segment() ? sdg_tuple_maker(p_, q_, r_) :
        (p_.is_segment() ? sdg_tuple_maker(q_, r_, p_) :
                           sdg_tuple_maker(r_, p_, q_) );

      const Site_2 & p1 = cpp11::get<0>(sites);
      const Site_2 & p2 = cpp11::get<1>(sites);
      const Site_2 & s  = cpp11::get<2>(sites);

      const RT d_fine = (CGAL::min)(CGAL::abs(scalediffdvtx),
                                    CGAL::abs(scalediffdvty));
      for (size_t i = 0; i < 2; ++i) {
        const Site_2 & cur = (i == 0) ? p1 : p2;
        const Point_2 pref = cur.point();
        const RT scalediffdvpx = ux_ - pref.x() * uz_;
        const RT scalediffdvpy = uy_ - pref.y() * uz_;

        Comparison_result sidecmp = EQUAL;
        const bool p_t_samex =
          CGAL::compare(scalediffdvpx, scalediffdvtx) == EQUAL;
        const bool p_t_on_same_ver_side =
          (p_t_samex) &&
          (CGAL::compare(CGAL::abs(scalediffdvpx), Rs1) == EQUAL) ;
        if (p_t_on_same_ver_side) {
          sidecmp = CGAL::compare(d_fine, CGAL::abs(scalediffdvpy));
        }
        const bool p_t_samey =
          CGAL::compare(scalediffdvpy, scalediffdvty) == EQUAL;
        const bool p_t_on_same_hor_side =
          (p_t_samey) &&
          (CGAL::compare(CGAL::abs(scalediffdvpy), Rs1) == EQUAL) ;
        if (p_t_on_same_hor_side) {
          sidecmp = CGAL::compare(d_fine, CGAL::abs(scalediffdvpx));
        }
        CGAL_SDG_DEBUG(std::cout << "vring test with p=" << cur
            << ", sidecmp=" << sidecmp << std::endl; );
        if (sidecmp == SMALLER) {
          return NEGATIVE;
        } else if (sidecmp == LARGER) {
          return POSITIVE;
        }
      }

      const bool is_s_hor = is_site_horizontal(s);
      const bool is_s_ver = is_site_vertical(s);

      const bool is_p1_endp_of_s = is_endpoint_of(p1, s);
      const bool is_p2_endp_of_s = is_endpoint_of(p2, s);

      // check for p, q with same coordinate and r non-hv segment
      if ((! (is_s_hor || is_s_ver)) &&
          (! (is_p1_endp_of_s || is_p2_endp_of_s))
         ) {
        CGAL_SDG_DEBUG(std::cout << "debug vring seg=" << s
            << " is non-axis parallel"
            << " and no points are its endpoints" << std::endl;);
        const bool pqsamex = scmpx(p1, p2) == EQUAL;
        bool pqsamey (false);
        if (pqsamex) {
          CGAL_SDG_DEBUG(std::cout << "debug vring points have same x, "
              << " might be on same Linf vertical side"
              << std::endl;);
        } else {
          pqsamey = scmpy(p1, p2) == EQUAL;
          if (pqsamey) {
            CGAL_SDG_DEBUG(std::cout << "debug vring points have same y, "
                << " might be on same Linf horizontal side" << std::endl;);
          }
        }
        if (pqsamex || pqsamey) {
          Line_2 lr = compute_supporting_line(s.supporting_site());
          Homogeneous_point_2 rref = compute_linf_projection_hom(lr, point());
          if ( pqsamex && (CGAL::compare(rref.x(), pointref.x()) == EQUAL) ) {
            RT scalediffdvry = uy_ - rref.y() * uz_;
            if (CGAL::sign(scalediffdvry) == CGAL::sign(scalediffdvty)) {
              if (CGAL::compare(CGAL::abs(scalediffdvtx),
                                CGAL::abs(scalediffdvty)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of pqsamex and rref same x case
          if ( pqsamey && (CGAL::compare(rref.y(), pointref.y()) == EQUAL) ) {
            RT scalediffdvrx = ux_ - rref.x() * uz_;
            if (CGAL::sign(scalediffdvrx) == CGAL::sign(scalediffdvtx)) {
              if (CGAL::compare(CGAL::abs(scalediffdvty),
                                CGAL::abs(scalediffdvtx)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of pqsamey and rref same y case
        } // end of case: pqsamex || pqsamey
      } // end of non-hv segment r case with p, q non-endpoints of r

      return ZERO;

      //FT radius_fine = linf_fine_radius(vv, p, q, r, type);

      //return CGAL::compare(d_fine, radius_fine);

    }

  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, PSS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    compute_v_if_not_computed();
    const Point_2 t = st.point();

    const Point_2 pref = p_ref().point();

    const RT xref = pref.x();
    const RT yref = pref.y();

    const RT vx = ux_ - xref * uz_;
    const RT vy = uy_ - yref * uz_;

    const RT Rs =
      (CGAL::max) ( CGAL::abs(vx), CGAL::abs(vy) );

    const RT tx = t.x() ;
    const RT ty = t.y() ;

    const RT scalediffdvtx = ux_ - tx * uz_;
    const RT scalediffdvty = uy_ - ty * uz_;

    const RT Rs1 =
      (CGAL::max)(
        CGAL::abs(scalediffdvtx),
        CGAL::abs(scalediffdvty) );

    const Sign s_Q = CGAL::sign(Rs1 - Rs);

    if (s_Q != ZERO) {
      return s_Q;
    } else {
      CGAL_SDG_DEBUG(std::cout
          << "debug vring refining in incircle_p_no_easy PSS pqr=("
          << p_ << ", " << q_ << ", " << r_ << "), "
          << "t=" << t
          << std::endl;);
      // here crude == ZERO, so
      // we might have to refine

      // tocheck

      const cpp11::tuple<
        const Site_2 &, const Site_2 &, const Site_2 &> sites =
         p_.is_point() ? sdg_tuple_maker(p_, q_, r_) :
        (q_.is_point() ? sdg_tuple_maker(q_, r_, p_) :
                         sdg_tuple_maker(r_, p_, q_) );

      const Site_2 & pt_site = cpp11::get<0>(sites);
      const Site_2 & s1 = cpp11::get<1>(sites);
      const Site_2 & s2 = cpp11::get<2>(sites);

      const bool is_s1src_s2 = is_endpoint_of(s1.source_site(), s2);
      const bool is_s1trg_s2 = is_endpoint_of(s1.target_site(), s2);

      if (is_s1src_s2 || is_s1trg_s2) {
        if ((is_site_h_or_v(s1) && (! is_site_h_or_v(s2))) ||
            (is_site_h_or_v(s2) && (! is_site_h_or_v(s1)))   )
        {
          CGAL_SDG_DEBUG(std::cout << "debug vring "
              << "s1, s2 candidates" << std::endl; );
          if (is_site_horizontal(s1) || is_site_horizontal(s2)) {
            Site_2 s1test = is_s1src_s2?
                       (s1.source_site()):
                       (s1.target_site());
            if (scmpx(s1test, st)
                == EQUAL)
            {
              // return NEGATIVE or ZERO
              Point_2 s1ref =
                      (is_s1src_s2?
                       s1.source_site(): s1.target_site())
                      .point();
              RT scalediffdvs1y = uy_ - s1ref.y() * uz_;
              Comparison_result test =
                CGAL::compare(
                    CGAL::abs(scalediffdvty),
                    CGAL::abs(scalediffdvs1y));
              return (test == SMALLER) ? NEGATIVE : ZERO;
            }
          } else { // one of q, r is vertical
            if (scmpy(is_s1src_s2?
                      s1.source_site(): s1.target_site(), st)
                == EQUAL)
            {
              // return NEGATIVE or ZERO
              CGAL_SDG_DEBUG(std::cout << "debug vring "
                  << "vertical case" << std::endl; );
              Point_2 s1ref =
                      (is_s1src_s2?
                       s1.source_site(): s1.target_site())
                      .point();
              RT scalediffdvs1x = ux_ - s1ref.x() * uz_;
              CGAL_SDG_DEBUG(std::cout << "debug vring "
                  << "scalediffdvs1x=" << scalediffdvs1x
                  << " scalediffdvtx=" << scalediffdvtx << std::endl; );
              Comparison_result test =
                CGAL::compare(
                    CGAL::abs(scalediffdvtx),
                    CGAL::abs(scalediffdvs1x));
              return (test == SMALLER) ? NEGATIVE : ZERO;
            }
          }
        }
      }

      const RT d_fine = (CGAL::min)(CGAL::abs(scalediffdvtx),
                                    CGAL::abs(scalediffdvty));
      const Point_2 pref = pt_site.point();
      const RT scalediffdvpx = ux_ - pref.x() * uz_;
      const RT scalediffdvpy = uy_ - pref.y() * uz_;
      Comparison_result sidecmp = EQUAL;
      const bool p_t_samex =
        CGAL::compare(scalediffdvpx, scalediffdvtx) == EQUAL;
      const bool p_t_on_same_ver_side =
        (p_t_samex) &&
        (CGAL::compare(CGAL::abs(scalediffdvpx), Rs1) == EQUAL) ;
      if (p_t_on_same_ver_side) {
        sidecmp = CGAL::compare(d_fine, CGAL::abs(scalediffdvpy));
      }
      const bool p_t_samey =
        CGAL::compare(scalediffdvpy, scalediffdvty) == EQUAL;
      const bool p_t_on_same_hor_side =
        (p_t_samey) &&
        (CGAL::compare(CGAL::abs(scalediffdvpy), Rs1) == EQUAL) ;
      if (p_t_on_same_hor_side) {
        sidecmp = CGAL::compare(d_fine, CGAL::abs(scalediffdvpx));
      }
      CGAL_SDG_DEBUG(std::cout << "debug: PSS temporary sidecmp = "
          << sidecmp << std::endl;);
      if (sidecmp == SMALLER) {
        return NEGATIVE;
      } else if (sidecmp == LARGER) {
        return POSITIVE;
      }

      if (p_t_samex || p_t_samey) {
        return ZERO;
      }

      if (! is_site_h_or_v(s1)) {
        // here s1 is non-axis parallel
        // therefore, it touches the square at a corner
        if (points_inside_touching_sides_v(
              s1, pt_site, s2, st, this->point())) {
          return NEGATIVE;
        }
      }

      if (! is_site_h_or_v(s2)) {
        // here s2 is non-axis parallel
        // therefore, it touches the square at a corner
        if (points_inside_touching_sides_v(
              s2, pt_site, s1, st, this->point())) {
          return NEGATIVE;
        }
      }

      CGAL_SDG_DEBUG(std::cout
          << "debug vring PSS P return final ZERO"
          << std::endl;);
      return ZERO;
    }

  }

  //--------------------------------------------------------------------------

  Sign incircle_p_no_easy(const Site_2& st, SSS_Type ) const
  {
    CGAL_precondition( st.is_point() );
    compute_v_if_not_computed();

    Point_2 t = st.point();

    Line_2 l = compute_supporting_line(p_.supporting_site());
    Homogeneous_point_2 hp = compute_linf_projection_hom(l, point());

    RT dup =
      (CGAL::max)(CGAL::abs(ux_ - hp.x() * uz_),
                  CGAL::abs(uy_ - hp.y() * uz_));

    RT scalediffdvtx = ux_ - t.x() * uz_;
    RT scalediffdvty = uy_ - t.y() * uz_;

    RT dut =
      (CGAL::max)(
        CGAL::abs(scalediffdvtx),
        CGAL::abs(scalediffdvty) );

    Sign crude_sign = CGAL::sign(dut - dup);
    if (crude_sign != ZERO) {
      return crude_sign;
    } else {
      CGAL_SDG_DEBUG(std::cout
          << "debug vring refining in incircle_p_no_easy SSS pqr=("
          << p_ << ", " << q_ << ", " << r_ << "), "
          << "t=" << t
          << std::endl;);

      const Site_2 * s1_ptr;
      const Site_2 * s2_ptr;

      const Site_2 * s1ptr_arr[] = {&p_, &q_, &r_};
      const Site_2 * s2ptr_arr[] = {&q_, &r_, &p_};

      for (int i = 0; i < 3; i++) {
        s1_ptr = s1ptr_arr[i];
        s2_ptr = s2ptr_arr[i];

        CGAL_SDG_DEBUG(std::cout << "debug vring check for candidates"
            << "(s1, s2) = " << *s1_ptr << ", " << *s2_ptr << std::endl; );

        bool is_s1src_s2 = is_endpoint_of((*s1_ptr).source_site(), *s2_ptr);
        bool is_s1trg_s2 = is_endpoint_of((*s1_ptr).target_site(), *s2_ptr);

        if (is_s1src_s2 || is_s1trg_s2) {
          if ((is_site_h_or_v(*s1_ptr) && (! is_site_h_or_v(*s2_ptr))) ||
              (is_site_h_or_v(*s2_ptr) && (! is_site_h_or_v(*s1_ptr)))   )
          {
            CGAL_SDG_DEBUG(std::cout << "debug vring "
                << "s1, s2 candidates" << std::endl; );
            if (is_site_horizontal(*s1_ptr) || is_site_horizontal(*s2_ptr)) {
              Site_2 s1test = is_s1src_s2?
                ((*s1_ptr).source_site()):
                ((*s1_ptr).target_site());
              if (scmpx(s1test, st)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                Point_2 s1ref =
                  (is_s1src_s2?
                   (*s1_ptr).source_site(): (*s1_ptr).target_site())
                  .point();
                RT scalediffdvs1y = uy_ - s1ref.y() * uz_;
                Comparison_result test =
                  CGAL::compare(
                      CGAL::abs(scalediffdvty),
                      CGAL::abs(scalediffdvs1y));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            } else { // one of q, r is vertical
              if (scmpy(is_s1src_s2?
                    (*s1_ptr).source_site(): (*s1_ptr).target_site(), st)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                CGAL_SDG_DEBUG(std::cout << "debug vring "
                    << "vertical case" << std::endl; );
                Point_2 s1ref =
                  (is_s1src_s2?
                   (*s1_ptr).source_site(): (*s1_ptr).target_site())
                  .point();
                RT scalediffdvs1x = ux_ - s1ref.x() * uz_;
                CGAL_SDG_DEBUG(std::cout << "debug vring "
                    << "scalediffdvs1x=" << scalediffdvs1x
                    << " scalediffdvtx=" << scalediffdvtx << std::endl; );
                Comparison_result test =
                  CGAL::compare(
                      CGAL::abs(scalediffdvtx),
                      CGAL::abs(scalediffdvs1x));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            }
          }
        }
      } // end for





      return ZERO;
    }
  }



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& t) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug: entering vring incircle_p with "
      << "v_type=" << v_type << " p="
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


  Sign incircle(const Line_2& l, PPP_Type) const
  {

    Point_2 pref = p_ref().point();
    Homogeneous_point_2 hp = compute_linf_projection_hom(l, point());

    CGAL_SDG_DEBUG(std::cout << "debug incircle l PPP: pref="
      << pref << std::endl;);

    RT dul = (CGAL::max)(
        CGAL::abs(ux_ - hp.x() * uz_),
        CGAL::abs(uy_ - hp.y() * uz_));

    RT dupref = (CGAL::max)(
        CGAL::abs(ux_ - pref.x() * uz_),
        CGAL::abs(uy_ - pref.y() * uz_));

    Comparison_result cr = CGAL::compare(dul, dupref);

    if ( cr == LARGER ) { return POSITIVE; }
    if ( cr == SMALLER ) { return NEGATIVE; }

    // here cr == EQUAL == ZERO, so
    // we might have to refine
    CGAL_SDG_DEBUG(std::cout
      << "debug vring refining in incircle l PPP pqr=("
      << p_ << ", " << q_ << ", " << r_ << "), "
      << "hp(x,y)=" << hp.x() << ' ' << hp.y()
      << ", l: " << l.a() << ' ' << l.b() << ' ' <<  l.c()
      << ", u(x,y,z)= " << ux_ << ' ' << uy_ << ' ' << uz_
      << std::endl;);

    Comparison_result other = linf_refine(l, hp);

    if (cr != other) {
      CGAL_SDG_DEBUG(std::cout
          << "incircle l PPP instead of 0 returning " << other
          << std::endl;);
    }

    return other;

  }

  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, PPS_Type) const
  {
    Point_2 pref = p_ref().point();

    CGAL_SDG_DEBUG(std::cout << "debug incircle l PPS: pref="
      << pref << std::endl;);

    RT vx = ux_ - pref.x() * uz_;
    RT vy = uy_ - pref.y() * uz_;

    RT dupref = (CGAL::max)(CGAL::abs(vx), CGAL::abs(vy));

    Homogeneous_point_2 hp = compute_linf_projection_hom(l, point());

    RT dul = (CGAL::max)(
        CGAL::abs(ux_ - hp.x() * uz_),
        CGAL::abs(uy_ - hp.y() * uz_));

    Sign cr = CGAL::sign(dul - dupref);

    if (cr != ZERO) {
      return cr;
    }

    // here cr == EQUAL == ZERO, so
    // we might have to refine
    CGAL_SDG_DEBUG(std::cout
      << "debug vring refining in incircle l PPS pqr=("
      << p_ << ", " << q_ << ", " << r_ << "), "
      << "hp(x,y)=" << hp.x() << ' ' << hp.y()
      << ", l: " << l.a() << ' ' << l.b() << ' ' <<  l.c()
      << ", u(x,y,z)= " << ux_ << ' ' << uy_ << ' ' << uz_
      << std::endl;);

    Comparison_result other = linf_refine(l, hp);

    if (cr != other) {
      CGAL_SDG_DEBUG(std::cout
          << "incircle l PPS instead of 0 returning " << other
          << std::endl;);
    }

    return other;
  }


  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, PSS_Type) const
  {
    Point_2 pref = p_ref().point();

    CGAL_SDG_DEBUG(std::cout << "debug incircle l PSS: pref="
      << pref << std::endl;);

    RT vx = ux_ - (pref.x() ) * uz_;
    RT vy = uy_ - (pref.y() ) * uz_;

    RT dupref = (CGAL::max)(CGAL::abs(vx), CGAL::abs(vy));

    Homogeneous_point_2 lhp = compute_linf_projection_hom(l, point());

    RT dul = (CGAL::max)(
        CGAL::abs(ux_ - lhp.x() * uz_),
        CGAL::abs(uy_ - lhp.y() * uz_));

    Sign cr = CGAL::sign(dul - dupref);

    if (cr != ZERO) {
      return cr;
    }

    // here cr == EQUAL == ZERO, so
    // we might have to refine
    CGAL_SDG_DEBUG(std::cout
      << "debug vring refining in incircle l PSS pqr=("
      << p_ << ", " << q_ << ", " << r_ << "), "
      << "hp(x,y)=" << lhp.x() << ' ' << lhp.y()
      << ", l: " << l.a() << ' ' << l.b() << ' ' <<  l.c()
      << ", u(x,y,z)= " << ux_ << ' ' << uy_ << ' ' << uz_
      << std::endl;);

    Comparison_result other = linf_refine(l, lhp);

    if (cr != other) {
      CGAL_SDG_DEBUG(std::cout
          << "incircle l PSS instead of 0 returning " << other
          << std::endl;);
    }

    return other;

  }


  //--------------------------------------------------------------------------

  Sign incircle(const Line_2& l, SSS_Type) const
  {
    Line_2 lref = compute_supporting_line(p_.supporting_site());
    Homogeneous_point_2 lrefhp =
      compute_linf_projection_hom(lref, point());
    RT dulref = (CGAL::max)(
        CGAL::abs(ux_ - lrefhp.x() * uz_),
        CGAL::abs(uy_ - lrefhp.y() * uz_));

    Homogeneous_point_2 lhp = compute_linf_projection_hom(l, point());
    RT dul = (CGAL::max)(
        CGAL::abs(ux_ - lhp.x() * uz_),
        CGAL::abs(uy_ - lhp.y() * uz_));

    Sign cr = CGAL::sign(dul - dulref);

    if (cr != ZERO) {
      return cr;
    }

    // here cr == EQUAL == ZERO, so
    // we might have to refine
    CGAL_SDG_DEBUG(std::cout
      << "debug vring refining in incircle l PSS pqr=("
      << p_ << ", " << q_ << ", " << r_ << "), "
      << "lhp(x,y)=" << lhp.x() << ' ' << lhp.y()
      << ", l: " << l.a() << ' ' << l.b() << ' ' <<  l.c()
      << ", u(x,y,z)= " << ux_ << ' ' << uy_ << ' ' << uz_
      << std::endl;);

    Comparison_result other = linf_refine(l, lhp);

    if (cr != other) {
      CGAL_SDG_DEBUG(std::cout
          << "incircle l SSS instead of 0 returning " << other
          << std::endl;);
    }

    return other;
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

    CGAL_SDG_DEBUG(std::cout
        << "debug vring incircle_s: about to return retval of"
        << " incircle_s_no_easy = " << retval << std::endl;);

    return retval;
  }

  template<class Type>
  Sign incircle_s_no_easy(const Site_2& t, Type type) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug vring fn incircle_s_no_easy pqrt= ("
        << p_ << ") (" << q_ << ") (" << r_ << ") (" << t << ")"
        << std::endl;);

    bool is_p_point = p_.is_point();
    bool is_q_point = q_.is_point();
    bool is_r_point = r_.is_point();

    unsigned int numpts_in_pqr =
      ((is_p_point)? 1 : 0) +
      ((is_q_point)? 1 : 0) +
      ((is_r_point)? 1 : 0)  ;

    bool is_p_tsrc(false);
    CGAL_assertion_code( bool has_p_endp_tsrc(false); )
    if ( is_p_point ) {
      if ( same_points(p_, t.source_site()) ) {
        is_p_tsrc = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // p is segment
      if (same_points(p_.source_site(), t.source_site()) ||
          same_points(p_.target_site(), t.source_site())   ) {
        has_p_endp_tsrc = true;
      }
    }
#endif

    bool is_q_tsrc(false);
    CGAL_assertion_code( bool has_q_endp_tsrc(false); )
    if ( is_q_point ) {
      if ( same_points(q_, t.source_site()) ) {
        is_q_tsrc = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // q is segment
      if (same_points(q_.source_site(), t.source_site()) ||
          same_points(q_.target_site(), t.source_site())   ) {
        has_q_endp_tsrc = true;
      }
    }
#endif

    bool is_r_tsrc(false);
    CGAL_assertion_code( bool has_r_endp_tsrc(false); )
    if ( is_r_point ) {
      if ( same_points(r_, t.source_site()) ) {
        is_r_tsrc = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // r is segment
      if (same_points(r_.source_site(), t.source_site()) ||
          same_points(r_.target_site(), t.source_site())   ) {
        has_r_endp_tsrc = true;
      }
    }
#endif

#ifndef CGAL_NO_ASSERTIONS
    unsigned int num_common_endp_tsrc =
      ((has_p_endp_tsrc)? 1 : 0) +
      ((has_q_endp_tsrc)? 1 : 0) +
      ((has_r_endp_tsrc)? 1 : 0)  ;
    CGAL_USE(num_common_endp_tsrc);
    CGAL_SDG_DEBUG(
    std::cout << "debug num_common_endp_tsrc="
      << num_common_endp_tsrc << std::endl;
    );
#endif

    unsigned int numendpts_of_t = 0;

    Sign d1, d2;
    if ( is_p_tsrc || is_q_tsrc || is_r_tsrc ) {
      d1 = ZERO;
      ++numendpts_of_t;
    } else {
      d1 = incircle_p(t.source_site());
    }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy d1="
        << d1 << " with tsrc=" << t.source_site() << std::endl;);

    if ( d1 == NEGATIVE ) { return NEGATIVE; }


    bool is_p_ttrg(false);
    CGAL_assertion_code( bool has_p_endp_ttrg(false); )
    if ( is_p_point ) {
      if ( same_points(p_, t.target_site()) ) {
        is_p_ttrg = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // p is segment
      if (same_points(p_.source_site(), t.target_site()) ||
          same_points(p_.target_site(), t.target_site())   ) {
        has_p_endp_ttrg = true;
      }
    }
#endif

    bool is_q_ttrg(false);
    CGAL_assertion_code( bool has_q_endp_ttrg(false); )
    if ( is_q_point ) {
      if ( same_points(q_, t.target_site()) ) {
        is_q_ttrg = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // q is segment
      if (same_points(q_.source_site(), t.target_site()) ||
          same_points(q_.target_site(), t.target_site())   ) {
        has_q_endp_ttrg = true;
      }
    }
#endif

    bool is_r_ttrg(false);
    CGAL_assertion_code( bool has_r_endp_ttrg(false); )
    if ( is_r_point ) {
      if ( same_points(r_, t.target_site()) ) {
        is_r_ttrg = true;
      }
    }
#ifndef CGAL_NO_ASSERTIONS
    else { // r is segment
      if (same_points(r_.source_site(), t.target_site()) ||
          same_points(r_.target_site(), t.target_site())   ) {
        has_r_endp_ttrg = true;
      }
    }
#endif

#ifndef CGAL_NO_ASSERTIONS
    unsigned int num_common_endp_ttrg =
      ((has_p_endp_ttrg)? 1 : 0) +
      ((has_q_endp_ttrg)? 1 : 0) +
      ((has_r_endp_ttrg)? 1 : 0)  ;
    CGAL_USE(num_common_endp_ttrg);
    CGAL_SDG_DEBUG(
    std::cout << "debug num_common_endp_ttrg="
      << num_common_endp_ttrg << std::endl;
    );
#endif

    if ( is_p_ttrg || is_q_ttrg || is_r_ttrg ) {
      d2 = ZERO;
      ++numendpts_of_t;
    } else {
      d2 = incircle_p(t.target_site());
    }
    if ( d2 == NEGATIVE ) { return NEGATIVE; }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy d2="
        << d2 << std::endl;);

    CGAL_assertion(numendpts_of_t < 2);

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy numendpts_of_t= "
      << numendpts_of_t << std::endl;);

    compute_v_if_not_computed();

    if (numendpts_of_t > 0) {
      bool is_t_horizontal = is_site_horizontal(t);
      bool is_t_vertical   = is_site_vertical(t);

      if (is_t_horizontal || is_t_vertical) {
        CGAL_assertion(numendpts_of_t == 1);

        // set endp to endpoint in {p,q,r}
        Site_2 endp;
        if ( is_p_tsrc || is_q_tsrc || is_r_tsrc ) {
          endp = t.source_site();
        } else {
          endp = t.target_site();
        }

        // numothers will be the number of segments
        // in {p,q,r} that have endp as an endpoint
        unsigned int numothers = 0;

        // a possible segment in {p,q,r} which has endpoint endp
        Site_2 other;

        // if there is a segment in {p,q,r}, try its endpoints
        if (numpts_in_pqr < 3) {

          if ((! is_p_point) && is_endpoint_of(endp, p_)) {
            numothers++;
            other = p_;
          }

          if ((! is_q_point) && is_endpoint_of(endp, q_)) {
            numothers++;
            other = q_;
          }

          if ((! is_r_point) && is_endpoint_of(endp, r_)) {
            numothers++;
            other = r_;
          }

        } // end of case: numpts_in_pqr < 3

        CGAL_assertion(numothers < 2);

        if (numothers == 1) {
          bool is_other_horizontal = is_site_horizontal(other);
          bool is_other_vertical = is_site_vertical(other);

          if ((is_t_horizontal && is_other_horizontal) ||
              (is_t_vertical && is_other_vertical)       ) {
            return POSITIVE;
          }
        } else {
          CGAL_assertion(numothers == 0);
          Point_2 vv(ux_, uy_,uz_);

          Comparison_result ptcmpxve =
            CGAL::compare(vv.x(), endp.point().x());
          Comparison_result ptcmpyve =
            CGAL::compare(vv.y(), endp.point().y());

          CGAL_SDG_DEBUG(std::cout << "debug vv = " << vv << std::endl;);

          if ( ( (ptcmpxve == EQUAL) && is_t_horizontal ) ||
               ( (ptcmpyve == EQUAL) && is_t_vertical   )    ) {
            return ZERO;
          }

        } // end of case numothers == 0
      }  // endif (is_t_horizontal || is_t_vertical)
    } // endif ((numendpts_of_t > 0) && (numpts_in_pqr < 3))

    bool is_tsrc_endp_of_p (false);
    bool is_tsrc_endp_of_q (false);
    bool is_tsrc_endp_of_r (false);
    bool is_ttrg_endp_of_p (false);
    bool is_ttrg_endp_of_q (false);
    bool is_ttrg_endp_of_r (false);

    if (! is_p_point) {
      is_tsrc_endp_of_p = same_points(t.source_site(), p_.source_site())
                       || same_points(t.source_site(), p_.target_site());
      is_ttrg_endp_of_p = same_points(t.target_site(), p_.source_site())
                       || same_points(t.target_site(), p_.target_site());
    }
    if (! is_q_point) {
      is_tsrc_endp_of_q = same_points(t.source_site(), q_.source_site())
                       || same_points(t.source_site(), q_.target_site());
      is_ttrg_endp_of_q = same_points(t.target_site(), q_.source_site())
                       || same_points(t.target_site(), q_.target_site());
    }
    if (! is_r_point) {
      is_tsrc_endp_of_r = same_points(t.source_site(), r_.source_site())
                       || same_points(t.source_site(), r_.target_site());
      is_ttrg_endp_of_r = same_points(t.target_site(), r_.source_site())
                       || same_points(t.target_site(), r_.target_site());
    }

    if (is_tsrc_endp_of_p && is_tsrc_endp_of_q) {
      if (test_star(t.source_site(), p_, q_, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_p && is_ttrg_endp_of_q) {
      if (test_star(t.target_site(), p_, q_, t)) {
        return NEGATIVE;
      }
    }
    if (is_tsrc_endp_of_q && is_tsrc_endp_of_r) {
      if (test_star(t.source_site(), q_, r_, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_q && is_ttrg_endp_of_r) {
      if (test_star(t.target_site(), q_, r_, t)) {
        return NEGATIVE;
      }
    }
    if (is_tsrc_endp_of_r && is_tsrc_endp_of_p) {
      if (test_star(t.source_site(), r_, p_, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_r && is_ttrg_endp_of_p) {
      if (test_star(t.target_site(), r_, p_, t)) {
        return NEGATIVE;
      }
    }

    Line_2 l = compute_supporting_line(t.supporting_site());
    Sign sl = incircle(l, type);

    CGAL_SDG_DEBUG(std::cout
        << "debug vring incircle_s_no_easy: incircle l returned "
        << sl << std::endl;);

    if ( sl == POSITIVE ) { return sl; }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy sl=" << sl <<
      " d1=" << d1 << " d2=" << d2 << std::endl;);

    CGAL_SDG_DEBUG(std::cout << "debug vring numpts_in_pqr="
        << numpts_in_pqr << std::endl;);

    // philaris: here we have a serious change related to L2
    if ( sl == ZERO && (d1 == ZERO || d2 == ZERO) ) {

      // if some site in {p,q,r} is a point and it is also:
      // an endpoint of t and an endpoint of another site in {p,q,r}

      Site_2 sqpnt, other_t, other_seg;

      if (compute_helper(p_, q_, r_, t, sqpnt, other_t, other_seg)) {

        CGAL_assertion(sqpnt.is_point());
        CGAL_assertion(other_t.is_point());
        CGAL_assertion(other_seg.is_point());

        Point_2 vv (ux_, uy_, uz_);

        CGAL_SDG_DEBUG(std::cout << "debug vring incircle_s_no_easy "
            << "compute_helper true, "
            << "  vv=" << vv << "  sqpnt= " << sqpnt
            << "  other_t=" << other_t
            << "  other_seg=" << other_seg
            << std::endl;);

        Line_2 lvs =
          compute_line_from_to(vv, sqpnt.point());

        Oriented_side os_t =
          oriented_side_of_line(lvs, other_t.point());
        Oriented_side os_s =
          oriented_side_of_line(lvs, other_seg.point());

        CGAL_assertion(os_s != ON_ORIENTED_BOUNDARY);

        if (os_t == os_s) {
          Line_2 lseg =
            compute_line_from_to(sqpnt.point(), other_seg.point());

          Oriented_side os_seg_vv =
           oriented_side_of_line(lseg, vv);
          Oriented_side os_seg_t =
           oriented_side_of_line(lseg, other_t.point());

          if (os_seg_t == os_seg_vv) {
            return NEGATIVE;
          } else {
            if (os_seg_t == ON_ORIENTED_BOUNDARY) {
              return ZERO;
            } else {
              return POSITIVE;
            }
          }
        } // end of case: os_t == os_s
      } // end of case where

      return ZERO;
    }

    Oriented_side os1 = oriented_side_linf(l, t.source(), type);
    Oriented_side os2 = oriented_side_linf(l, t.target(), type);

    CGAL_SDG_DEBUG(std::cout << "debug vring incircle_s_no_easy: os1="
        << os1 << " os2="
        << os2 << std::endl;);

    if ( sl == ZERO ) {
      if (os1 == ON_ORIENTED_BOUNDARY || os2 == ON_ORIENTED_BOUNDARY) {
	return ZERO;
      }
      return ( os1 == os2 ) ? POSITIVE : ZERO;
    }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s_no_easy non-zero sl: os1="
      << os1 << " os2=" << os2 << std::endl;);

    return (os1 == os2) ? POSITIVE : NEGATIVE;
  }

  inline
  bool
  compute_helper(const Site_2& p, const Site_2& q, const Site_2& r,
      const Site_2& t,
      Site_2& sqpnt, Site_2& other_of_t, Site_2& other_of_seg)
  const
  {
    CGAL_assertion(t.is_segment());

    const bool is_p_point = p.is_point();
    const bool is_q_point = q.is_point();
    const bool is_r_point = r.is_point();

    const unsigned int numpts =
      ((is_p_point)? 1 : 0) +
      ((is_q_point)? 1 : 0) +
      ((is_r_point)? 1 : 0)  ;

    CGAL_SDG_DEBUG(std::cout << "debug vring compute_helper #pts="
        << numpts << std::endl;);

    if (numpts == 3) {
      return false;
    }

    // here and on, there are 1 or 2 points in {p,q,r}


    bool is_p_tsrc(false);
    bool is_p_endp_of_t(false);
    if (is_p_point) {
      is_p_tsrc = same_points(p, t.source_site());
      const bool is_p_ttrg = same_points(p, t.target_site());
      is_p_endp_of_t = is_p_tsrc || is_p_ttrg;

      if (is_p_endp_of_t) {
        sqpnt = p;
      }
    }

    bool is_q_tsrc(false);
    bool is_q_endp_of_t(false);
    if (is_q_point) {
      is_q_tsrc = same_points(q, t.source_site());
      const bool is_q_ttrg = same_points(q, t.target_site());
      is_q_endp_of_t = is_q_tsrc || is_q_ttrg;
      if (is_q_endp_of_t) {
        sqpnt = q;
      }
    }

    bool is_r_tsrc(false);
    bool is_r_endp_of_t(false);

    if (is_r_point) {
      is_r_tsrc = same_points(r, t.source_site());
      const bool is_r_ttrg = same_points(r, t.target_site());
      is_r_endp_of_t = is_r_tsrc || is_r_ttrg;
      if (is_r_endp_of_t) {
        sqpnt = r;
      }
    }

    const unsigned int numendpts_of_t =
      ((is_p_endp_of_t)? 1 : 0) +
      ((is_q_endp_of_t)? 1 : 0) +
      ((is_r_endp_of_t)? 1 : 0)  ;

    CGAL_SDG_DEBUG(std::cout << "debug compute_helper #endpts_of_t=" <<
      numendpts_of_t << std::endl;);

    if (numendpts_of_t == 0) {

      bool have_common_p_tsrc(false),
           have_common_p_ttrg(false),
           have_common_p_t(false);

      if (! is_p_point) {
        CGAL_assertion( ! same_segments(p, t) );
        const bool is_psrc_tsrc =
          same_points(p.source_site(), t.source_site());
        const bool is_ptrg_tsrc =
          same_points(p.target_site(), t.source_site());
        const bool is_psrc_ttrg =
          same_points(p.source_site(), t.target_site());
        const bool is_ptrg_ttrg =
          same_points(p.target_site(), t.target_site());
        have_common_p_tsrc = is_psrc_tsrc || is_ptrg_tsrc;
        have_common_p_ttrg = is_psrc_ttrg || is_ptrg_ttrg;
        have_common_p_t = have_common_p_tsrc || have_common_p_ttrg;
      }

      bool have_common_q_tsrc(false),
           have_common_q_ttrg(false),
           have_common_q_t(false);

      if (! is_q_point) {
        CGAL_assertion( ! same_segments(q, t) );
        const bool is_qsrc_tsrc = same_points(q.source_site(), t.source_site());
        const bool is_qtrg_tsrc = same_points(q.target_site(), t.source_site());
        const bool is_qsrc_ttrg = same_points(q.source_site(), t.target_site());
        const bool is_qtrg_ttrg = same_points(q.target_site(), t.target_site());
        have_common_q_tsrc = is_qsrc_tsrc || is_qtrg_tsrc;
        have_common_q_ttrg = is_qsrc_ttrg || is_qtrg_ttrg;
        have_common_q_t = have_common_q_tsrc || have_common_q_ttrg;
      }

      bool have_common_r_tsrc(false),
           have_common_r_ttrg(false),
           have_common_r_t(false);

      if (! is_r_point) {
        CGAL_assertion( ! same_segments(r, t) );
        const bool is_rsrc_tsrc = same_points(r.source_site(), t.source_site());
        const bool is_rtrg_tsrc = same_points(r.target_site(), t.source_site());
        const bool is_rsrc_ttrg = same_points(r.source_site(), t.target_site());
        const bool is_rtrg_ttrg = same_points(r.target_site(), t.target_site());
        have_common_r_tsrc = is_rsrc_tsrc || is_rtrg_tsrc;
        have_common_r_ttrg = is_rsrc_ttrg || is_rtrg_ttrg;
        have_common_r_t = have_common_r_tsrc || have_common_r_ttrg;
      }

      const unsigned int numcommon =
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

      const unsigned int numcommon_tsrc =
      ((have_common_p_tsrc)? 1 : 0) +
      ((have_common_q_tsrc)? 1 : 0) +
      ((have_common_r_tsrc)? 1 : 0)  ;

      const unsigned int numcommon_ttrg =
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

      if (have_common_p_t && have_common_q_t) {
        compute_helper_two_seg(p, q, sqpnt, other_of_seg);
      } else if (have_common_q_t && have_common_r_t) {
        compute_helper_two_seg(q, r, sqpnt, other_of_seg);
      } else if (have_common_r_t && have_common_p_t) {
        compute_helper_two_seg(r, p, sqpnt, other_of_seg);
      } else {
        CGAL_assertion(false);
      }

      return true;
    }

    // philaris: tocheck
    CGAL_assertion( numendpts_of_t == 1 );

    if (is_p_tsrc || is_q_tsrc || is_r_tsrc) {
      other_of_t = t.target_site();
    } else {
      other_of_t = t.source_site();
    }

    if (is_p_endp_of_t) {
      if (q.is_segment()) {
        bool is_p_qsrc = same_points(p, q.source_site());
        bool is_p_qtrg = same_points(p, q.target_site());
        if (is_p_qsrc || is_p_qtrg) {
          other_of_seg = is_p_qsrc ? q.target_site() : q.source_site();
          return true;
        }
      }
      if (r.is_segment()) {
        bool is_p_rsrc = same_points(p, r.source_site());
        bool is_p_rtrg = same_points(p, r.target_site());
        if (is_p_rsrc || is_p_rtrg) {
          other_of_seg = is_p_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

    } // end of case: is_p_endp_of_t

    if (is_q_endp_of_t) {
      if (r.is_segment()) {
        bool is_q_rsrc = same_points(q, r.source_site());
        bool is_q_rtrg = same_points(q, r.target_site());
        if (is_q_rsrc || is_q_rtrg) {
          other_of_seg = is_q_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

      if (p.is_segment()) {
        bool is_q_psrc = same_points(q, p.source_site());
        bool is_q_ptrg = same_points(q, p.target_site());
        if (is_q_psrc || is_q_ptrg) {
          other_of_seg = is_q_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

    } // end of case: is_q_endp_of_t

    if (is_r_endp_of_t) {
      if (p.is_segment()) {
        bool is_r_psrc = same_points(r, p.source_site());
        bool is_r_ptrg = same_points(r, p.target_site());
        if (is_r_psrc || is_r_ptrg) {
          other_of_seg = is_r_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

      if (q.is_segment()) {
        bool is_r_qsrc = same_points(r, q.source_site());
        bool is_r_qtrg = same_points(r, q.target_site());
        if (is_r_qsrc || is_r_qtrg) {
          other_of_seg = is_r_qsrc ? q.target_site() : q.source_site();
          return true;
        }
      }

    } // end of case: is_r_endp_of_t

    CGAL_SDG_DEBUG(std::cout << "debug compute_helper about to return false"
        << std::endl;);
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

    CGAL_SDG_DEBUG(std::cout
        << "debug compute_helper_two_seg entering with "
        << a << " and " << b << " having common "
        << common_site << std::endl;);

    if (is_site_h_or_v(a)) {
      if ( same_points(common_site, b.source_site()) ) {
        other_of_seg = b.target_site();
      } else {
        other_of_seg = b.source_site();
      }
    } else {
      CGAL_assertion(is_site_h_or_v(b));

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
      //CGAL_SDG_DEBUG(std::cout << "debug p_ref pps_idx="
      //    << pps_idx << std::endl;);

      if ( pps_idx == 0 ) {
        CGAL_assertion( p_.is_point());
        return p_;
      }

      if ( pps_idx == 1 ) {
        CGAL_assertion( q_.is_point());
        return q_;
      }

      //CGAL_SDG_DEBUG(std::cout << "debug p_ref about to return r="
      //    << r_ << std::endl;);

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
    CGAL_assertion( is_v_computed );
    return ux_;
  }

  FT hy() const {
    CGAL_assertion( is_v_computed );
    return uy_;
  }

  FT hw() const {
    CGAL_assertion( is_v_computed );
    return uz_;
  }

  FT radius() const {
    switch (v_type) {
    case PPP:    case PPS:    case PSS:
      {
	Point_2 pref = p_ref().point();
	//FT absdx = CGAL::abs(x() - pref.x());
	//FT absdy = CGAL::abs(y() - pref.y());
        return (CGAL::max)( CGAL::abs(x() - pref.x()),
		            CGAL::abs(y() - pref.y()) );
      }
      break;
    case SSS:
      {
	Line_2 l = compute_supporting_line(p_.supporting_site());
	Homogeneous_point_2 q = compute_linf_projection_hom(l, point());

	FT dx = CGAL::abs(x() - q.x());
	FT dy = CGAL::abs(y() - q.y());
	return (CGAL::max)(dx, dy);
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
    compute_v_if_not_computed();
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
    : p_(p), q_(q), r_(r), is_v_computed(false)
  {
    CGAL_SDG_DEBUG(std::cout << "Voronoi_vertex_ring_C2()" << std::endl;);
    analyze_vertex(p, q, r);
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
    CGAL_assertion( is_v_computed);
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
    compute_v_if_not_computed();
    Orientation o = orientation(l);

    if ( o == COLLINEAR ) { return ON_ORIENTED_BOUNDARY; }
    return ( o == LEFT_TURN ) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }


  // L_inf refinement
  inline
  Comparison_result
  linf_refine( const Line_2& l, Homogeneous_point_2& lrefhp ) const
  {
    CGAL_assertion( is_v_computed );
    Point_2 vv ( ux_, uy_, uz_ );

    bool is_l_h_or_v = is_line_h_or_v(l);

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    FT difxvl = vv.x() - lrefhp.x();
    FT difyvl = vv.y() - lrefhp.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);

    // lref is corner if and only if the line is not axis-parallel
    CGAL_assertion( (cmplabsxy == EQUAL) == (! is_l_h_or_v) );

    Oriented_side oslvv (ON_ORIENTED_BOUNDARY);
    if ((p_.is_segment() || q_.is_segment() || r_.is_segment()) &&
        is_l_h_or_v) {
      oslvv = oriented_side_of_line(l, vv);
      CGAL_assertion(oslvv != ON_ORIENTED_BOUNDARY);
    }

    // The following boolean variables are used to:
    // compute whether lref agrees with the x and y coordinates
    // of some of the points in p, q, r
    bool corner_agree_pt_x(false);
    bool corner_agree_pt_y(false);

    if (p_.is_point()) {
      Point_2 pp = p_.point();
      FT difxvp = vv.x() - pp.x();
      FT difyvp = vv.y() - pp.y();
      FT absdifxvp = CGAL::abs(difxvp);
      FT absdifyvp = CGAL::abs(difyvp);
      Comparison_result cmppabsxy = CGAL::compare(absdifxvp, absdifyvp);
      if (! ( (cmplabsxy == SMALLER) && (cmppabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvp) == EQUAL) {
          compare_p = CGAL::compare(absdifyvl, absdifyvp);
          corner_agree_pt_x = true;
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmppabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvp) == EQUAL) {
          CGAL_assertion(compare_p == EQUAL);
          compare_p = CGAL::compare(absdifxvl, absdifxvp);
          corner_agree_pt_y = true;
        }
      }
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslpsrc =
          oriented_side_of_line(l, p_.source_site().point());
        Oriented_side oslptrg =
          oriented_side_of_line(l, p_.target_site().point());
        if (((oslpsrc != oslvv) && (oslptrg != oslvv)) &&
            ((oslpsrc != ON_ORIENTED_BOUNDARY) ||
             (oslptrg != ON_ORIENTED_BOUNDARY)   )         ) {
          compare_p = SMALLER;
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
      if (! ( (cmplabsxy == SMALLER) && (cmpqabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvq) == EQUAL) {
          compare_q = CGAL::compare(absdifyvl, absdifyvq);
          corner_agree_pt_x = true;
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmpqabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvq) == EQUAL) {
          CGAL_assertion(compare_q == EQUAL);
          compare_q = CGAL::compare(absdifxvl, absdifxvq);
          corner_agree_pt_y = true;
        }
      }
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslqsrc =
          oriented_side_of_line(l, q_.source_site().point());
        Oriented_side oslqtrg =
          oriented_side_of_line(l, q_.target_site().point());
        if (((oslqsrc != oslvv) && (oslqtrg != oslvv)) &&
            ((oslqsrc != ON_ORIENTED_BOUNDARY) ||
             (oslqtrg != ON_ORIENTED_BOUNDARY)   )         ) {
          compare_q = SMALLER;
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
      if (! ( (cmplabsxy == SMALLER) && (cmprabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvr) == EQUAL) {
          compare_r = CGAL::compare(absdifyvl, absdifyvr);
          corner_agree_pt_x = true;
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmprabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvr) == EQUAL) {
          CGAL_assertion(compare_r == EQUAL);
          compare_r = CGAL::compare(absdifxvl, absdifxvr);
          corner_agree_pt_y = true;
        }
      }
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslrsrc =
          oriented_side_of_line(l, r_.source_site().point());
        Oriented_side oslrtrg =
          oriented_side_of_line(l, r_.target_site().point());
        if (((oslrsrc != oslvv) && (oslrtrg != oslvv)) &&
            ((oslrsrc != ON_ORIENTED_BOUNDARY) ||
             (oslrtrg != ON_ORIENTED_BOUNDARY)   )         ) {
          compare_r = SMALLER;
        }
      }
    }

    CGAL_SDG_DEBUG(std::cout << "debug linf_refine compare p q r = "
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) ||
        (compare_q == SMALLER) ||
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    if (corner_agree_pt_x && corner_agree_pt_y) {
      CGAL_assertion(! is_l_h_or_v);
      const unsigned int count_larger =
        (compare_p ? 1 : 0) +
        (compare_q ? 1 : 0) +
        (compare_r ? 1 : 0) ;
      if (count_larger >= 2) {
        // two points (among p, q, r) hide the line l from the vertex vv
        return LARGER;
      }
    }
    return EQUAL;

  }

  inline
  Comparison_result
  linf_refinement( Homogeneous_point_2& lrefhp ) const
  {
    CGAL_assertion( is_v_computed );
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
      if (! ( (cmplabsxy == SMALLER) && (cmppabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvp) == EQUAL) {
          compare_p = CGAL::compare(absdifyvl, absdifyvp);
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmppabsxy == LARGER ) ))
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
      if (! ( (cmplabsxy == SMALLER) && (cmpqabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvq) == EQUAL) {
          compare_q = CGAL::compare(absdifyvl, absdifyvq);
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmpqabsxy == LARGER ) ))
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
      if (! ( (cmplabsxy == SMALLER) && (cmprabsxy == SMALLER) ))
      {
        if (CGAL::compare(difxvl, difxvr) == EQUAL) {
          compare_r = CGAL::compare(absdifyvl, absdifyvr);
        }
      }
      if (! ( (cmplabsxy == LARGER ) && (cmprabsxy == LARGER ) ))
      {
        if (CGAL::compare(difyvl, difyvr) == EQUAL) {
          CGAL_assertion(compare_r == EQUAL);
          compare_r = CGAL::compare(absdifxvl, absdifxvr);
        }
      }
    }

    CGAL_SDG_DEBUG(std::cout << "debug compare p q r = "
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) ||
        (compare_q == SMALLER) ||
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    /*
    if ((compare_p == LARGER) ||
        (compare_q == LARGER) ||
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

  mutable bool is_v_computed;

  mutable RT ux_, uy_, uz_;
};


} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_RING_C2_H
