// Copyright (c) 2015  Università della Svizzera italiana.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H

#include <CGAL/license/Segment_Delaunay_graph_Linf_2.h>


#include <fstream>

#include <CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_x_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Compare_y_2.h>
#include <CGAL/Side_of_bounded_square_2.h>
#include <CGAL/Orientation_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>

namespace CGAL {

namespace SegmentDelaunayGraphLinf_2 {

template<class K>
class Voronoi_vertex_sqrt_field_new_C2
  : public Basic_predicates_C2<K>
{
public:
  typedef Basic_predicates_C2<K> Base;

  using Base::compute_supporting_line;
  using Base::oriented_side_of_line;
  using Base::opposite_line;
  using Base::compute_linf_projection_hom;
  using Base::compute_linf_projection_nonhom;
  using Base::compute_linf_perpendicular;
  using Base::compute_line_from_to;
  using Base::compute_horizontal_projection;
  using Base::compute_vertical_projection;
  using Base::has_positive_slope;
  using Base::have_same_slope;
  using Base::is_site_horizontal;
  using Base::is_site_vertical;
  using Base::is_site_h_or_v;
  using Base::is_line_h_or_v;
  using Base::test_star;
  using Base::compute_neg_45_line_at;
  using Base::compute_pos_45_line_at;
  using Base::compute_hor_line_at;
  using Base::compute_ver_line_at;
  using Base::are_in_same_open_halfspace_of;
  using Base::horseg_y_coord;
  using Base::verseg_x_coord;
  using Base::hvseg_coord;
  using Base::coord_at;
  using Base::touch_same_side;
  using Base::is_orth_dist_smaller_than_pt_dist;
  using Base::compute_intersection_of_lines;
  using Base::orient_lines_linf;
  using Base::are_parallel_lines;
  using Base::direction;
  using Base::compute_line_dir;
  using Base::parallel_bis;
  using Base::dir_from_lines;
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

  typedef enum {PPP = 0, PPS, PSS, SSS} vertex_t;
  struct PPP_Type {};
  struct PPS_Type {};
  struct PSS_Type {};
  struct SSS_Type {};

  typedef typename Base::Point_2             Point_2;
  typedef typename Base::Segment_2           Segment_2;
  typedef typename Base::Line_2              Line_2;
  typedef typename Base::Site_2              Site_2;
  typedef typename Base::Direction_2         Direction_2;
  typedef typename Base::FT                  FT;
  typedef typename Base::RT                  RT;

  typedef typename Base::Homogeneous_point_2 Homogeneous_point_2;

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Bounded_side        Bounded_side;
  typedef typename Base::Sign                Sign;
  typedef typename Base::Compute_scalar_product_2 Compute_scalar_product_2;

  typedef typename Base::Polychainline_2     Polychainline_2;

  typedef typename Base::Bearing Bearing;

private:
  typedef SegmentDelaunayGraph_2::Are_same_points_C2<K>   Are_same_points_2;
  typedef SegmentDelaunayGraph_2::Are_same_segments_C2<K> Are_same_segments_2;
  typedef Side_of_bounded_square_2<K>    Side_of_bounded_square_2_Type;
  typedef Orientation_Linf_2<K>          Orientation_Linf_points_2;
  typedef Bisector_Linf<K>               Bisector_Linf_Type;

  typedef SegmentDelaunayGraph_2::Compare_x_2<K> Compare_x_2_Sites_Type;
  typedef SegmentDelaunayGraph_2::Compare_y_2<K> Compare_y_2_Sites_Type;

  typedef typename K::Compare_x_2 Compare_x_2_Points_Type;
  typedef typename K::Compare_y_2 Compare_y_2_Points_Type;

  Are_same_points_2     same_points;
  Are_same_segments_2   same_segments;
  Side_of_bounded_square_2_Type    side_of_bounded_square;
  Compare_x_2_Sites_Type           scmpx;
  Compare_y_2_Sites_Type           scmpy;
  Compare_x_2_Points_Type          cmpx;
  Compare_y_2_Points_Type          cmpy;
  Orientation_Linf_points_2 or_linf;
  Bisector_Linf_Type bisector_linf;
  Bisector_Linf_Type linf_bisect_direction;

private:
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Voronoi vertex computation
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  // the Voronoi vertex of three points
  //--------------------------------------------------------------------------

  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
             const PPP_Type&) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
                       sr.is_point() );

    //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew entering compute_vv" << std::endl;);

    // the following check is not really needed in this
    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    Point_2 p = sp.point(), q = sq.point(), r = sr.point();

    return compute_vv_points(p, q, r);
  }

  inline
  void
  compute_vv_points(
      const Point_2 & p, const Point_2 & q, const Point_2 & r)
  const
  {
    //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew (p q r) = " <<
    //  p << ' ' << q << ' ' << r << std::endl;);

    FT x_min, x_max, y_min, y_max;
    FT x_center, y_center;
    FT half(0.5);
    FT two(2);

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

      //CGAL_SDG_DEBUG(std::cout << "debug set y_center=" <<
      //  y_center << std::endl;);

      Comparison_result cmpxrothers = CGAL::compare(r.x(), p.x());
      if (cmpxrothers == SMALLER) {
        //CGAL_SDG_DEBUG(std::cout << "debug r is left of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  && (cmpyrq == LARGER)) ||
            ((cmpyrp == SMALLER) && (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two*y_center - r.y();
            is_set_y_min = true;
            //CGAL_SDG_DEBUG(std::cout << "debug set y_min=" <<
            //  y_min << std::endl;);
          } else {
            y_max = two*y_center - r.y();
            is_set_y_max = true;
            //CGAL_SDG_DEBUG(std::cout << "debug set y_max=" <<
            //  y_max << std::endl;);
          }
        }
      } else if (cmpxrothers == LARGER) {
        //CGAL_SDG_DEBUG(std::cout << "debug r is right of p, q" << std::endl;);
        Comparison_result cmpyrp = CGAL::compare(r.y(), p.y());
        Comparison_result cmpyrq = CGAL::compare(r.y(), q.y());
        if (((cmpyrp == LARGER)  && (cmpyrq == LARGER)) ||
            ((cmpyrp == SMALLER) && (cmpyrq == SMALLER))
           ) {
          // do fix
          if (cmpyrp == LARGER) {
            y_min = two*y_center - r.y();
            is_set_y_min = true;
            //CGAL_SDG_DEBUG(std::cout << "debug set y_min=" <<
            //  y_min << std::endl;);
          } else {
            y_max = two*y_center - r.y();
            is_set_y_max = true;
            //CGAL_SDG_DEBUG(std::cout << "debug set y_max=" <<
            //  y_max << std::endl;);
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
      x_center = half * (p.x() + q.x());
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
        if (((cmpxrp == LARGER)  && (cmpxrq == LARGER)) ||
            ((cmpxrp == SMALLER) && (cmpxrq == SMALLER))
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
          y_center = half * (p.y() + r.y());
          //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
          Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
          if ((cmpyqp == LARGER) && (cmpyqr == LARGER)) {
            y_min = two*y_center - q.y();
            is_set_y_min = true;
          }
          if ((cmpyqp == SMALLER) && (cmpyqr == SMALLER)) {
            y_max = two*y_center - q.y();
            is_set_y_max = true;
          }
        } else {
          y_center = half * (q.y() + r.y());
          Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
          Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
          if ((cmpypq == LARGER) && (cmpypr == LARGER)) {
            y_min = two*y_center - p.y();
            is_set_y_min = true;
          }
          if ((cmpypq == SMALLER) && (cmpypr == SMALLER)) {
            y_max = two*y_center - p.y();
            is_set_y_max = true;
          }
        }
        is_set_y_center = true;
      }
    } else {
      // here r.x() = x_min
      // r.x() = p.x() || r.x() = q.x()
      if (CGAL::compare(r.x(), p.x()) == EQUAL) {
        //CGAL_SDG_DEBUG(std::cout << "debug r.x = p.x" << std::endl;);
        // r.x() = p.x()
        y_center = half * (p.y() + r.y());
        //Comparison_result cmpyqp = CGAL::compare(q.y(),p.y());
        Comparison_result cmpyqr = CGAL::compare(q.y(),r.y());
        if ((cmpyqp == LARGER) && (cmpyqr == LARGER)) {
          //CGAL_SDG_DEBUG(std::cout << "debug q is above p, r" << std::endl;);
          y_min = two*y_center - q.y();
          is_set_y_min = true;
        }
        if ((cmpyqp == SMALLER) && (cmpyqr == SMALLER)) {
          //CGAL_SDG_DEBUG(std::cout << "debug q is below p, r" << std::endl;);
          y_max = two*y_center - q.y();
          is_set_y_max = true;
        }
      } else {
        // r.x() = q.x()
        //CGAL_SDG_DEBUG(std::cout << "debug r.x = q.x" << std::endl;);
        y_center = half * (q.y() + r.y());
        Comparison_result cmpypq = CGAL::compare(p.y(),q.y());
        Comparison_result cmpypr = CGAL::compare(p.y(),r.y());
        if ((cmpypq == LARGER) && (cmpypr == LARGER)) {
          //CGAL_SDG_DEBUG(std::cout << "debug p is above q, r" << std::endl;);
          y_min = two*y_center - p.y();
          is_set_y_min = true;
        }
        if ((cmpypq == SMALLER) && (cmpypr == SMALLER)) {
          //CGAL_SDG_DEBUG(std::cout << "debug p is below q, r" << std::endl;);
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
          x_center = half * (p.x() + r.x());
          //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
          Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
          if ((cmpxqp == LARGER) && (cmpxqr == LARGER)) {
            x_min = two*x_center - q.x();
            is_set_x_min = true;
          }
          if ((cmpxqp == SMALLER) && (cmpxqr == SMALLER)) {
            x_max = two*x_center - q.x();
            is_set_x_max = true;
          }
        } else {
          x_center = half * (q.x() + r.x());
          Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
          Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
          if ((cmpxpq == LARGER) && (cmpxpr == LARGER)) {
            x_min = two*x_center - p.x();
            is_set_x_min = true;
          }
          if ((cmpxpq == SMALLER) && (cmpxpr == SMALLER)) {
            x_max = two*x_center - p.x();
            is_set_x_max = true;
          }
        }
        is_set_x_center = true;
      }
    } else {
      // here r.y() = y_min
      // r.y() = p.y() || r.y() = q.y()
      if (CGAL::compare(r.y(), p.y()) == EQUAL) {
        x_center = half * (p.x() + r.x());
        //Comparison_result cmpxqp = CGAL::compare(q.x(),p.x());
        Comparison_result cmpxqr = CGAL::compare(q.x(),r.x());
        if ((cmpxqp == LARGER) && (cmpxqr == LARGER)) {
          x_min = two*x_center - q.x();
          is_set_x_min = true;
        }
        if ((cmpxqp == SMALLER) && (cmpxqr == SMALLER)) {
          x_max = two*x_center - q.x();
          is_set_x_max = true;
        }
      } else {
        x_center = half * (q.x() + r.x());
        Comparison_result cmpxpq = CGAL::compare(p.x(),q.x());
        Comparison_result cmpxpr = CGAL::compare(p.x(),r.x());
        if ((cmpxpq == LARGER) && (cmpxpr == LARGER)) {
          x_min = two*x_center - p.x();
          is_set_x_min = true;
        }
        if ((cmpxpq == SMALLER) && (cmpxpr == SMALLER)) {
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
        //CGAL_SDG_DEBUG(std::cout << "rectangle has to be made fatter" << std::endl;);
        // make rectangle fatter
        if (is_set_x_center) {
          //CGAL_SDG_DEBUG(std::cout << "x_center already set" << std::endl;);
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
          //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew grow right" << std::endl;);
          x_max = x_min + y_max - y_min;
        } else
        { // grow rectangle to the left
          //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew grow left" << std::endl;);
          x_min = x_max - y_max + y_min;
        }
        break;
      case LARGER:
        //CGAL_SDG_DEBUG(std::cout << "debug rectangle has to be made taller" << std::endl;);
        // make rectangle taller
        if (is_set_y_center) {
          // grow in both sides
          //CGAL_SDG_DEBUG(std::cout << "debug y_center already set" << std::endl;);
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
          //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew grow upwards" << std::endl;);
          y_max = y_min + x_max - x_min;
        } else
        { // grow rectangle downwards
          //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew grow downwards" << std::endl;);
          y_min = y_max - x_max + x_min;
        }
        break;
      case EQUAL:
        // do nothing
        break;
    }

    FT ux, uy, uz;

    ux = x_min + x_max;
    uy = y_min + y_max;
    uz = FT(2) ;

    //CGAL_SDG_DEBUG(std::cout << "debug vsqrnew sets vv = "
    //          << ux/uz << ' ' << uy/uz << std::endl;);

    vv = Point_2(ux / uz, uy / uz);
  }


  //--------------------------------------------------------------------------
  // the Voronoi vertex of two points and a segment
  //--------------------------------------------------------------------------

  inline void
  compute_pps_endp_hv(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r, const bool is_r_horizontal) const
  {
    CGAL_USE(r);
    const Site_2 & A = p_endp_r ? p : q;
    const Site_2 & B = p_endp_r ? q : p;
    const FT Apar = is_r_horizontal ? A.point().x() : A.point().y();
    const FT Aort = is_r_horizontal ? A.point().y() : A.point().x();
    const FT Bpar = is_r_horizontal ? B.point().x() : B.point().y();
    const FT Bort = is_r_horizontal ? B.point().y() : B.point().x();
    const FT dpar = Apar - Bpar;
    const FT dort = Aort - Bort;
    const FT absdpar = CGAL::abs(dpar);
    FT vx_, vy_;
    FT & vpar = is_r_horizontal ? vx_ : vy_;
    FT & vort = is_r_horizontal ? vy_ : vx_;

    if (2*absdpar < CGAL::abs(dort)) {
      vpar = Apar;
      vort = Aort - (dort)/(FT(2));
    } else {
      vpar = Apar;
      vort = Aort - CGAL::sign(dort)*absdpar;
    }
    vv = Point_2(vx_, vy_);
  }

  inline void
  compute_pps_endp_slope(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r, const bool pos_slope) const
  {
    CGAL_USE(r);
    const Site_2 & A = p_endp_r ? p : q;
    const Site_2 & B = p_endp_r ? q : p;
    const FT Ax = A.point().x();
    const FT Ay = A.point().y();
    const FT Bx = B.point().x();
    const FT By = B.point().y();
    const FT dx = Ax - Bx;
    const FT dy = Ay - By;
    const FT absdx = CGAL::abs(dx);
    const FT absdy = CGAL::abs(dy);
    FT x_, y_;
    if (absdx > absdy) {
      x_ = FT(2)*Ax - dx;
      y_ = FT(2)*Ay - FT(pos_slope? -1 : +1)*dx;
    } else {
      x_ = FT(2)*Ax - FT(pos_slope? -1 : +1)*dy;
      y_ = FT(2)*Ay - dy;
    }
    vv = Point_2(x_/FT(2), y_/FT(2));
  }

  inline void
  compute_pps_endp(const Site_2& p, const Site_2& q, const Site_2& r,
                   const bool p_endp_r) const
  {
    const bool is_r_horizontal = is_site_horizontal(r);
    if (is_r_horizontal || is_site_vertical(r)) {
      return compute_pps_endp_hv(p, q, r, p_endp_r, is_r_horizontal);
    } else {
      const bool pos_slope = has_positive_slope(r);
      return compute_pps_endp_slope(p, q, r, p_endp_r, pos_slope);
    }
  }

  inline void
  compute_pps_nonendp_hv_samecoord(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_r_horizontal) const
  {
    const FT ppar = is_r_horizontal ? p.point().x() : p.point().y();
    const FT port = is_r_horizontal ? p.point().y() : p.point().x();
    const FT qort = is_r_horizontal ? q.point().y() : q.point().x();
    FT vx_, vy_;
    FT & vpar = is_r_horizontal ? vx_ : vy_;
    FT & vort = is_r_horizontal ? vy_ : vx_;
    const FT segort = (is_r_horizontal)?
      horseg_y_coord(r) : verseg_x_coord(r);
    const FT sumort = port + qort;
    vort = sumort/FT(2);
    const int vhsign = is_r_horizontal ? +1 : -1;
    const int distsign = CGAL::abs(segort-qort) < CGAL::abs(segort-port) ?
      +1: -1;
    vpar = ppar - FT(vhsign*distsign)*(segort-vort);
    vv = Point_2(vx_, vy_);
    CGAL_SDG_DEBUG(std::cout << "debug: PPS returns with vv=" << vv << std::endl;);
  }

  /* compute pps vertex when the points p, q are not endpoints of
   * segment r, and when segment r is axis-parallel
   */
  inline void
  compute_pps_nonendp_hv(const Site_2& p, const Site_2& q, const Site_2& r,
                         const bool is_r_horizontal) const
  {
    if ((is_r_horizontal       && (scmpx(p, q) == EQUAL)) ||
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

    if (perpcomp == EQUAL) {
      const RT pqdist = is_r_horizontal ?
        CGAL::abs(pp.x()-qq.x()) : CGAL::abs(pp.y()-qq.y());
      const RT signrdist = (is_r_horizontal ? pp.y() : pp.x()) - coordr;
      Comparison_result comp = CGAL::compare(pqdist, CGAL::abs(signrdist));
      vv = Point_2(is_r_horizontal ?
                     (pp.x() + qq.x()) :
                     (comp == LARGER) ?
                       RT(2)*coordr + CGAL::sign(signrdist)*pqdist :
                       coordr + pp.x(),
                   is_r_horizontal ?
                     (comp == LARGER) ?
                       RT(2)*coordr + CGAL::sign(signrdist)*pqdist :
                       coordr + pp.y() :
                     (pp.y() + qq.y()),
                   RT(2));
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
      vv = Point_2(is_r_horizontal ?
                   (RT(2) * closep.x() - (is_p_farthest? -1 : +1) * sdistf) :
                   (RT(2)*coordr + sdistf),
                   is_r_horizontal ?
                   (RT(2)*coordr + sdistf) :
                   (RT(2) * closep.y() + (is_p_farthest? -1 : +1) * sdistf),
                   RT(2));
    } else {
      vv = Point_2(is_r_horizontal ?
                   (pp.x() + qq.x()) :
                   (RT(2)*coordr + CGAL::sign(sdistf)*pqdist),
                   is_r_horizontal ?
                   (RT(2)*coordr + CGAL::sign(sdistf)*pqdist) :
                   (pp.y() + qq.y()),
                   RT(2));
    }
    return;
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
      const FT pcoord = pos_slope ? p.point().x() : p.point().y();
      const FT lineval = coord_at(l, pcoord, pos_slope);
      const Point_2 corner = pos_slope?
        Point_2(pcoord, lineval) : Point_2(lineval, pcoord);
      const FT sidelen = (CGAL::max)(CGAL::abs(corner.x() - q.point().x()),
                                     CGAL::abs(corner.y() - q.point().y()));
      vv = Point_2(FT(2)*corner.x() + signla*sidelen,
                   FT(2)*corner.y() + signlb*sidelen,
                   FT(2));
      return;
    }
    if (second_comp == second_value) {
      const FT qcoord = pos_slope ? q.point().y() : q.point().x();
      const FT lineval = coord_at(l, qcoord, ! pos_slope);
      const Point_2 corner = pos_slope?
        Point_2(lineval, qcoord) : Point_2(qcoord, lineval);
      const FT sidelen = (CGAL::max)(CGAL::abs(corner.x() - p.point().x()),
                                     CGAL::abs(corner.y() - p.point().y()));
      vv = Point_2(FT(2)*corner.x() + signla*sidelen,
                   FT(2)*corner.y() + signlb*sidelen,
                   FT(2));
      return;
    }

    CGAL_assertion((first_comp  == -first_value ) &&
                   (second_comp == -second_value)    );

    const FT px = p.point().x();
    const FT py = p.point().y();
    const FT qx = q.point().x();
    const FT qy = q.point().y();
    const FT pqdist = (CGAL::max)(CGAL::abs(px - qx), CGAL::abs(py - qy));

    CGAL_SDG_DEBUG(std::cout
        << "debug: vsqrt pqdist=" << pqdist << std::endl;);

    const FT & pcoord = pos_slope ? px : py;
    const FT plineval = coord_at(l, pcoord, pos_slope);
    const FT & pothercoord = pos_slope ? py : px;
    const FT plen = CGAL::abs(plineval -  pothercoord);
    CGAL_SDG_DEBUG(std::cout
        << "debug: vsqrt plen=" << plen << std::endl;);
    if (CGAL::compare(pqdist, plen) != SMALLER) {
      // here, appropriate projection of p on supporting line of segment r
      // is shorter than Linf p, q distance
      const Point_2 corner = pos_slope?
        Point_2(pcoord, plineval) : Point_2(plineval, pcoord);
      vv = Point_2(FT(2)*corner.x() + signla*pqdist,
                   FT(2)*corner.y() + signlb*pqdist,
                   FT(2));
      return;
    }

    const FT & qcoord = pos_slope ? qy : qx;
    const FT qlineval = coord_at(l, qcoord, ! pos_slope);
    const FT & qothercoord = pos_slope ? qx : qy;
    const FT qlen = CGAL::abs(qlineval -  qothercoord);
    CGAL_SDG_DEBUG(std::cout
        << "debug: vsqrt qlen=" << qlen << std::endl;);
    if (CGAL::compare(pqdist, qlen) != SMALLER) {
      // here, appropriate projection of q on supporting line of segment r
      // is shorter than Linf p, q distance
      const Point_2 corner = pos_slope?
        Point_2(qlineval, qcoord) : Point_2(qcoord, qlineval);
      vv = Point_2(FT(2)*corner.x() + signla*pqdist,
                   FT(2)*corner.y() + signlb*pqdist,
                   FT(2));
      return;
    }

    CGAL_assertion((pqdist < plen) && (pqdist < qlen));

    // here, compute corner opposite of corner on line of segment r
    const Point_2 opposite_corner = pos_slope ?
      Point_2(qx, py) : Point_2(px, qy);
    CGAL_SDG_DEBUG(std::cout << "debug: vsqrt opposite_corner="
        << opposite_corner << std::endl;);

    const Point_2 corner =
      compute_linf_projection_nonhom(l, opposite_corner);

    vv = Point_2(corner.x() + opposite_corner.x(),
                 corner.y() + opposite_corner.y(),
                 FT(2));
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
      vv = (samexpq) ? Point_2(vsamecoord, vdiffcoord) :
                       Point_2(vdiffcoord, vsamecoord) ;
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
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
             const PPS_Type& type) const
  {
    CGAL_precondition( sp.is_point() && sq.is_point() &&
                       sr.is_segment() );
    CGAL_USE(type);

    CGAL_SDG_DEBUG(std::cout
        << "debug: compute_vv PPS entering p=" << sp
        << " q=" << sq << " r=" << sr << std::endl;);

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    const bool p_endp_r = is_endpoint_of(sp, sr);
    const bool q_endp_r = is_endpoint_of(sq, sr);

    if (p_endp_r || q_endp_r) {
      return compute_pps_endp(sp, sq, sr, p_endp_r);
    }
    CGAL_assertion(are_in_same_open_halfspace_of(sp, sq, sr));
    return compute_pps_nonendp(sp, sq, sr);
  }

  inline
  void compute_vv_bisectors(
      const Site_2& sp, const Site_2& sq, const Site_2& sr,
      const PPS_Type&) const
  {
    Polychainline_2 bpq = bisector_linf(sp, sq);
    CGAL_SDG_DEBUG(std::cout << "debug: bpq p=" << sp << " q=" << sq << std::endl;);
    CGAL_SDG_DEBUG(std::cout << "debug: bpq =" << bpq << std::endl;);

    bool samexpq = (scmpx(sp, sq) == EQUAL);
    bool sameypq = (scmpy(sp, sq) == EQUAL);

    CGAL_assertion(! (samexpq && sameypq));

    bool samecoordpq = samexpq || sameypq ;

    Polychainline_2 goodbisector;
    if (is_endpoint_of(sp, sr)) {
      goodbisector = bisector_linf(sr, sp);
      CGAL_SDG_DEBUG(std::cout << "debug: brp r=" << sr << " p=" << sp << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: brp res=" << goodbisector << std::endl;);
    } else if (is_endpoint_of(sq, sr)) {
      goodbisector = bisector_linf(sq, sr);
      CGAL_SDG_DEBUG(std::cout << "debug: bqr q=" << sq << " r=" << sr << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: bqr res=" << goodbisector << std::endl;);
    } else if (samecoordpq) {
      CGAL_SDG_DEBUG(std::cout << "debug PPS samecoordpq" << std::endl;);

      // check which of points p, q is closer to segment r

      bool use_bqr;

      Line_2 l = compute_supporting_line(sr);

      if (((CGAL::sign(l.a()) == ZERO) && sameypq) ||
          ((CGAL::sign(l.b()) == ZERO) && samexpq)   )  {
        // here l is horizontal or vertical and parallel to pq;
        // bqr or brp are equally good
        use_bqr = true;
      } else {
        // here l and segment are neither hor. nor ver.
        Point_2 proj;
        FT projft, pft, qft;
        if (samexpq) {
          // compute vertical projection
          proj = compute_vertical_projection(l, sp.point());
          projft = proj.y();
          pft = sp.point().y();
          qft = sq.point().y();


        } else {
          CGAL_assertion(sameypq);
          // compute horizontal projection
          proj = compute_horizontal_projection(l, sp.point());
          projft = proj.x();
          pft = sp.point().x();
          qft = sq.point().x();
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
        goodbisector = bisector_linf(sq, sr);
        CGAL_SDG_DEBUG(std::cout << "debug: bqr q=" << sq << " r=" << sr << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug: bqr res=" << goodbisector << std::endl;);
      } else {
        goodbisector = bisector_linf(sr, sp);
        CGAL_SDG_DEBUG(std::cout << "debug: brp r=" << sr << " p=" << sp << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug: brp res=" << goodbisector << std::endl;);
      }
    } else {
      goodbisector = bisector_linf(sq, sr);
      CGAL_SDG_DEBUG(std::cout << "debug: bqr q=" << sq << " r=" << sr << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: bqr res=" << goodbisector << std::endl;);
    }

    vv = bpq.first_intersection_point_with(goodbisector);
    CGAL_SDG_DEBUG(std::cout << "debug: PPS returns with vv=" << vv << std::endl;);
  }


  //--------------------------------------------------------------------------
  // the Voronoi vertex of a point and two segments
  // also: the Voronoi vertex of a point and two lines
  //--------------------------------------------------------------------------


  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
             const PSS_Type&) const
  {
#ifdef CGAL_PROFILE
    // In case CGAL profile is called then output the sites in case of
    // a filter failure
    if ( Algebraic_structure_traits<FT>::Is_exact::value ) {
      std::ofstream ofs("vv-failure-log.cin", std::ios_base::app);
      ofs.precision(16);
      ofs << sp << std::endl;
      ofs << sq << std::endl;
      ofs << sr << std::endl;
      ofs << "=======" << std::endl;
      ofs.close();
    }
#endif

    CGAL_precondition( sp.is_point() && sq.is_segment() &&
                       sr.is_segment() );

    CGAL_SDG_DEBUG(std::cout
        << "debug: computevv PSS entering p=" << sp
        << " q=" << sq << " r=" << sr << std::endl;);

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    const bool pq = is_endpoint_of(sp, sq);
    const bool pr = is_endpoint_of(sp, sr);

    Point_2 pp = sp.point();

    if ( pq && pr ) {
      CGAL_SDG_DEBUG(std::cout
          << "debug field_new setting vv = pp = "
          << pp << std::endl;);
      vv = pp;
      return;
    }
    const bool is_q_hor = is_site_horizontal(sq);
    const bool is_q_ver = is_site_vertical(sq);
    const bool is_r_hor = is_site_horizontal(sr);
    const bool is_r_ver = is_site_vertical(sr);
    const bool is_q_hv = is_q_hor || is_q_ver;
    const bool is_r_hv = is_r_hor || is_r_ver;
    if (is_q_hv && is_r_hv) {
      compute_pss_both_hv(sp, sq, sr, is_q_hor, is_r_hor, pq, pr);
    } else {
      if (pq || pr) {
        compute_pss_endp(sp, sq, sr,
            is_q_hv, is_q_hor, pq, is_r_hv, is_r_hor, pr);
      } else {
        compute_pss_nonendp(sp, sq, sr,
            is_q_hv, is_q_hor, is_r_hv, is_r_hor);
      }
    }
  }

  // PSS case when not both segments are axis-parallel and p is
  // not an endpoint of any of the segments q and r
  inline void
  compute_pss_nonendp(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hv, const bool is_q_hor,
      const bool is_r_hv, const bool is_r_hor) const
  {
    CGAL_USE(is_q_hv);
    CGAL_USE(is_q_hor);
    CGAL_USE(is_r_hv);
    CGAL_USE(is_r_hor);
    const Line_2 lq = orient_line_nonendp(p, q);
    const Line_2 lr = orient_line_nonendp(p, r);
    const Bearing bq = bearing(lq);
    const Bearing br = bearing(lr);
    const Bearing bdiff = bearing_diff(bq, br);
    CGAL_assertion( bdiff != 0 );
    CGAL_assertion( bdiff != 7 );

    if (bdiff == 1) {
      compute_pss_corner_and_pt(p, lq, lr, bq, br);
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
    CGAL_assertion( oriented_side_of_line(lq, this->point()) != ZERO );
    CGAL_assertion( oriented_side_of_line(lr, this->point()) != ZERO );
  }

  inline void
  compute_pss_lines_side(const Site_2& p,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bside) const
  {
    CGAL_precondition(bside % 2 == 1);
    const bool side_ver = (bside % 4 == 1);
    FT pcoord = (side_ver) ? p.point().x() : p.point().y();
    FT qcoord = coord_at(lq, pcoord, side_ver);
    FT rcoord = coord_at(lr, pcoord, side_ver);
    FT sidelen = CGAL::abs(qcoord-rcoord);
    const int sgn = (bside < 4) ? -1 : +1;
    vv = side_ver ?
      Point_2(pcoord + sgn*sidelen/FT(2), (qcoord+rcoord)/FT(2)) :
      Point_2((qcoord+rcoord)/FT(2), pcoord + sgn*sidelen/FT(2)) ;
  }

  inline void
  compute_pss_side_p_known(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br) const
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
    vv = center_from_same_side_corners(rcorner, qcorner, bside);
  }

  inline void
  compute_pss_ortho_wedge(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br) const
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
    vv = center_from_opposite_corners(Point_2(hx, hy, hw), corner);
  }

  inline void
  compute_pss_nonhv_consecutive(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br) const
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
      const Bearing bqr) const
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
      vv = Point_2(RT(2)*xs + (yp - ys), yp + ys, RT(2));
    } else if (CGAL::compare(yp, yr) == ((bqr == 1) ? LARGER : SMALLER)) {
      // p close to r
      const FT xs = coord_at(lr, yp, false);
      const FT ys = coord_at(lq, xs, true);
      vv = Point_2(RT(2)*xs + (ys - yp), yp + ys, RT(2));
    } else {
      // p on opposite side of two lines (or on its corners)
      vv = Point_2(xp + x, yq + yr, RT(2));
    }
  }

  inline void
  compute_pss_y_consecutive(
      const Site_2& p, const Site_2& q, const Site_2& r,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br,
      const Bearing bqr) const
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
      vv = Point_2(xp + xs, RT(2)*ys + (xs - xp), RT(2));
    } else if (CGAL::compare(xp, xr) == ((bqr == 3) ? SMALLER : LARGER)) {
      // p close to r
      const FT ys = coord_at(lr, xp, true);
      const FT xs = coord_at(lq, ys, false);
      vv = Point_2(xp + xs, RT(2)*ys + (xp - xs), RT(2));
    } else {
      // p on opposite side of two lines (or on its corners)
      vv = Point_2(xq + xr, yp + y, RT(2));
    }
  }

  inline void
  compute_pss_corner_and_pt(const Site_2& p,
      const Line_2& lq, const Line_2 & lr,
      const Bearing bq, const Bearing br) const
  {
    RT cx, cy, cw;
    compute_intersection_of_lines(lq, lr, cx, cy, cw);
    const Bearing cb = (bq % 2 == 0) ? bq : br;
    vv = center_from_corner_and_pt(Point_2(cx, cy, cw), cb, p.point());
  }

  // PSS case when not both segments are axis-parallel and p is
  // an endpoint of one of the segments
  inline void
  compute_pss_endp(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hv, const bool is_q_hor, const bool pq,
      const bool is_r_hv, const bool is_r_hor, const bool pr) const
  {
    CGAL_precondition(pq || pr);
    CGAL_USE(pr);
    const Line_2 lendp = orient_line_endp(p, (pq ? q : r), pq);
    const Line_2 lnon = orient_line_nonendp(p, (pq ? r : q));
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
    RT ux, uy, uz;
    compute_intersection_of_lines(llbis, lperp, ux, uy, uz);
    vv = Point_2(ux, uy, uz);
    CGAL_assertion( oriented_side_of_line(lendp, this->point()) != ZERO );
    CGAL_assertion( oriented_side_of_line(lnon, this->point()) != ZERO );
  }

  // both segments are axis-parallel
  inline void
  compute_pss_both_hv(const Site_2& p, const Site_2& q, const Site_2& r,
      const bool is_q_hor, const bool is_r_hor,
      const bool pq, const bool pr) const
  {
    CGAL_precondition(! (pq && pr));
    if (is_q_hor == is_r_hor) {
      // parallel segments
      const RT q_coord = hvseg_coord(q, is_q_hor);
      const RT r_coord = hvseg_coord(r, is_r_hor);
      RT ux_, uy_, uz_;
      RT & upar = is_q_hor ? ux_ : uy_;
      RT & uort = is_q_hor ? uy_ : ux_;
      upar = RT(2)*(is_q_hor ? p.point().x() : p.point().y())
        + (( pq || pr ) ? RT(0) :
                          RT(is_q_hor ? +1 : -1)*(r_coord - q_coord));
      uort = q_coord + r_coord;
      uz_ = RT(2);
      vv = Point_2(ux_, uy_, uz_);
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
      const bool pq, const bool pr) const
  {
    CGAL_precondition(is_q_hor != is_r_hor);
    if (pq || pr) {
      const RT q_coord = hvseg_coord(q, is_q_hor);
      const RT r_coord = hvseg_coord(r, is_r_hor);
      const bool is_touched_hor = pq ? is_q_hor : is_r_hor;
      const RT coord_c = is_touched_hor ? p.point().x() : p.point().y();
      const RT radius = CGAL::abs(coord_c - (pq ? r_coord : q_coord));
      RT ux_, uy_, uz_;
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
      vv = Point_2(ux_, uy_, uz_);
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
      const bool pq, const bool pr) const
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
    RT ux_, uy_, uz_;
    if (cmp == LARGER) {
      ux_ = is_q_hor ? r_coord + p_coord_r : q_coord + p_coord_q;
      uy_ = is_q_hor ? (RT(2)*q_coord + CGAL::sign(sdistq)*dx) :
                       (RT(2)*r_coord + CGAL::sign(sdistr)*dx) ;
    } else if (cmp == SMALLER) {
      uy_ = is_r_hor ? r_coord + p_coord_r : q_coord + p_coord_q;
      ux_ = is_r_hor ? (RT(2)*q_coord + CGAL::sign(sdistq)*dy) :
                       (RT(2)*r_coord + CGAL::sign(sdistr)*dy) ;
    } else {
      ux_ = is_q_hor ? r_coord + p_coord_r : q_coord + p_coord_q;
      uy_ = is_q_hor ? q_coord + p_coord_q : r_coord + p_coord_r;
    }
    uz_ = RT(2);
    vv = Point_2(ux_, uy_, uz_);
  }

  inline void
  compute_vv_bisectors(
      const Site_2& sp, const Site_2& sq, const Site_2& sr,
      const PSS_Type&) const
  {
    const bool pq = is_endpoint_of(sp, sq);
    Polychainline_2 goodbisector;
    if (pq) {
      goodbisector = bisector_linf(sp, sq);
      CGAL_SDG_DEBUG(std::cout
          << "debug: computevv bpq p=" << sp << " q=" << sq << std::endl;);
      CGAL_SDG_DEBUG(std::cout
          << "debug: computevv bpq =" << goodbisector << std::endl;);
    } else {
      goodbisector = bisector_linf(sr, sp);
      CGAL_SDG_DEBUG(std::cout
          << "debug: computevv brp r=" << sr << " p=" << sp << std::endl;);
      CGAL_SDG_DEBUG(std::cout << "debug: computevv brp ="
          << goodbisector << std::endl;);
    }

    Polychainline_2 bqr = bisector_linf(sq, sr);
    CGAL_SDG_DEBUG(std::cout << "debug: computevv bqr q="
        << sq << " r=" << sr << std::endl;);
    CGAL_SDG_DEBUG(std::cout << "debug: computevv bqr ="
        << bqr << std::endl;);

    vv = goodbisector.first_intersection_point_with(bqr);
    CGAL_SDG_DEBUG(std::cout
        << "debug: computevv PSS vv=" << vv << std::endl;);
  }



  //--------------------------------------------------------------------------
  // the Voronoi vertex of three segments
  //--------------------------------------------------------------------------

  void
  compute_vv(const Site_2& sp, const Site_2& sq, const Site_2& sr,
             const SSS_Type&) const
  {
    CGAL_precondition( sp.is_segment() && sq.is_segment() &&
                       sr.is_segment() );

    if ( is_vv_computed ) { return; }
    is_vv_computed = true;

    const bool is_psrc_q = is_endpoint_of(sp.source_site(), sq);
    const bool is_psrc_r = is_endpoint_of(sp.source_site(), sr);
    const bool is_ptrg_q = is_endpoint_of(sp.target_site(), sq);
    const bool is_ptrg_r = is_endpoint_of(sp.target_site(), sr);

    if (is_psrc_q && is_psrc_r) {
      vv = sp.source();
    } else if (is_ptrg_q && is_ptrg_r) {
      vv = sp.target();
    } else {
      // here, not all segments have a common point
      const bool is_p_hor = is_site_horizontal(sp);
      const bool is_q_hor = is_site_horizontal(sq);
      const bool is_r_hor = is_site_horizontal(sr);

      const bool is_p_hv = is_p_hor || is_site_vertical(sp);
      const bool is_q_hv = is_q_hor || is_site_vertical(sq);
      const bool is_r_hv = is_r_hor || is_site_vertical(sr);

      if (is_p_hv && is_q_hv && is_r_hv) {
        return compute_vv_sss_hv(sp, sq, sr, is_p_hor, is_q_hor, is_r_hor);
      }

      Line_2 lines[3];
      orient_lines_linf(sp, sq, sr, lines);

      const bool have_common_pq = is_psrc_q || is_ptrg_q;
      const bool have_common_rp = is_psrc_r || is_ptrg_r;

#ifndef CGAL_NO_ASSERTIONS
      const bool is_qsrc_r = is_endpoint_of(sq.source_site(), sr);
      const bool is_qtrg_r = is_endpoint_of(sq.target_site(), sr);
      const bool have_common_qr = is_qsrc_r || is_qtrg_r;

      const unsigned int num_common =
        ((have_common_pq) ? 1 : 0) +
        ((have_common_qr) ? 1 : 0) +
        ((have_common_rp) ? 1 : 0)  ;

      // num_common can be 3 if the three segments create a triangle;
      // trivial assertion
      CGAL_assertion(num_common <= 3);

      const unsigned int num_hv =
        ((is_p_hv) ? 1 : 0) +
        ((is_q_hv) ? 1 : 0) +
        ((is_r_hv) ? 1 : 0)  ;

      // this is a trivial assertion
      CGAL_assertion(num_hv <= 3);

      CGAL_SDG_DEBUG(
        std::cout
            << "debug: vsqr num_common=" << num_common
            << " pq=" << have_common_pq
            << " qr=" << have_common_qr
            << " rp=" << have_common_rp
            << " num_hv=" << num_hv
            << std::endl;
      );
#endif

      bool bpqset(false);
      bool bqrset(false);
      bool brpset(false);

      Line_2 bpq;
      if ((is_p_hv && is_q_hv && have_common_pq) ) {
        const Point_2 xpq = is_psrc_q ? sp.source() : sp.target();
        Direction_2 dirbpq = dir_from_lines(lines[0], lines[1]);
        bpq = compute_line_dir(xpq, dirbpq);
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr bpq p=" << sp << " q=" << sq << std::endl;);
        //CGAL_SDG_DEBUG(std::cout
        //    << "debug: vsqr bpq =" << bpq << std::endl;);
        bpqset = true;
      }

      Line_2 bqr;
      if ((! bpqset) || (is_q_hv && is_r_hv) ||
          (! have_common_rp)) {
        bqr = bisector_linf_line(sq, sr, lines[1], lines[2]);
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr bqr q=" << sq << " r=" << sr << std::endl;);
        //CGAL_SDG_DEBUG(std::cout
        //    << "debug: vsqr bqr =" << bqr << std::endl;);
        bqrset = true;
      }

      Line_2 brp;
      if ((! (bpqset && bqrset))) {
        if (are_parallel_lines(lines[0], lines[2])) {
          brp = parallel_bis(lines[0], lines[2]);
        } else {
          Point_2 xrp;
          if (have_common_rp) {
            xrp = is_psrc_r ? sp.source() : sp.target();
          } else {
            RT hx, hy, hz;
            compute_intersection_of_lines(lines[0], lines[2], hx, hy, hz);
            xrp = Point_2(hx, hy, hz);
          }
          Direction_2 dirbrp = dir_from_lines(lines[2], lines[0]);
          brp = compute_line_dir(xrp, dirbrp);
        }
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr brp r=" << sr << " p=" << sp << std::endl;);
        //CGAL_SDG_DEBUG(std::cout
        //    << "debug: vsqr brp =" << brp << std::endl;);
        brpset = true;
      }

      CGAL_assertion((bpqset && bqrset) || (bqrset && brpset)
          || (brpset && bpqset));

      RT ux, uy, uz;
      if (bpqset && bqrset) {
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr SSS using bpq bqr" << std::endl;);
        compute_intersection_of_lines(bpq, bqr, ux, uy, uz);
      } else if (bqrset && brpset) {
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr SSS using bqr brp" << std::endl;);
        compute_intersection_of_lines(bqr, brp, ux, uy, uz);
      } else {
        CGAL_SDG_DEBUG(std::cout
            << "debug: vsqr SSS using brp bpq" << std::endl;);
        compute_intersection_of_lines(brp, bpq, ux, uy, uz);
      }
      vv = Point_2(ux, uy, uz);
      CGAL_SDG_DEBUG(std::cout
          << "debug: vsqr SSS vv=" << vv << std::endl;);
      CGAL_assertion( oriented_side_of_line(lines[0], this->point()) != ZERO );
      CGAL_assertion( oriented_side_of_line(lines[1], this->point()) != ZERO );
      CGAL_assertion( oriented_side_of_line(lines[2], this->point()) != ZERO );
    }
  }


  // SSS: all sites are axis-parallel
  inline void
  compute_vv_sss_hv(const Site_2 & p, const Site_2 & q, const Site_2 & r,
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
    const RT umid = prevc + nextc;
    const RT udis = RT(2)*hvseg_coord(odd, is_odd_hor) +
      RT(are_common_hor ? +1 : -1) * (prevc - nextc);
    vv = is_odd_hor ? Point_2(umid, udis, RT(2)) :
                      Point_2(udis, umid, RT(2)) ;
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Voronoi Linf radius computation
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  template<class Type>
  inline
  FT
  linf_radius(const Point_2& vv,
                 const Site_2& p, const Site_2& /*q*/, const Site_2& /*r*/,
                 const Type&) const
  {
    CGAL_precondition( p.is_point() );

    Point_2 pp = p.point();
    FT dx = CGAL::abs(vv.x() - pp.x());
    FT dy = CGAL::abs(vv.y() - pp.y());

    return (CGAL::max)(dx, dy);
  }

  inline
  FT
  linf_radius_pps(const Point_2& vv,
                  const Site_2& p, const Site_2& q,
                  const bool is_dx_max,
                  const FT & diffdvtx, const FT & diffdvty
                 ) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );

    const Sign sgnmaxdiff = CGAL::sign( is_dx_max ? diffdvtx : diffdvty );

    const Point_2 qq = q.point();
    const FT diffdvqx = vv.x() - qq.x();
    const FT diffdvqy = vv.y() - qq.y();
    const FT absdvqx = CGAL::abs(diffdvqx);
    const FT absdvqy = CGAL::abs(diffdvqy);
    const bool q_dx_max = CGAL::compare(absdvqx, absdvqy) == LARGER;
    if (is_dx_max == q_dx_max) {
      if (CGAL::sign( q_dx_max ? diffdvqx : diffdvqy ) == sgnmaxdiff) {
        return q_dx_max ? absdvqx : absdvqy;
      }
    }

    const Point_2 pp = p.point();
    const FT diffdvpx = vv.x() - pp.x();
    const FT diffdvpy = vv.y() - pp.y();
    const FT absdvpx = CGAL::abs(diffdvpx);
    const FT absdvpy = CGAL::abs(diffdvpy);
    return (CGAL::max)(absdvpx, absdvpy);
  }

  inline
  FT
  linf_radius(const Point_2& vv,
                 const Site_2& p, const Site_2& q, const Site_2& r,
                 const SSS_Type&) const
  {
    CGAL_assertion( p.is_segment() && q.is_segment() && r.is_segment() );
    CGAL_USE(q);
    CGAL_USE(r);

    Line_2 l = compute_supporting_line(p.supporting_site());
    Homogeneous_point_2 pref = compute_linf_projection_hom(l, vv);

    FT dx = CGAL::abs(vv.x() - pref.x());
    FT dy = CGAL::abs(vv.y() - pref.y());
    return (CGAL::max)(dx, dy);
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Voronoi Linf fine radius computation
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  inline
  FT
  linf_fine_radius(const Point_2& vv,
                   const Site_2& p, const Site_2& q, const Site_2& r,
                   const PPP_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_point() && r.is_point());

    Point_2 pp = p.point();
    FT minp = (CGAL::min)(CGAL::abs(vv.x() - pp.x()),
                          CGAL::abs(vv.y() - pp.y()));

    Point_2 qq = q.point();
    FT minq = (CGAL::min)(CGAL::abs(vv.x() - qq.x()),
                          CGAL::abs(vv.y() - qq.y()));

    Point_2 rr = r.point();
    FT minr = (CGAL::min)(CGAL::abs(vv.x() - rr.x()),
                          CGAL::abs(vv.y() - rr.y()));

    return (CGAL::max)(minp, (CGAL::max)(minq, minr));
  }

  inline
  FT
  linf_fine_radius(const Point_2& vv,
                   const Site_2& p, const Site_2& q, const Site_2& r,
                   const PPS_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_point() && r.is_segment());

    Point_2 pp = p.point();
    FT minp = (CGAL::min)(CGAL::abs(vv.x() - pp.x()),
                        CGAL::abs(vv.y() - pp.y()));

    Point_2 qq = q.point();
    FT minq = (CGAL::min)(CGAL::abs(vv.x() - qq.x()),
                        CGAL::abs(vv.y() - qq.y()));

    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);
    FT drx = CGAL::abs(vv.x() - rref.x());
    FT dry = CGAL::abs(vv.y() - rref.y());
    FT minr = (CGAL::min)(drx, dry);

    return (CGAL::max)(minp, (CGAL::max)(minq, minr));
  }

  inline
  FT
  linf_fine_radius(const Point_2& vv,
                   const Site_2& p, const Site_2& q, const Site_2& r,
                   const PSS_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_segment() && r.is_segment());

    Point_2 pp = p.point();
    FT minp = (CGAL::min)(CGAL::abs(vv.x() - pp.x()),
                        CGAL::abs(vv.y() - pp.y()));

    Line_2 lq = compute_supporting_line(q.supporting_site());
    Homogeneous_point_2 qref = compute_linf_projection_hom(lq, vv);
    FT dqx = CGAL::abs(vv.x() - qref.x());
    FT dqy = CGAL::abs(vv.y() - qref.y());
    FT minq = (CGAL::min)(dqx, dqy);

    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);
    FT drx = CGAL::abs(vv.x() - rref.x());
    FT dry = CGAL::abs(vv.y() - rref.y());
    FT minr = (CGAL::min)(drx, dry);

    return (CGAL::max)(minp, (CGAL::max)(minq, minr));
  }


  inline
  FT
  linf_fine_radius(const Point_2& vv,
                   const Site_2& p, const Site_2& q, const Site_2& r,
                   const SSS_Type&) const
  {
    CGAL_assertion( p.is_segment() && q.is_segment() && r.is_segment() );

    Line_2 lp = compute_supporting_line(p.supporting_site());
    Homogeneous_point_2 pref = compute_linf_projection_hom(lp, vv);
    FT dpx = CGAL::abs(vv.x() - pref.x());
    FT dpy = CGAL::abs(vv.y() - pref.y());
    FT minp = (CGAL::min)(dpx, dpy);

    Line_2 lq = compute_supporting_line(q.supporting_site());
    Homogeneous_point_2 qref = compute_linf_projection_hom(lq, vv);
    FT dqx = CGAL::abs(vv.x() - qref.x());
    FT dqy = CGAL::abs(vv.y() - qref.y());
    FT minq = (CGAL::min)(dqx, dqy);

    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);
    FT drx = CGAL::abs(vv.x() - rref.x());
    FT dry = CGAL::abs(vv.y() - rref.y());
    FT minr = (CGAL::min)(drx, dry);

    return (CGAL::max)(minp, (CGAL::max)(minq, minr));
  }

  // L_inf refinement for non-axis parallel lines
  template<class TypePSSorSSS> // for PSS and SSS return false
  inline
  Comparison_result
  linf_refine_nonhv(const Point_2& ,
              const Site_2& , const Site_2& , const Site_2& ,
              const Line_2& , Homogeneous_point_2& ,
              const TypePSSorSSS &
             ) const
  {
    return EQUAL;
  }

  // L_inf refinement for non-axis parallel lines PPS case
  inline
  Comparison_result
  linf_refine_nonhv(const Point_2& vv,
              const Site_2& p, const Site_2& q, const Site_2& r,
              const Line_2& l, Homogeneous_point_2& lref,
              const PPS_Type &
             ) const
  {
    CGAL_assertion(p.is_point());
    CGAL_assertion(q.is_point());
    CGAL_assertion(r.is_segment());
    CGAL_USE(vv);
    CGAL_USE(lref);
    if (points_inside_touching_sides_v(l, p, r, q, this->point())) {
      return LARGER;
    }
    return EQUAL;
  }

  // L_inf refinement for non-axis parallel lines PPP case
  inline
  Comparison_result
  linf_refine_nonhv(const Point_2& vv,
              const Site_2& p, const Site_2& q, const Site_2& r,
              const Line_2& l, Homogeneous_point_2& lref,
              const PPP_Type &
             ) const
  {
    CGAL_assertion(p.is_point());
    CGAL_assertion(q.is_point());
    CGAL_assertion(r.is_point());
    CGAL_USE(vv);
    CGAL_USE(lref);
    if (points_inside_touching_sides_v(l, p, q, r, this->point())) {
      return LARGER;
    }
    if (points_inside_touching_sides_v(l, q, r, p, this->point())) {
      return LARGER;
    }
    if (points_inside_touching_sides_v(l, r, p, q, this->point())) {
      return LARGER;
    }
    return EQUAL;
  }

  // L_inf refinement
  template<class Type>
  inline
  Comparison_result
  linf_refine(const Point_2& vv,
              const Site_2& p, const Site_2& q, const Site_2& r,
              const Line_2& l, Homogeneous_point_2& lref,
              const Type & type
             ) const
  {
    const bool is_l_h_or_v = is_line_h_or_v(l);

    if (! is_l_h_or_v) {
      return linf_refine_nonhv(vv, p, q, r, l, lref, type);
    }

    FT difxvl = vv.x() - lref.x();
    FT difyvl = vv.y() - lref.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);
    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    CGAL_SDG_DEBUG(std::cout << "debug linf_refine cmplabsxy = "
      << cmplabsxy << std::endl;);

    if ((cmplabsxy == EQUAL) && is_l_h_or_v) {
      return POSITIVE;
    }

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    Oriented_side oslvv (ON_ORIENTED_BOUNDARY);
    if (p.is_segment() || q.is_segment() || r.is_segment()) {
      oslvv = oriented_side_of_line(l, vv);
      CGAL_assertion(oslvv != ON_ORIENTED_BOUNDARY);
    }

    if (p.is_point()) {
      Point_2 pp = p.point();
      FT difxvp = vv.x() - pp.x();
      FT difyvp = vv.y() - pp.y();
      FT absdifxvp = CGAL::abs(difxvp);
      FT absdifyvp = CGAL::abs(difyvp);
      Comparison_result cmppabsxy = CGAL::compare(absdifxvp, absdifyvp);
      CGAL_SDG_DEBUG(std::cout << "debug linf_refine cmplabsxy = "
          << cmplabsxy << " cmppabsxy=" << cmppabsxy << std::endl;);
      if (cmplabsxy != EQUAL) {
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
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslpsrc =
          oriented_side_of_line(l, p.source_site().point());
        Oriented_side oslptrg =
          oriented_side_of_line(l, p.target_site().point());
        if (((oslpsrc != oslvv) && (oslptrg != oslvv)) &&
            ((oslpsrc != ON_ORIENTED_BOUNDARY) ||
             (oslptrg != ON_ORIENTED_BOUNDARY)   )         ) {
          compare_p = SMALLER;
        }
      }
    }

    if (q.is_point()) {
      Point_2 qq = q.point();
      FT difxvq = vv.x() - qq.x();
      FT difyvq = vv.y() - qq.y();
      FT absdifxvq = CGAL::abs(difxvq);
      FT absdifyvq = CGAL::abs(difyvq);
      Comparison_result cmpqabsxy = CGAL::compare(absdifxvq, absdifyvq);
      CGAL_SDG_DEBUG(std::cout << "debug linf_refine cmplabsxy = "
          << cmplabsxy << " cmpqabsxy=" << cmpqabsxy << std::endl;);
      if (cmplabsxy != EQUAL) {
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
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslqsrc =
          oriented_side_of_line(l, q.source_site().point());
        Oriented_side oslqtrg =
          oriented_side_of_line(l, q.target_site().point());
        if (((oslqsrc != oslvv) && (oslqtrg != oslvv)) &&
            ((oslqsrc != ON_ORIENTED_BOUNDARY) ||
             (oslqtrg != ON_ORIENTED_BOUNDARY)   )         ) {
          compare_q = SMALLER;
        }
      }
    }

    if (r.is_point()) {
      Point_2 rr = r.point();
      FT difxvr = vv.x() - rr.x();
      FT difyvr = vv.y() - rr.y();
      FT absdifxvr = CGAL::abs(difxvr);
      FT absdifyvr = CGAL::abs(difyvr);
      Comparison_result cmprabsxy = CGAL::compare(absdifxvr, absdifyvr);
      CGAL_SDG_DEBUG(std::cout << "debug linf_refine cmplabsxy = "
          << cmplabsxy << " cmprabsxy=" << cmprabsxy << std::endl;);
      if (cmplabsxy != EQUAL) {
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
    } else {
      if (is_l_h_or_v) {
        Oriented_side oslrsrc =
          oriented_side_of_line(l, r.source_site().point());
        Oriented_side oslrtrg =
          oriented_side_of_line(l, r.target_site().point());
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

  inline
  Comparison_result
  linf_refinement(const Point_2& vv,
                  const Site_2& p, const Site_2& q, const Site_2& r,
                  const Line_2& l, Homogeneous_point_2& lref,
                  const PPP_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_point() && r.is_point());
    CGAL_USE(l);

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    FT difxvl = vv.x() - lref.x();
    FT difyvl = vv.y() - lref.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);

    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    Point_2 pp = p.point();
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

    Point_2 qq = q.point();
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

    Point_2 rr = r.point();
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

    CGAL_SDG_DEBUG(std::cout << "debug compare PPP p q r = "
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) ||
        (compare_q == SMALLER) ||
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    if ((compare_p == LARGER) ||
        (compare_q == LARGER) ||
        (compare_r == LARGER)   ) {
      // tocheck
      return LARGER;
    }
    return EQUAL;
  }


  inline
  Comparison_result
  linf_refinement(const Point_2& vv,
                  const Site_2& p, const Site_2& q, const Site_2& r,
                  const Line_2& l, Homogeneous_point_2& lref,
                  const PPS_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_point() && r.is_segment());
    CGAL_USE(l);
    CGAL_USE(r);

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    FT difxvl = vv.x() - lref.x();
    FT difyvl = vv.y() - lref.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);

    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    Point_2 pp = p.point();
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

    Point_2 qq = q.point();
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

    /* do nothing for line */
    /*
    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);

    Point_2 rr (rref.x(), rref.y());
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
    */

    CGAL_SDG_DEBUG(std::cout << "debug compare PPS p q r = "
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) ||
        (compare_q == SMALLER) ||
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    if ((compare_p == LARGER) ||
        (compare_q == LARGER) ||
        (compare_r == LARGER)   ) {
      // tocheck
      return LARGER;
    }
    return EQUAL;
  }


  inline
  Comparison_result
  linf_refinement(const Point_2& vv,
                  const Site_2& p, const Site_2& q, const Site_2& r,
                  const Line_2& l, Homogeneous_point_2& lref,
                  const PSS_Type&) const
  {
    CGAL_precondition(
        p.is_point() && q.is_segment() && r.is_segment());
    CGAL_USE(l);
    CGAL_USE(q);
    CGAL_USE(r);

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    FT difxvl = vv.x() - lref.x();
    FT difyvl = vv.y() - lref.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);

    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    CGAL_SDG_DEBUG(std::cout << "debug vv=" << vv << std::endl;);

    Point_2 pp = p.point();
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
        CGAL_SDG_DEBUG(std::cout << "debug difyvl==difyvp" << std::endl;);
        CGAL_assertion(compare_p == EQUAL);
        compare_p = CGAL::compare(absdifxvl, absdifxvp);
      }
    }

    /* do nothing for line */
    /*
    Line_2 lq = compute_supporting_line(q.supporting_site());
    Homogeneous_point_2 qref = compute_linf_projection_hom(lq, vv);

    Point_2 qq (qref.x(), qref.y());
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
    */

    /* do nothing for line */
    /*
    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);

    Point_2 rr (rref.x(), rref.y());
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
    */

    CGAL_SDG_DEBUG(std::cout << "debug compare PSS p q r = "
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
      return LARGER;
    }
    */
    return EQUAL;
  }

  inline
  Comparison_result
  linf_refinement(const Point_2& vv,
                  const Site_2& p, const Site_2& q, const Site_2& r,
                  const Line_2& l, Homogeneous_point_2& lref,
                  const SSS_Type&) const
  {
    CGAL_precondition(
        p.is_segment() && q.is_segment() && r.is_segment());
    CGAL_USE(vv);
    CGAL_USE(l);
    CGAL_USE(lref);
    CGAL_USE(p);
    CGAL_USE(q);
    CGAL_USE(r);

    Comparison_result compare_p(EQUAL);
    Comparison_result compare_q(EQUAL);
    Comparison_result compare_r(EQUAL);

    /*
    FT difxvl = vv.x() - lref.x();
    FT difyvl = vv.y() - lref.y();
    FT absdifxvl = CGAL::abs(difxvl);
    FT absdifyvl = CGAL::abs(difyvl);
    Comparison_result cmplabsxy = CGAL::compare(absdifxvl, absdifyvl);
    */

    // philaris: (cmplabsxy == EQUAL) means that lref is
    // one of the corners of the square with center vv

    /* do nothing for line */
    /*
    Line_2 lp = compute_supporting_line(p.supporting_site());
    Homogeneous_point_2 pref = compute_linf_projection_hom(lp, vv);

    Point_2 pp (pref.x(), pref.y());
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
    */

    /* do nothing for line */
    /*
    Line_2 lq = compute_supporting_line(q.supporting_site());
    Homogeneous_point_2 qref = compute_linf_projection_hom(lq, vv);

    Point_2 qq (qref.x(), qref.y());
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
    */

    /* do nothing for line */
    /*
    Line_2 lr = compute_supporting_line(r.supporting_site());
    Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);

    Point_2 rr (rref.x(), rref.y());
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
    */

    CGAL_SDG_DEBUG(std::cout << "debug compare SSS p q r = "
      << compare_p << " " << compare_q << " " << compare_r << std::endl;);

    if ((compare_p == SMALLER) ||
        (compare_q == SMALLER) ||
        (compare_r == SMALLER)   ) {
      return SMALLER;
    }
    if ((compare_p == LARGER) ||
        (compare_q == LARGER) ||
        (compare_r == LARGER)   ) {
      return LARGER;
    }
    return EQUAL;
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // the incircle test --- start
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------


  // given the Voronoi vertex vv or p, q and r, returns the result of
  // the incircle test when the query object t is a point
  template<class Type>
  Sign incircle_p(const Point_2& vv,
                  const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_point() );

    CGAL_SDG_DEBUG(std::cout
        << "debug vsqr incircle_p entering pqr = ("
        << p << ", " << q << ", " << r << "), " << "t=" << t
        << " with known vv=" << vv << std::endl;);

    CGAL_assertion(r.is_segment()); // the PPP case is handled elsewhere

    // p segment implies q segment
    CGAL_assertion(p.is_point() || q.is_segment());

    const Point_2 tt = t.point();
    const FT diffdvtx = vv.x() - tt.x();
    const FT diffdvty = vv.y() - tt.y();

    CGAL_SDG_DEBUG(std::cout << "debug diffdvtx=" << diffdvtx
      << " diffdvty=" << diffdvty << std::endl;);

    const FT absdvtx = CGAL::abs(diffdvtx);
    const FT absdvty = CGAL::abs(diffdvty);
    const FT d = (CGAL::max)(absdvtx, absdvty);
    const FT radius = linf_radius(vv, p, q, r, type);

    Comparison_result crude = CGAL::compare(d, radius);

    CGAL_SDG_DEBUG(std::cout << "debug d=" << d
      << " radius=" << radius
      << " comparison=" << crude << std::endl;);

    if (crude != ZERO) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout << "debug vsqr refining in incircle_p pqr=("
        << p << ", " << q << ", " << r << "), "
        << "t=" << t
        << std::endl;);
      // here crude == ZERO, so
      // we might have to refine

      Sign retval (ZERO);

      bool is_p_hor (false);
      bool is_p_ver (false);
      bool is_q_hor (false);
      bool is_q_ver (false);
      bool is_r_hor (false);
      bool is_r_ver (false);
      bool is_t_endp_of_q (false);
      bool is_t_endp_of_r (false);
      bool is_p_endp_of_q (false);
      bool is_p_endp_of_r (false);
      bool is_q_endp_of_r (false);
      bool is_t_endp_of_p (false);

      if (p.is_segment()) {
        is_p_hor = is_site_horizontal(p);
        is_p_ver = is_site_vertical(p);
        is_t_endp_of_p = is_endpoint_of(t,p);
      }
      if (q.is_segment()) {
        is_q_hor = is_site_horizontal(q);
        is_q_ver = is_site_vertical(q);
        is_t_endp_of_q = is_endpoint_of(t,q);
        if (p.is_point()) {
          is_p_endp_of_q = is_endpoint_of(p,q);
        }
      }

      // r is segment anyway
      is_r_hor = is_site_horizontal(r);
      is_r_ver = is_site_vertical(r);
      is_t_endp_of_r = is_endpoint_of(t,r);
      if (p.is_point()) {
        is_p_endp_of_r = is_endpoint_of(p,r);
      }
      if (q.is_point()) {
        is_q_endp_of_r = is_endpoint_of(q,r);
      }

      // check if t is endpoint of a hor/ver segment;
      // this code is only reached in validity tests
      if (is_t_endp_of_p && (is_p_hor || is_p_ver)) {
        retval = (CGAL::sign(is_p_hor ? diffdvtx : diffdvty) != ZERO ) ?
          POSITIVE : ZERO;
        CGAL_SDG_DEBUG(std::cout << "debug vsqr t on p hor/ver retval="
            << retval << std::endl;);
        return retval;
      }
      if (is_t_endp_of_q && (is_q_hor || is_q_ver)) {
        retval = (CGAL::sign(is_q_hor ? diffdvtx : diffdvty) != ZERO ) ?
          POSITIVE : ZERO;
        CGAL_SDG_DEBUG(std::cout << "debug vsqr t on q hor/ver retval="
            << retval << std::endl;);
        return retval;
      }
      if (is_t_endp_of_r && (is_r_hor || is_r_ver)) {
        retval = (CGAL::sign(is_r_hor ? diffdvtx : diffdvty) != ZERO ) ?
          POSITIVE : ZERO;
        CGAL_SDG_DEBUG(std::cout << "debug vsqr t on r hor/ver retval="
            << retval << std::endl;);
        return retval;
      }

      FT d_fine = (CGAL::min)(absdvtx, absdvty);

      CGAL_SDG_DEBUG(std::cout
          << "debug d=" << d << " d_fine=" << d_fine << std::endl;);

      Point_2 pref, qref, rref;

      FT diffdvpx;
      FT diffdvpy;
      FT diffdvqx;
      FT diffdvqy;
      FT diffdvrx;
      FT diffdvry;


      const bool is_p_hv = is_p_hor || is_p_ver;
      const bool is_q_hv = is_q_hor || is_q_ver;
      const bool is_r_hv = is_r_hor || is_r_ver;

      if (p.is_segment() && q.is_segment()) {
        const bool is_psrc_q = is_endpoint_of(p.source_site(), q);
        const bool is_ptrg_q = is_endpoint_of(p.target_site(), q);
        if (is_psrc_q || is_ptrg_q) {
          if ((is_p_hv && (! is_q_hv)) ||
              (is_q_hv && (! is_p_hv))   ) {
            CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                << "p, q candidates" << std::endl; );
            if (is_p_hor || is_q_hor) {
              if (scmpx(is_psrc_q? p.source_site(): p.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                pref = (is_psrc_q?
                          p.source_site(): p.target_site()).point();
                diffdvpy = vv.y() - pref.y();
                Comparison_result test =
                  CGAL::compare(absdvty, CGAL::abs(diffdvpy));
                return (test == SMALLER) ? NEGATIVE : ZERO;

              }
            } else { // one of p, q is vertical
              if (scmpy(is_psrc_q? p.source_site(): p.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                pref = (is_psrc_q?
                          p.source_site(): p.target_site()).point();
                diffdvpx = vv.x() - pref.x();
                CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                    << "diffdvpx=" << diffdvpx
                    << " absdvtx=" << absdvtx << std::endl; );
                Comparison_result test =
                  CGAL::compare(absdvtx, CGAL::abs(diffdvpx));
                return (test == SMALLER) ? NEGATIVE : ZERO;

              }
            }
          }
        }
      }

      if (q.is_segment()) {
        const bool is_qsrc_r = is_endpoint_of(q.source_site(), r);
        const bool is_qtrg_r = is_endpoint_of(q.target_site(), r);
        if (is_qsrc_r || is_qtrg_r) {
          if ((is_site_h_or_v(q) && (! is_site_h_or_v(r))) ||
              (is_site_h_or_v(r) && (! is_site_h_or_v(q)))   ) {
            CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                << "q, r candidates" << std::endl; );
            if (is_q_hor || is_r_hor) {
              if (scmpx(is_qsrc_r? q.source_site(): q.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                qref = (is_qsrc_r?
                          q.source_site(): q.target_site()).point();
                diffdvqy = vv.y() - qref.y();
                Comparison_result test =
                  CGAL::compare(absdvty, CGAL::abs(diffdvqy));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            } else { // one of q, r is vertical
              if (scmpy(is_qsrc_r? q.source_site(): q.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                    << "vertical case" << std::endl; );
                qref = (is_qsrc_r?
                          q.source_site(): q.target_site()).point();
                diffdvqx = vv.x() - qref.x();
                CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                    << "diffdvqx=" << diffdvqx
                    << " absdvtx=" << absdvtx << std::endl; );
                Comparison_result test =
                  CGAL::compare(absdvtx, CGAL::abs(diffdvqx));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            }
          }
        }
      }

      if (p.is_segment()) {
        const bool is_rsrc_p = is_endpoint_of(r.source_site(), p);
        const bool is_rtrg_p = is_endpoint_of(r.target_site(), p);
        if (is_rsrc_p || is_rtrg_p) {
          if ((is_site_h_or_v(r) && (! is_site_h_or_v(p))) ||
              (is_site_h_or_v(p) && (! is_site_h_or_v(r)))   ) {
            CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                << "r, p candidates" << std::endl; );
            if (is_r_hor || is_p_hor) {
              if (scmpx(is_rsrc_p? r.source_site(): r.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                rref = (is_rsrc_p?
                          r.source_site(): r.target_site()).point();
                diffdvry = vv.y() - rref.y();
                Comparison_result test =
                  CGAL::compare(absdvty, CGAL::abs(diffdvry));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            } else { // one of r, p is vertical
              if (scmpy(is_rsrc_p? r.source_site(): r.target_site(), t)
                  == EQUAL)
              {
                // return NEGATIVE or ZERO
                CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                    << "vertical case" << std::endl; );
                rref = (is_rsrc_p?
                          r.source_site(): r.target_site()).point();
                diffdvrx = vv.x() - rref.x();
                CGAL_SDG_DEBUG(std::cout << "debug vsqr "
                    << "diffdvrx=" << diffdvrx
                    << " absdvtx=" << absdvtx << std::endl; );
                Comparison_result test =
                  CGAL::compare(absdvtx, CGAL::abs(diffdvrx));
                return (test == SMALLER) ? NEGATIVE : ZERO;
              }
            }
          }
        }
      }

      // check if p, t are endpoints of different segments
      // among q, r
      bool pt_endps_of_diff_qr =
        (is_t_endp_of_q && is_p_endp_of_r) ||
        (is_t_endp_of_r && is_p_endp_of_q)   ;

      CGAL_SDG_DEBUG( std::cout << "debug pt_endps_of_diff_qr = "
          << pt_endps_of_diff_qr << std::endl;);

      if (p.is_point()) {
        pref = p.point();
        diffdvpx = vv.x() - pref.x();
        diffdvpy = vv.y() - pref.y();
        if (v_type == PSS) {
          Comparison_result test (EQUAL);
          // check if p and t lie on the same side of the Linf-square
          if (  (CGAL::compare(diffdvpx, diffdvtx) == EQUAL) ) {
            CGAL_SDG_DEBUG(std::cout << "debug on same vertical side "
                << " p=" << p << " t=" << t << std::endl;);
            FT absdvpy = CGAL::abs(diffdvpy);
            CGAL_SDG_DEBUG(std::cout << "debug vsqr absdvty=" << absdvty
                << " absdvpy=" << absdvpy << std::endl;);
            CGAL_SDG_DEBUG(std::cout << "debug vsqr abs diff ty py ="
                << absdvty - absdvpy << std::endl;);
            if (pt_endps_of_diff_qr) {
              test = EQUAL;
            } else {
              test = CGAL::compare(absdvty, absdvpy);
            }
          } else if (CGAL::compare(diffdvpy, diffdvty) == EQUAL) {
            CGAL_SDG_DEBUG(std::cout << "debug on same horizontal side "
                << " p=" << p << " t=" << t << std::endl;);
            FT absdvpx = CGAL::abs(diffdvpx);
            if (pt_endps_of_diff_qr) {
              test = EQUAL;
            } else {
              test = CGAL::compare(absdvtx, absdvpx);
            }
          }
          CGAL_SDG_DEBUG(std::cout << "debug test=" << test << std::endl;);

          if (test == SMALLER) {
            return NEGATIVE;
          } else if (test == LARGER) {
            return POSITIVE;
          }

          if ((v_type == PSS) && (! is_q_hv) && (! is_r_hv)) {
            if (test == EQUAL) {
              if (points_inside_touching_sides_v(q, p, r, t, vv)) {
                return NEGATIVE;
              }
              if (points_inside_touching_sides_v(r, p, q, t, vv)) {
                return NEGATIVE;
              }
              CGAL_SDG_DEBUG(std::cout
                  << "debug equivalent points and two non-hv segments,"
                  << " thus return zero" << std::endl;);
              return ZERO;
            }
          }

        }
      } else {
        // tocheck and tofix
        CGAL_assertion(
            p.is_segment() && q.is_segment());
        return ZERO;
      }

      CGAL_assertion(p.is_point());

      CGAL_SDG_DEBUG(std::cout << "debug diffdvpx=" << diffdvpx
        << " diffdvpy=" << diffdvpy << std::endl;);

      if (CGAL::compare(diffdvpx, diffdvtx) == EQUAL) {
        CGAL_SDG_DEBUG(std::cout << "debug diffdvpx="
            << "diffdvtx=" << diffdvpx << std::endl;);
        if (CGAL::compare(CGAL::abs(diffdvpx), d) == EQUAL) {
          if (pt_endps_of_diff_qr) {
            retval = ZERO;
          } else {
            retval = CGAL::compare(d_fine, CGAL::abs(diffdvpy));
            CGAL_SDG_DEBUG(std::cout << "debug d_fine=" << d_fine
              << " absdiffdvpy=" << CGAL::abs(diffdvpy)
              << " comparison=" << retval << std::endl;);
          }
        }
      }
      if (CGAL::compare(diffdvpy, diffdvty) == EQUAL) {
        if (CGAL::compare(CGAL::abs(diffdvpy), d) == EQUAL) {
          if (pt_endps_of_diff_qr) {
            retval = ZERO;
          } else {
            retval = CGAL::compare(d_fine, CGAL::abs(diffdvpx));
            CGAL_SDG_DEBUG(std::cout << "debug d_fine=" << d_fine
              << " absdiffdvpx=" << CGAL::abs(diffdvpx)
              << " comparison=" << retval << std::endl;);
          }
        }
      }
      if (retval == SMALLER) {
        return NEGATIVE;
      }
      if (retval == LARGER) {
        return POSITIVE;
      }

      if (q.is_point()) {
        qref = q.point();
        diffdvqx = vv.x() - qref.x();
        diffdvqy = vv.y() - qref.y();
      } else {
        // tocheck and tofix
        // here q and r are segments and p is point
        CGAL_assertion(p.is_point());
        if ((is_q_hv && ! is_r_hv) ||
            (is_r_hv && ! is_q_hv)   ) {
          Line_2 lnap;
          Homogeneous_point_2 sref;
          bool samex;
          if (is_q_hv && ! is_r_hv) {
            CGAL_SDG_DEBUG(std::cout << "debug q:ap r:non-ap"
                << std::endl;);
            if (is_p_endp_of_q) {
              samex = is_q_ver ? true : false;
              lnap = compute_supporting_line(r.supporting_site());
            } else {
              return ZERO;
            }
          } else
          { // here, we have: (is_r_hv && ! is_q_hv)
            CGAL_SDG_DEBUG(std::cout << "debug q:non-ap r:ap"
                << std::endl;);
            if (is_p_endp_of_r) {
              samex = is_r_ver ? true : false;
              lnap = compute_supporting_line(q.supporting_site());
            } else {
              return ZERO;
            }
          }
          sref = compute_linf_projection_hom(lnap, vv);
          CGAL_SDG_DEBUG(std::cout << "debug sref ="
               << sref.x() << ' ' << sref.y() << std::endl;);
          if (samex) {
            CGAL_SDG_DEBUG(std::cout << "debug samex case" << std::endl;);
            FT diffdvsy = vv.y() - sref.y();
            if (CGAL::sign(diffdvsy) == CGAL::sign(diffdvty)) {
              if (CGAL::compare(CGAL::abs(diffdvtx),
                                CGAL::abs(diffdvty)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of samex case
          else { // samey case
            CGAL_SDG_DEBUG(std::cout << "debug samey case" << std::endl;);
            FT diffdvsx = vv.x() - sref.x();
            if (CGAL::sign(diffdvsx) == CGAL::sign(diffdvtx)) {
              if (CGAL::compare(CGAL::abs(diffdvty),
                                CGAL::abs(diffdvtx)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of samey case
        }
        Line_2 lq = compute_supporting_line(q.supporting_site());
        //qref = compute_linf_projection_nonhom(lq, vv);
        Homogeneous_point_2 hqref = compute_linf_projection_hom(lq, vv);
        diffdvqx = vv.x() - hqref.x();
        diffdvqy = vv.y() - hqref.y();
      }

      if (q.is_point()) {
        CGAL_SDG_DEBUG(std::cout << "debug diffdvqx=" << diffdvqx
            << " diffdvqy=" << diffdvqy << std::endl;);

        if (CGAL::compare(diffdvqx, diffdvtx) == EQUAL) {
          CGAL_SDG_DEBUG(std::cout << "debug diffdvqx="
              << " diffdvtx=" << diffdvtx
              << std::endl;);
          if (CGAL::compare(CGAL::abs(diffdvqx), d) == EQUAL) {
            retval = CGAL::compare(d_fine, CGAL::abs(diffdvqy));
            CGAL_SDG_DEBUG(std::cout << "debug d_fine=" << d_fine
                << " absdiffdvqy=" << CGAL::abs(diffdvqy)
                << " comparison=" << retval << std::endl;);
          }
        }
        if (CGAL::compare(diffdvqy, diffdvty) == EQUAL) {
          if (CGAL::compare(CGAL::abs(diffdvqy), d) == EQUAL) {
            retval = CGAL::compare(d_fine, CGAL::abs(diffdvqx));
            CGAL_SDG_DEBUG(std::cout << "debug d_fine=" << d_fine
                << " absdiffdvqx=" << CGAL::abs(diffdvqx)
                << " comparison=" << retval << std::endl;);
          }
        }
        if (retval == SMALLER) {
          return NEGATIVE;
        }

        if (retval == LARGER) {
          return POSITIVE;
        }
      }

      if (q.is_segment() && (! is_q_hv)) {
        if (CGAL::compare(d_fine, d) == SMALLER) {
          CGAL_assertion(p.is_point());
          if (points_inside_touching_sides_v(q, p, r, t, vv)) {
            return NEGATIVE;
          }
          if (! is_r_hv) {
            if (points_inside_touching_sides_v(r, p, q, t, vv)) {
              return NEGATIVE;
            }
          }
        }
      }

      // check for p, q with same coordinate and r non-hv segment
      if ((! (is_r_hor || is_r_ver)) &&
          (! (is_p_endp_of_r || is_q_endp_of_r))
         ) {
        CGAL_SDG_DEBUG(std::cout << "debug r is non-axis parallel"
            << " and neither p nor q endpoints" << std::endl;);
        CGAL_SDG_DEBUG(std::cout << "debug d=" << d <<
            " d_fine=" << d_fine << std::endl;);
        bool pqsamex = CGAL::compare(diffdvpx, diffdvqx) == EQUAL;
        bool pqsamey (false);
        if (pqsamex) {
          CGAL_SDG_DEBUG(std::cout << "debug p, q have same x, "
              << "might be on same Linf vertical side"
              << std::endl;);
        } else {
          pqsamey = CGAL::compare(diffdvpy, diffdvqy) == EQUAL;
          if (pqsamey) {
            CGAL_SDG_DEBUG(std::cout << "debug p, q have same y, "
                << "might be on same Linf horizontal side" << std::endl;);
          }
        }
        if (pqsamex || pqsamey) {
          Line_2 lr = compute_supporting_line(r.supporting_site());
          Homogeneous_point_2 rref = compute_linf_projection_hom(lr, vv);
          if ( pqsamex && (CGAL::compare(pref.x(), rref.x()) == EQUAL) ) {
            FT diffdvry = vv.y() - rref.y();
            if (CGAL::sign(diffdvry) == CGAL::sign(diffdvty)) {
              if (CGAL::compare(CGAL::abs(diffdvtx),
                                CGAL::abs(diffdvty)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of pqsamex and rref same x case
          if ( pqsamey && (CGAL::compare(pref.y(), rref.y()) == EQUAL) ) {
            FT diffdvrx = vv.x() - rref.x();
            if (CGAL::sign(diffdvrx) == CGAL::sign(diffdvtx)) {
              if (CGAL::compare(CGAL::abs(diffdvty),
                                CGAL::abs(diffdvtx)) == SMALLER) {
                return NEGATIVE;
              }
            }
          } // end of pqsamey and rref same y case
        } // end of case: pqsamex or pqsamey
      } // end of non-hv segment r case with p, q non-endpoints of r

      CGAL_SDG_DEBUG(std::cout
          << "debug in refinement return final zero" << std::endl;);

      return ZERO;

      //FT radius_fine = linf_fine_radius(vv, p, q, r, type);

      //return CGAL::compare(d_fine, radius_fine);

    } // end of case of crude == ZERO
  }




  //--------------------------------------------------------------------------
  // the first three objects are points and the query object is also a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PPP_Type&) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_point() );
    CGAL_precondition( t.is_point() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_p entering" << std::endl;);
    CGAL_SDG_DEBUG(std::cout << "debug incircle_p (p q r t)= "
      << p << ' ' << q << ' ' << r << ' ' << t
      << std::endl;);

    Bounded_side bs =
      side_of_bounded_square(p.point(), q.point(), r.point(), t.point());

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
  // the first three objects are two points and a segment and the query
  // object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PPS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() );

    // easy degeneracies --- start

    // if t is one of p or q then we know the result which is ZERO
    if (  same_points(p, t) || same_points(q, t)  ) {
      return ZERO;
    }

    const bool endp_t_of_r = is_endpoint_of(t, r);
    if (endp_t_of_r) {
      if (! is_site_h_or_v(r)) {
        // For a non-axis-parallel segment, the possible return values
        // are either ZERO or POSITIVE. In the non-axis-parallel case,
        // the point t can be at the corner of the Linf-square and thus the
        // possibility of ZERO cannot be excluded. However, this code seems
        // to be reached only in validations and POSITIVE seems to be
        // acceptable.
        return POSITIVE;
      } else {
        // here r is axis-parallel:
        // if p or q are endpoints of r, return POSITIVE
        if (is_endpoint_of(p, r) || is_endpoint_of(q, r)) {
          return POSITIVE;
        }
      }
    } else {
      if (is_on_hv_seg_line(t, r)) { return POSITIVE; }
    }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the first three objects are a point and two segments and the query
  // object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PSS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() );

    // easy degeneracies --- start

    // if p is a common point for q and r, then the Voronoi vertex of
    // p, q, and r is p. Since t cannot be the same as p, the result
    // is POSITIVE
    if ( is_endpoint_of(p, q) && is_endpoint_of(p, r) ) {
      return POSITIVE;
    }

    // if p and t are the same point, then return ZERO
    if ( same_points(p, t) ) { return ZERO; }

    // if t is an endpoint of either q or r, then t has to be outside
    // the Voronoi circle and thus the result is POSITIVE
    // philaris: this is only reached in validity checks
    if ( is_endpoint_of(t, q) || is_endpoint_of(t, r) ) {
      return POSITIVE;
    }


    // philaris: addition for Linf

    if (is_endpoint_of(p, q)) {
      if (is_site_horizontal(q)) {
        if (scmpy(p,t) == EQUAL) {
          Site_2 other =
            same_points(p, q.source_site()) ?
            q.target_site() : q.source_site();
          if (scmpx(other, p) == scmpx(p, t)) {
            return POSITIVE;
          }
        }
      }
      if (is_site_vertical(q)) {
        if (scmpx(p,t) == EQUAL) {
          Site_2 other =
            same_points(p, q.source_site()) ?
            q.target_site() : q.source_site();
          if (scmpy(other, p) == scmpy(p, t)) {
            return POSITIVE;
          }
        }
      }
    }

    if (is_endpoint_of(p, r)) {
      if (is_site_horizontal(r)) {
        if (scmpy(p,t) == EQUAL) {
          Site_2 other =
            same_points(p, r.source_site()) ?
            r.target_site() : r.source_site();
          if (scmpx(other, p) == scmpx(p, t)) {
            return POSITIVE;
          }
        }
      }
      if (is_site_vertical(r)) {
        if (scmpx(p,t) == EQUAL) {
          Site_2 other =
            same_points(p, r.source_site()) ?
            r.target_site() : r.source_site();
          if (scmpy(other, p) == scmpy(p, t)) {
            return POSITIVE;
          }
        }
      }
    }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  // the first three objects are segments and the query object is a point
  //--------------------------------------------------------------------------

  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const SSS_Type& type) const
  {
    CGAL_precondition( p.is_segment() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_point() );

    // easy degeneracies --- start

    // if t is an endpoint of p, q or r, then t has to lie outside the
    // Voronoi circle of p, q and r and thus the result is positive
    if ( is_endpoint_of(t, p) || is_endpoint_of(t, q) ||
         is_endpoint_of(t, r) ) {
      return POSITIVE;
    }

    // easy degeneracies --- end

    compute_vv(p, q, r, type);
    return incircle_p(vv, p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the incircle test when the query object is a line
  //--------------------------------------------------------------------------

  template<class Type>
  inline
  Sign
  incircle_xxxl(const Point_2& vv,
                const Site_2& p, const Site_2& q, const Site_2& r,
                const Line_2& l, const Type& type) const
  {
    FT radius = linf_radius(vv, p, q, r, type);

    // compute Linf distance of vv from l
    Homogeneous_point_2 lref = compute_linf_projection_hom(l, vv);
    FT absdvlx = CGAL::abs(vv.x() - lref.x());
    FT absdvly = CGAL::abs(vv.y() - lref.y());
    FT d = (CGAL::max)(absdvlx, absdvly);

    Comparison_result crude = CGAL::compare(d, radius);

    if (crude != ZERO) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout << "debug vsqr refining in xxxl pqr=("
        << p << ", " << q << ", " << r << "), "
        << "lref=" << lref.x() << ' ' << lref.y()
        << ", l: " << l.a() << ' ' << l.b() << ' ' <<  l.c()
        << std::endl;);
      // here crude == ZERO, so
      // we might have to refine

      Comparison_result other =
        linf_refine(vv, p, q, r, l, lref, type);

      if (crude != other) {
        CGAL_SDG_DEBUG(std::cout << "xxxl instead of 0 returning " << other <<
          std::endl;);
      }

      return other;

    }
  }


  // philaris: I might have to change that,
  // it seems in L2 that an L2-perpendicular to l is computed
  // that goes through vv
  inline
  Oriented_side
  oriented_side_l2(const Point_2& vv, const Line_2& l, const Point_2& p) const
  {
    Line_2 l1(l.b(), -l.a(), l.a() * vv.y() - l.b() * vv.x());

    return oriented_side_of_line(l1, p);
  }

  // philaris: the linf version of the oriented_side predicate
  inline
  Oriented_side
  oriented_side_linf(const Point_2& vv, const Line_2& l, const Point_2& p) const
  {
    CGAL_SDG_DEBUG(std::cout << "debug oriented_side_linf " << std::endl;);

    Line_2 l1 = compute_linf_perpendicular(l, vv);

    return oriented_side_of_line(l1, p);
  }


  //--------------------------------------------------------------------------
  // generic incircle test when the query object is a segment
  //--------------------------------------------------------------------------

  // first check is at least one of the endpoints of s is in conflict
  // with the Voronoi circle; in this case return NEGATIVE
  // otherwise test against the supporting line l of t; if l does not
  // conflict with the Voronoi circle return POSITIVE, otherwise check
  // if the endpoints of t are on the same oriented side of the line
  // perpendicular to l, passing through the Voronoi vertex of p, q,
  // and r, and respond accordingly: if they are on the same side
  // there is no conflict, otherwise there is a conflict.
  template<class Type>
  Sign incircle_xxxs(const Site_2& p, const Site_2& q, const Site_2& r,
                     const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug fn incircle_xxxs pqrt= ("
        << p << ") ("
        << q << ") (" << r << ") (" << t << ")" << std::endl;);

    bool is_p_point = p.is_point();
    bool is_q_point = q.is_point();
    bool is_r_point = r.is_point();

    unsigned int numpts_in_pqr =
      ((is_p_point)? 1 : 0) +
      ((is_q_point)? 1 : 0) +
      ((is_r_point)? 1 : 0)  ;

    bool is_p_tsrc(false);
    if ( is_p_point && same_points(p, t.source_site()) ) {
      is_p_tsrc = true;
    }
    bool is_q_tsrc(false);
    if ( is_q_point && same_points(q, t.source_site()) ) {
      is_q_tsrc = true;
    }
    bool is_r_tsrc(false);
    if ( is_r_point && same_points(r, t.source_site()) ) {
      is_r_tsrc = true;
    }

    unsigned int numendpts_of_t = 0;

    Sign d1, d2;
    if ( is_p_tsrc || is_q_tsrc || is_r_tsrc ) {
      d1 = ZERO;
      ++numendpts_of_t;
    } else {
      d1 = incircle_p(p, q, r, t.source_site(), type);
    }

    if (  certainly(d1 == NEGATIVE)  ) { return NEGATIVE; }
    if (  !is_certain(d1 == NEGATIVE)  ) { return indeterminate<Sign>(); }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_xxxs d1=" << d1 << std::endl;);

    bool is_p_ttrg(false);
    if ( is_p_point && same_points(p, t.target_site()) ) {
      is_p_ttrg = true;
    }
    bool is_q_ttrg(false);
    if ( is_q_point && same_points(q, t.target_site()) ) {
      is_q_ttrg = true;
    }
    bool is_r_ttrg(false);
    if ( is_r_point && same_points(r, t.target_site()) ) {
      is_r_ttrg = true;
    }

    if ( is_p_ttrg || is_q_ttrg || is_r_ttrg ) {
      d2 = ZERO;
      ++numendpts_of_t;
    } else {
      d2 = incircle_p(p, q, r, t.target_site(), type);
    }

    if (  certainly( d2 == NEGATIVE )  ) { return NEGATIVE; }
    if (  !is_certain( d2 == NEGATIVE )  ) { return indeterminate<Sign>(); }

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_xxxs d2=" << d2 << std::endl;);

    CGAL_assertion(numendpts_of_t < 2);

    CGAL_SDG_DEBUG(std::cout << "debug incircle_xxxs numendpts_of_t= "
      << numendpts_of_t << std::endl;);

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
          if ((! is_p_point) && is_endpoint_of(endp, p)) {
            numothers++;
            other = p;
          }

          if ((! is_q_point) && is_endpoint_of(endp, q)) {
            numothers++;
            other = q;
          }

          if ((! is_r_point) && is_endpoint_of(endp, r)) {
            numothers++;
            other = r;
          }
        }

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
          compute_vv(p, q, r, type);

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
    } // endif ((numendpts_of_t > 0)

    bool is_tsrc_endp_of_p (false);
    bool is_tsrc_endp_of_q (false);
    bool is_tsrc_endp_of_r (false);
    bool is_ttrg_endp_of_p (false);
    bool is_ttrg_endp_of_q (false);
    bool is_ttrg_endp_of_r (false);

    if (! is_p_point) {
      is_tsrc_endp_of_p = same_points(t.source_site(), p.source_site())
                       || same_points(t.source_site(), p.target_site());
      is_ttrg_endp_of_p = same_points(t.target_site(), p.source_site())
                       || same_points(t.target_site(), p.target_site());
    }
    if (! is_q_point) {
      is_tsrc_endp_of_q = same_points(t.source_site(), q.source_site())
                       || same_points(t.source_site(), q.target_site());
      is_ttrg_endp_of_q = same_points(t.target_site(), q.source_site())
                       || same_points(t.target_site(), q.target_site());
    }
    if (! is_r_point) {
      is_tsrc_endp_of_r = same_points(t.source_site(), r.source_site())
                       || same_points(t.source_site(), r.target_site());
      is_ttrg_endp_of_r = same_points(t.target_site(), r.source_site())
                       || same_points(t.target_site(), r.target_site());
    }

    if (is_tsrc_endp_of_p && is_tsrc_endp_of_q) {
      if (test_star(t.source_site(), p, q, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_p && is_ttrg_endp_of_q) {
      if (test_star(t.target_site(), p, q, t)) {
        return NEGATIVE;
      }
    }
    if (is_tsrc_endp_of_q && is_tsrc_endp_of_r) {
      if (test_star(t.source_site(), q, r, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_q && is_ttrg_endp_of_r) {
      if (test_star(t.target_site(), q, r, t)) {
        return NEGATIVE;
      }
    }
    if (is_tsrc_endp_of_r && is_tsrc_endp_of_p) {
      if (test_star(t.source_site(), r, p, t)) {
        return NEGATIVE;
      }
    }
    if (is_ttrg_endp_of_r && is_ttrg_endp_of_p) {
      if (test_star(t.target_site(), r, p, t)) {
        return NEGATIVE;
      }
    }

    bool same_slope_at_corner(false);
    if (numendpts_of_t > 0) {
      CGAL_assertion(numendpts_of_t == 1);
      CGAL_assertion(is_p_point);
      if (! is_r_point) { // there is at least one segment in p, q, r
        if (   (is_tsrc_endp_of_r && (is_p_tsrc || is_q_tsrc))
            || (is_ttrg_endp_of_r && (is_p_ttrg || is_q_ttrg))
           ) {
          if (have_same_slope(r, t)) {
            same_slope_at_corner = true;
          }
        }
        if (   (is_tsrc_endp_of_q && (is_p_tsrc))
            || (is_ttrg_endp_of_q && (is_p_ttrg))
           ) {
          if (have_same_slope(q, t)) {
            same_slope_at_corner = true;
          }
        }
      }
    }

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_xxxs: same_slope_at_corner="
        << same_slope_at_corner << std::endl;);

    Line_2 l = compute_supporting_line(t.supporting_site());
    compute_vv(p, q, r, type);
    Sign sl(ZERO);
    if (! same_slope_at_corner) {
      sl = incircle_xxxl(vv, p, q, r, l, type);
    }

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_xxxs: incircle_xxxl returned sl="
        << sl << std::endl;);

    if (  certainly( sl == POSITIVE )  ) { return sl; }
    if (  !is_certain( sl == POSITIVE )  ) { return indeterminate<Sign>(); }

    CGAL_SDG_DEBUG(std::cout << "debug incircle_xxxs sl=" << sl <<
      " d1=" << d1 << " d2=" << d2 << std::endl;);

    CGAL_SDG_DEBUG(std::cout
        << "debug numpts_in_pqr=" << numpts_in_pqr << std::endl;);

    // philaris: here we have a serious change related to L2
    if ( sl == ZERO && (d1 == ZERO || d2 == ZERO) ) {

      // if some site in {p,q,r} is a point and it is also:
      // an endpoint of t and an endpoint of another site in {p,q,r}

      // or if t has a common endpoint with a segment in {p,q,r}

      Site_2 sqpnt, other_t, other_seg;

      if (compute_helper(p, q, r, t, sqpnt, other_t, other_seg)) {

        CGAL_assertion(sqpnt.is_point());
        CGAL_assertion(other_t.is_point());
        CGAL_assertion(other_seg.is_point());

        CGAL_SDG_DEBUG(std::cout
          << "debug incircle_xxxs compute_helper true, "
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

    Oriented_side os1 = oriented_side_linf(vv, l, t.source());
    Oriented_side os2 = oriented_side_linf(vv, l, t.target());

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_xxxs: os1=" << os1 << " os2="
        << os2 << std::endl;);

    if ( sl == ZERO ) {
      if (os1 == ON_ORIENTED_BOUNDARY || os2 == ON_ORIENTED_BOUNDARY) {
        return ZERO;
      }
      return ( os1 == os2 ) ? POSITIVE : ZERO;
    }

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_xxxs non-zero sl=" << sl << " : os1="
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

    CGAL_SDG_DEBUG(std::cout
        << "debug compute_helper #pts=" << numpts << std::endl;);

    if (numpts == 3) {
      return false;
    }

    // here and on, there are 0, 1 or 2 points in {p,q,r}


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
           have_common_p_t(false),
           is_ptrg_endp_of_t(false);

      if (! is_p_point) {
        CGAL_assertion( ! same_segments(p, t) );
        const bool is_psrc_tsrc = same_points(p.source_site(), t.source_site());
        const bool is_ptrg_tsrc = same_points(p.target_site(), t.source_site());
        const bool is_psrc_ttrg = same_points(p.source_site(), t.target_site());
        const bool is_ptrg_ttrg = same_points(p.target_site(), t.target_site());
        have_common_p_tsrc = is_psrc_tsrc || is_ptrg_tsrc;
        have_common_p_ttrg = is_psrc_ttrg || is_ptrg_ttrg;
        have_common_p_t = have_common_p_tsrc || have_common_p_ttrg;
        is_ptrg_endp_of_t = is_ptrg_tsrc || is_ptrg_ttrg;
      }

      bool have_common_q_tsrc(false),
           have_common_q_ttrg(false),
           have_common_q_t(false),
           is_qtrg_endp_of_t(false);

      if (! is_q_point) {
        CGAL_assertion( ! same_segments(q, t) );
        const bool is_qsrc_tsrc = same_points(q.source_site(), t.source_site());
        const bool is_qtrg_tsrc = same_points(q.target_site(), t.source_site());
        const bool is_qsrc_ttrg = same_points(q.source_site(), t.target_site());
        const bool is_qtrg_ttrg = same_points(q.target_site(), t.target_site());
        have_common_q_tsrc = is_qsrc_tsrc || is_qtrg_tsrc;
        have_common_q_ttrg = is_qsrc_ttrg || is_qtrg_ttrg;
        have_common_q_t = have_common_q_tsrc || have_common_q_ttrg;
        is_qtrg_endp_of_t = is_qtrg_tsrc || is_qtrg_ttrg;
      }

      bool have_common_r_tsrc(false),
           have_common_r_ttrg(false),
           have_common_r_t(false),
           is_rtrg_endp_of_t(false);

      if (! is_r_point) {
        CGAL_assertion( ! same_segments(r, t) );
        const bool is_rsrc_tsrc = same_points(r.source_site(), t.source_site());
        const bool is_rtrg_tsrc = same_points(r.target_site(), t.source_site());
        const bool is_rsrc_ttrg = same_points(r.source_site(), t.target_site());
        const bool is_rtrg_ttrg = same_points(r.target_site(), t.target_site());
        have_common_r_tsrc = is_rsrc_tsrc || is_rtrg_tsrc;
        have_common_r_ttrg = is_rsrc_ttrg || is_rtrg_ttrg;
        have_common_r_t = have_common_r_tsrc || have_common_r_ttrg;
        is_rtrg_endp_of_t = is_rtrg_tsrc || is_rtrg_ttrg;
      }

      const unsigned int numcommon =
        ((have_common_p_t)? 1 : 0) +
        ((have_common_q_t)? 1 : 0) +
        ((have_common_r_t)? 1 : 0)  ;

      CGAL_SDG_DEBUG(std::cout
          << "debug compute_helper #numcommon="
          << numcommon << std::endl;);

      CGAL_assertion(numcommon < 3);

      if (numcommon == 0) {
        return false;
      }

      // here, numcommon equals 1 or 2

      const unsigned int numcommon_tsrc =
        ((have_common_p_tsrc)? 1 : 0) +
        ((have_common_q_tsrc)? 1 : 0) +
        ((have_common_r_tsrc)? 1 : 0)  ;

      const unsigned int numcommon_ttrg =
        ((have_common_p_ttrg)? 1 : 0) +
        ((have_common_q_ttrg)? 1 : 0) +
        ((have_common_r_ttrg)? 1 : 0)  ;


      if (numcommon == 1) {
        if (numcommon_tsrc > 0) {
          // here, numcommon_tsrc == 1
          sqpnt = t.source_site();
          other_of_t = t.target_site();
        } else {
          // here, numcommon_ttrg == 1
          sqpnt = t.target_site();
          other_of_t = t.source_site();
        }

        if (have_common_p_t) {
          other_of_seg = (is_ptrg_endp_of_t) ?
                         p.source_site() :
                         p.target_site();
        } else if (have_common_q_t) {
          other_of_seg = (is_qtrg_endp_of_t) ?
                         q.source_site() :
                         q.target_site();
        } else if (have_common_r_t) {
          other_of_seg = (is_rtrg_endp_of_t) ?
                         r.source_site() :
                         r.target_site();
        } else {
          CGAL_assertion(false);
        }

        return true;

      }

      CGAL_assertion( numcommon == 2 );

      CGAL_assertion( numcommon_tsrc + numcommon_ttrg == 2 );

      if ( numcommon_tsrc == numcommon_ttrg )  { // both equal 1
        CGAL_SDG_DEBUG(std::cout
            << "debug compute_helper #numcommon tsrc ttrg equal"
            << std::endl;);
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
        const bool is_p_qsrc = same_points(p, q.source_site());
        if (is_p_qsrc || same_points(p, q.target_site())) {
          other_of_seg = is_p_qsrc ? q.target_site() : q.source_site();
          return true;
        }
      }
      if (r.is_segment()) {
        const bool is_p_rsrc = same_points(p, r.source_site());
        if (is_p_rsrc || same_points(p, r.target_site())) {
          other_of_seg = is_p_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

    } // end of case: is_p_endp_of_t

    if (is_q_endp_of_t) {
      if (r.is_segment()) {
        const bool is_q_rsrc = same_points(q, r.source_site());
        if (is_q_rsrc || same_points(q, r.target_site())) {
          other_of_seg = is_q_rsrc ? r.target_site() : r.source_site();
          return true;
        }
      }

      if (p.is_segment()) {
        const bool is_q_psrc = same_points(q, p.source_site());
        if (is_q_psrc || same_points(q, p.target_site())) {
          other_of_seg = is_q_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

    } // end of case: is_q_endp_of_t

    if (is_r_endp_of_t) {
      if (p.is_segment()) {
        const bool is_r_psrc = same_points(r, p.source_site());
        if (is_r_psrc || same_points(r, p.target_site())) {
          other_of_seg = is_r_psrc ? p.target_site() : p.source_site();
          return true;
        }
      }

      if (q.is_segment()) {
        const bool is_r_qsrc = same_points(r, q.source_site());
        if (is_r_qsrc || same_points(r, q.target_site())) {
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
  // the first three objects are points and the query object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PPP_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_point() );
    CGAL_precondition( t.is_segment() );

    // easy degeneracies --- start

    // check if the endpoints of t are two of p, q and r
    unsigned int n_ends = 0;
    bool end_pt = is_endpoint_of(p, t);
    bool end_qt = is_endpoint_of(q, t);
    bool end_rt = is_endpoint_of(r, t);
    if ( end_pt ) ++n_ends;
    if ( end_qt ) ++n_ends;
    if ( end_rt ) ++n_ends;

    CGAL_assertion( n_ends < 3 );
    if ( n_ends == 2 ) { return NEGATIVE; }

    /*
#ifndef CGAL_DISABLE_AM_CODE
    // code added in previous version by Andreas + Monique -- start
    Site_2 const *pp1 = nullptr;
    if ( end_pt ) pp1 = &p;
    else if ( end_qt ) pp1 = &q;
    else if ( end_rt ) pp1 = &r;
    if ( pp1 != nullptr ) {
      // As the Voronoi circle and the segment t touch in p1,
      // it is enough to check that the center and the non-touching
      // point of the segment
      // are not in the same halfspace defined by the tangent line through p1
      Point_2 p1 = pp1->point();
      Point_2 p2 = other_site(*pp1, t).point();
      compute_vv(p, q, r, type);

      // philaris: this has to be changed

      CGAL_SDG_DEBUG(std::cout << "debug: unreachable" << std::endl;);

      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p1, p2 - p1) );
    }
    // code added in previous version by Andreas + Monique -- end
#endif // CGAL_DISABLE_AM_CODE
    */

    // easy degeneracies --- end

    return incircle_xxxs(p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  // the first three objects are two points and a segment and the
  // query object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PPS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_point() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    // easy degeneracies --- start

    // philaris: removed some of easy degeneracies

    // check if the endpoints of t are p and q
    bool end_pt = is_endpoint_of(p, t);
    bool end_qt = is_endpoint_of(q, t);

    if ( end_pt && end_qt ) { return NEGATIVE; }

    /*

#ifndef CGAL_DISABLE_AM_CODE
    // code added in previous version by Andreas + Monique -- start
    Site_2 const *pp1 = &p, *pp2 = &q;
    if ( !end_qt ) { std::swap(pp1, pp2); }

    if ( is_endpoint_of(*pp2, t) ) {
      Point_2 p1 = other_site(*pp2, t).point();
      Point_2 p2 = pp2->point();
      compute_vv(p, q, r, type);

      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p2, p1 - p2) );
    }
    // code added in previous version by Andreas + Monique -- end
#endif // CGAL_DISABLE_AM_CODE
    */

    // easy degeneracies --- end

    Sign retval = incircle_xxxs(p, q, r, t, type);

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s: about to return retval of"
      << "incircle_xxxs = " << retval << std::endl;);

    return retval;
  }

  //--------------------------------------------------------------------------
  // the first three objects are a point and two segments and the query
  // object is a segment
  //--------------------------------------------------------------------------

  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const PSS_Type& type) const
  {
    CGAL_precondition( p.is_point() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s PSS (pqrt) = "
      << "(" << p << ") (" << q << ") (" << r << ") "
      << "(" << t << ")" << std::endl;);

    // easy degeneracies --- start

    bool is_p_endp_of_q = is_endpoint_of(p, q);
    bool is_p_endp_of_r = is_endpoint_of(p, r);
    bool is_p_endp_of_t = is_endpoint_of(p, t);

    CGAL_SDG_DEBUG(std::cout
        << "debug incircle_s p is endp of q,r,t = "
        << is_p_endp_of_q << is_p_endp_of_r << is_p_endp_of_t
        << std::endl;);

    // check if p is a common endpoint of q and r, in which case the
    // Voronoi circle degenerates to p
    if ( is_p_endp_of_q && is_p_endp_of_r ) {
      CGAL_SDG_DEBUG(std::cout
          << "debug incircle_s voronoi circle degenerates to p="
          << p << std::endl;);
      // case 1: the new segment is not adjacent to the center of the
      //         degenerate Voronoi circle, i.e., not adjacent to p
      if ( ! is_p_endp_of_t ) { return POSITIVE; }

      // check if t has the same support as either q or r
      if ( same_segments(q.supporting_site(), t.supporting_site()) ) {
        return ZERO;
      }

      if ( same_segments(r.supporting_site(), t.supporting_site()) ) {
        return ZERO;
      }

      CGAL_SDG_DEBUG(std::cout
          << "debug incircle_s q, r, t have p as endpoint"
          << std::endl;);

      Point_2 r_ = r.source(), q_ = q.source(), t_ = t.source();

      if ( same_points(q.source_site(), p) ) { q_ = q.target(); }
      if ( same_points(r.source_site(), p) ) { r_ = r.target(); }
      if ( same_points(t.source_site(), p) ) { t_ =  t.target(); }

      Point_2 p_ = p.point();

      CGAL_SDG_DEBUG(std::cout << "debug incircle_s"
        << "  p_=" << p_ << "  q_= " << q_
        << "  r_=" << r_ << "  t_= " << t_
        << std::endl;);

      if ( CGAL::orientation(p_, q_, t_) == LEFT_TURN &&
           CGAL::orientation(p_, r_, t_) == RIGHT_TURN ) {
        return NEGATIVE;
      }
      return ZERO;
    }

    // here, voronoi circle does not degenerate to a point

    // philaris: remove
    /*
#ifndef CGAL_DISABLE_M_CODE
    // code added by Menelaos -- begin

    // in the code that follows we check whether one endpoint of the
    // query segment t is the same as the point p of a PSS circle. in
    // this case the result is known by taking the other point of t
    // and checking against the tangent to the Voronoi circle at p.
    if ( is_endpoint_of(p, t) ) {
      Point_2 p1 = p.point();
      Point_2 p2 = other_site(p, t).point();
      compute_vv(p, q, r, type);
      Compute_scalar_product_2 csp;
      return -CGAL::sign( csp(vv - p1, p2 - p1) );
    }
    // code added by Menelaos -- end
#endif // CGAL_DISABLE_M_CODE
    */

    // philaris: with assertions
    // check if t has the same support as either q or r
    if ( same_segments(q.supporting_site(), t.supporting_site()) ) {
      CGAL_SDG_DEBUG(std::cout
          << "debug q=" << q << " t=" << t << std::endl; );
      // philaris: the following assertion is too strong
      //CGAL_assertion(
      //    same_points(q.source_site(), t.source_site()) ||
      //    same_points(q.source_site(), t.target_site()) ||
      //    same_points(q.target_site(), t.source_site()) ||
      //    same_points(q.target_site(), t.target_site())   );
      // philaris: the following assertion is too strong and does not
      // work for inexact arithmetic;
      // philaris: I also remove it because of validity tests
      //CGAL_assertion(is_p_endp_of_q && is_p_endp_of_t);
      return POSITIVE;
    }
    if ( same_segments(r.supporting_site(), t.supporting_site()) ) {
      CGAL_SDG_DEBUG(std::cout
          << "debug r=" << r << " t=" << t
          << " have support " << r.supporting_site()  << std::endl; );

      // philaris: r and t share a point, which must be p

      // philaris: the assertion below fails in some cases,
      //           maybe because of inexact arithmetic
      // philaris: I also remove it because of validity tests
      //CGAL_assertion(is_p_endp_of_r && is_p_endp_of_t);

      // philaris: I also remove the following assertion
      //CGAL_assertion(
      //    same_points(r.source_site(), t.source_site()) ||
      //    same_points(r.source_site(), t.target_site()) ||
      //    same_points(r.target_site(), t.source_site()) ||
      //    same_points(r.target_site(), t.target_site())   );
      return POSITIVE;
    }

    if ( is_p_endp_of_q && is_p_endp_of_t ) {
      CGAL_SDG_DEBUG(std::cout << "debug incircle_s orientation"
        << " is_p_endp_of_q && is_p_endp_of_t"
        << std::endl;);
      Point_2 qother =
        same_points(p, q.source_site()) ? q.target() : q.source();
      Point_2 tother =
        same_points(p, t.source_site()) ? t.target() : t.source();
      Orientation o = CGAL::orientation(qother, p.point(), tother);

      CGAL_SDG_DEBUG(std::cout << "debug vsqrt "
          << "q^ = " << qother << ", p = " << p.point()
          << ", t^ = " << tother
          << " o = " << o << std::endl; );

      if (o != RIGHT_TURN) {
        return POSITIVE;
      } else {
        if (is_site_h_or_v(q) || is_site_h_or_v(t)) {
          return NEGATIVE;
        } else {
          bool has_q_pos_slope = has_positive_slope(q);
          if (has_q_pos_slope) {
            CGAL_SDG_DEBUG(std::cout << "debug vsqrt "
                << "q has positive slope" << std::endl; );

            return cmpy(qother, p.point()) == cmpy(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          } else {
            CGAL_SDG_DEBUG(std::cout << "debug vsqrt "
                << "q has has negative slope" << std::endl; );
            return cmpx(qother, p.point()) == cmpx(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          }
        }
      }
      //return (o == RIGHT_TURN)? NEGATIVE : POSITIVE;
    } // end of case (is_p_endp_of_q && is_p_endp_of_t) {

    if ( is_p_endp_of_r && is_p_endp_of_t ) {
      CGAL_SDG_DEBUG(std::cout << "debug incircle_s orientation"
        << " is_p_endp_of_r && is_p_endp_of_t"
        << std::endl;);
      Point_2 rother =
        same_points(p, r.source_site()) ? r.target() : r.source();
      Point_2 tother =
        same_points(p, t.source_site()) ? t.target() : t.source();
      Orientation o = CGAL::orientation(rother, p.point(), tother);

      if (o != LEFT_TURN) {
        return POSITIVE;
      } else {
        if (is_site_h_or_v(r) || is_site_h_or_v(t)) {
          return NEGATIVE;
        } else {
          bool has_r_pos_slope = has_positive_slope(r);
          if (has_r_pos_slope) {
            return cmpx(rother, p.point()) == cmpx(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          } else {
            return cmpy(rother, p.point()) == cmpy(tother, p.point()) ?
                   NEGATIVE : POSITIVE;
          }
        }
      }
    } // end of case (is_p_endp_of_r && is_p_endp_of_t)

    // easy degeneracies --- end

    return incircle_xxxs(p, q, r, t, type);
  }


  //--------------------------------------------------------------------------
  // the first three objects are segments and the query object is also
  // a segment
  //--------------------------------------------------------------------------

  inline
  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t, const SSS_Type& type) const
  {
    CGAL_precondition( p.is_segment() );
    CGAL_precondition( q.is_segment() );
    CGAL_precondition( r.is_segment() );
    CGAL_precondition( t.is_segment() );

    return incircle_xxxs(p, q, r, t, type);
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // the incircle test --- end
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------



  vertex_t
  compute_type(const Site_2& s1, const Site_2& s2, const Site_2& s3) const
  {
    int npts = 0;
    if ( s1.is_point() ) ++npts;
    if ( s2.is_point() ) ++npts;
    if ( s3.is_point() ) ++npts;

    switch ( npts ) {
    case 0:
      return SSS;
      break;
    case 1:
      return PSS;
      break;
    case 2:
      return PPS;
      break;
    default:
      return PPP;
    }
  }


  Sign incircle_p(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t) const
  {
    CGAL_precondition( t.is_point() );

    switch ( v_type ) {
    case PPP:
      return incircle_p(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
        return incircle_p(q, r, p, t, PPS_Type());
      } else if ( q.is_segment() ) {
        return incircle_p(r, p, q, t, PPS_Type());
      } else {
        return incircle_p(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
        return incircle_p(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
        return incircle_p(q, r, p, t, PSS_Type());
      } else {
        return incircle_p(r, p, q, t, PSS_Type());
      }
      return incircle_p(p, q, r, t, PSS_Type());
    default: // case SSS
      return incircle_p(p, q, r, t, SSS_Type());
    }
  }



  Sign incircle_s(const Site_2& p, const Site_2& q, const Site_2& r,
                  const Site_2& t) const
  {
    CGAL_precondition( t.is_segment() );

    CGAL_SDG_DEBUG(std::cout << "debug incircle_s (pqrt) = "
      << "(" << p << ") (" << q << ") (" << r << ") "
      << "(" << t << ")" << std::endl;);

    switch ( v_type ) {
    case PPP:
      return incircle_s(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
        return incircle_s(q, r, p, t, PPS_Type());
      } else if ( q_.is_segment() ) {
        return incircle_s(r, p, q, t, PPS_Type());
      } else {
        return incircle_s(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
        return incircle_s(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
        return incircle_s(q, r, p, t, PSS_Type());
      } else {
        return incircle_s(r, p, q, t, PSS_Type());
      }
    default: // case SSS
      return incircle_s(p, q, r, t, SSS_Type());
    }
  }


public:
  Voronoi_vertex_sqrt_field_new_C2(const Site_2& p,
                                   const Site_2& q,
                                   const Site_2& r)
    : p_(p), q_(q), r_(r), is_vv_computed(false)
  {
    v_type = compute_type(p, q, r);
  }


  inline bool is_degenerate_Voronoi_circle() const
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


  Point_2 degenerate_point() const
  {
    CGAL_precondition( is_degenerate_Voronoi_circle() );
    if ( p_.is_point() ) return p_.point();
    if ( q_.is_point() ) return q_.point();
    return r_.point();
  }


  Point_2 point() const
  {
    if ( is_degenerate_Voronoi_circle() ) {
      return degenerate_point();
    }

    if ( !is_vv_computed ) {
      switch ( v_type ) {
      case PPP:
        compute_vv(p_, q_, r_, PPP_Type());
        break;
      case PPS:
        if ( p_.is_segment() ) {
          compute_vv(q_, r_, p_, PPS_Type());
        } else if ( q_.is_segment() ) {
          compute_vv(r_, p_, q_, PPS_Type());
        } else {
          compute_vv(p_, q_, r_, PPS_Type());
        }
        break;
      case PSS:
        if ( p_.is_point() ) {
          compute_vv(p_, q_, r_, PSS_Type());
        } else if ( q_.is_point() ) {
          compute_vv(q_, r_, p_, PSS_Type());
        } else {
          compute_vv(r_, p_, q_, PSS_Type());
        }
        break;
      default: // case SSS:
        compute_vv(p_, q_, r_, SSS_Type());
        break;
      }
    }

    return vv;
  }


  inline Sign incircle(const Site_2& t) const
  {

    CGAL_SDG_DEBUG(std::cout
        << "debug field_new incircle t=" << t << std::endl;);

    if ( t.is_point() ) {
      return incircle_p(p_, q_, r_, t);
    }
    CGAL_SDG_DEBUG(std::cout
        << "debug about to run incircle_s (pqrt) ="
        << "(" << p_ << ") (" << q_ << ") (" << r_ << ") "
        << "(" << t << ")" << std::endl;);
    return incircle_s(p_, q_, r_, t);
  }


  inline Sign operator()(const Site_2& p, const Site_2& q, const Site_2& r,
                         const Site_2& t) const
  {
    if ( t.is_point() ) {
      return incircle_p(p, q, r, t);
    }
    return incircle_s(p, q, r, t);
  }

private:

  template<class Type>
  Sign incircle_p_no_easy(const Point_2& vv,
                          const Site_2& p, const Site_2& q, const Site_2& r,
                          const Site_2& t, const Type& type) const
  {
    CGAL_precondition( t.is_point() );

    FT radius = linf_radius(vv, p, q, r, type);

    Point_2 tt = t.point();

    FT absdvtx = CGAL::abs(vv.x() - tt.x());
    FT absdvty = CGAL::abs(vv.y() - tt.y());

    FT d = (CGAL::max)(absdvtx, absdvty);

    Comparison_result crude = CGAL::compare(d, radius);

    if (crude != ZERO) {
      return crude;
    } else {
      CGAL_SDG_DEBUG(std::cout << "debug vsqr refining in noeasy pqr=("
        << p << ", " << q << ", " << r << "), "
        << "t=" << t
        << std::endl;);
      // here crude == ZERO, so
      // we might have to refine

      FT radius_fine = linf_fine_radius(vv, p, q, r, type);

      FT d_fine = (CGAL::min)(absdvtx, absdvty);

      return CGAL::compare(d_fine, radius_fine);

    }
  }


  Sign incircle_p_no_easy(const Site_2& p, const Site_2& q, const Site_2& r,
                          const Site_2& t) const
  {
    Sign s(ZERO);
    switch ( v_type ) {
    case PPP:
      s = incircle_p(p, q, r, t, PPP_Type());
      break;
    case PPS:
      PPS_Type pps;
      if ( p.is_segment() ) {
        compute_vv(q, r, p, pps);
        s = incircle_p_no_easy(vv, q, r, p, t, pps);
      } else if ( q.is_segment() ) {
        compute_vv(r, p, q, pps);
        s = incircle_p_no_easy(vv, r, p, q, t, pps);
      } else {
        compute_vv(p, q, r, pps);
        s = incircle_p_no_easy(vv, p, q, r, t, pps);
      }
      break;
    case PSS:
      PSS_Type pss;
      if ( p.is_point() ) {
        compute_vv(p, q, r, pss);
        s = incircle_p_no_easy(vv, p, q, r, t, pss);
      } else if ( q.is_point() ) {
        compute_vv(q, r, p, pss);
        s = incircle_p_no_easy(vv, q, r, p, t, pss);
      } else {
        compute_vv(r, p, q, pss);
        s = incircle_p_no_easy(vv, r, p, q, t, pss);
      }
      break;
    case SSS:
      SSS_Type sss;
      compute_vv(p, q, r, sss);
      s = incircle_p_no_easy(vv, p, q, r, t, sss);
      break;
    }

    return s;
  }


  Sign incircle_s_no_easy(const Site_2& p, const Site_2& q, const Site_2& r,
                          const Site_2& t) const
  {
    switch ( v_type ) {
    case PPP:
      return incircle_xxxs(p, q, r, t, PPP_Type());
    case PPS:
      if ( p.is_segment() ) {
        return incircle_xxxs(q, r, p, t, PPS_Type());
      } else if ( q_.is_segment() ) {
        return incircle_xxxs(r, p, q, t, PPS_Type());
      } else {
        return incircle_xxxs(p, q, r, t, PPS_Type());
      }
    case PSS:
      if ( p.is_point() ) {
        return incircle_xxxs(p, q, r, t, PSS_Type());
      } else if ( q.is_point() ) {
        return incircle_xxxs(q, r, p, t, PSS_Type());
      } else {
        return incircle_xxxs(r, p, q, t, PSS_Type());
      }
    default: // case SSS:
      return incircle_xxxs(p, q, r, t, SSS_Type());
    }
  }




public:
  inline Sign incircle_no_easy(const Site_2& t) const
  {
    Sign s;

    if ( t.is_point() ) {
      s = incircle_p_no_easy(p_, q_, r_, t);
    } else {
      s = incircle_s_no_easy(p_, q_, r_, t);
    }

    return s;
  }


  Orientation orientation(const Line_2& l) const
  {
    switch ( v_type ) {
    case PPP:
      compute_vv(p_, q_, r_, PPP_Type());
      break;
    case PPS:
      if ( p_.is_segment() ) {
        compute_vv(q_, r_, p_, PPS_Type());
      } else if ( q_.is_segment() ) {
        compute_vv(r_, p_, q_, PPS_Type());
      } else {
        compute_vv(p_, q_, r_, PPS_Type());
      }
      break;
    case PSS:
      if ( p_.is_point() ) {
        compute_vv(p_, q_, r_, PSS_Type());
      } else if ( q_.is_point() ) {
        compute_vv(q_, r_, p_, PSS_Type());
      } else {
        compute_vv(r_, p_, q_, PSS_Type());
      }
      break;
    default: // case SSS:
      compute_vv(p_, q_, r_, SSS_Type());
      break;
    }

    return CGAL::sign(l.a() * vv.x() + l.b() * vv.y() + l.c());
  }


  inline Oriented_side oriented_side(const Line_2& l) const
  {
    return orientation(l);
  }

private:
  // the defining sites of the Voronoi vertex
  const Site_2& p_, &q_, &r_;

  // indicates whether the Voronoi vertex has been computed
  mutable bool is_vv_computed;

  // the type of the Voronoi vertex
  vertex_t v_type;

  // the computed Voronoi vertex is cached in this variable
  mutable Point_2 vv;
};





} //namespace SegmentDelaunayGraphLinf_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_VORONOI_VERTEX_SQRT_FIELD_NEW_C2_H
