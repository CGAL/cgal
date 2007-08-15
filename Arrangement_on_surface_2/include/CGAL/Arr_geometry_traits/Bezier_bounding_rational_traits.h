// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
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
// Author(s)     : Iddo Hanniel <iddoh@cs.technion.ac.il>

#ifndef CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H
#define CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H

/*! \file
 * Definition of the Bezier_bounding_rational_traits<Kernel> class.
 */

#include <CGAL/functional.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h> 
#include <CGAL/Arr_geometry_traits/de_Casteljau_2.h> 

#include <deque>
#include <list>
#include <math.h>

CGAL_BEGIN_NAMESPACE

template <class _NT> 
struct _Bez_point_bound ;

template <class _NT> 
struct _Bez_point_bbox ;

/*!\ class Bezier_bounding_rational_traits
 * A traits class for performing bounding operations on Bezier curves,
 * assuming the NT is rational.
 */
template <typename _Kernel>
class Bezier_bounding_rational_traits
{
public:
  typedef _Kernel Kernel; // Used in polygon and point/vector operations
  typedef typename Kernel::FT NT; // The number type used to represent control points and t-parameters.
  typedef typename Kernel::Point_2 Point_2; // The control point type.
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef std::deque<Point_2>                   Control_point_vec;

  typedef  _Bez_point_bound<NT>                 Bez_point_bound;
  typedef  _Bez_point_bbox<NT>                  Bez_point_bbox;

  // Auxilary structure for output of intersection points.
  struct Bound_pair
  {
    Bez_point_bound bound1;
    Bez_point_bound bound2;
    Bez_point_bbox bbox;
  
    Bound_pair(const Bez_point_bound& b1, const Bez_point_bound& b2,
      const Bez_point_bbox& bb) : bound1(b1), bound2(b2), bbox(bb) {}
    Bound_pair() : bound1(), bound2(), bbox() {}
  };

  typedef std::list<Bound_pair> Bound_pair_lst;

  //TODO - add here typedefs for vertical tangency output (pair<..>).

  Kernel  kernel;
  NT      _can_refine_bound;  

  // a ctr that enables to change can_refine_bound.
  // iddo: important! we need to make the arr traits receive this traits as a reference,
  // since it now contains state.
  // The line below currently causes a runtime bug, so for now.
  //Bezier_bounding_rational_traits(double can_refine_bound = std::pow(2,-53))
                                  //0.00000000000000000000001)
  Bezier_bounding_rational_traits (double can_refine_bound =
                                   0.00000000000000011)
    : _can_refine_bound(can_refine_bound)
  {}
  
  // iddo: can_refine() function can be based on Euclidean or parametric space.
  // For rationals it is a simple parametric-space check.
  bool can_refine(const Control_point_vec& , const NT& l, const NT& r)
  {
    // iddo: might be made more efficient if needed (based on denom...)
    if (r-l < _can_refine_bound)
      return false;
    return true;  
  }

  // Can be made a template function <Control_point_vec, NT>
  /*
  void DeCasteljau(const Control_point_vec& cp, const NT& t0,
    Control_point_vec& lcp, Control_point_vec& rcp) const
  {
    const NT t1 = NT(1) - t0; 

    lcp.push_back(cp.front());
    rcp.push_front(cp.back());

    Control_point_vec aux(cp.begin(), cp.end());
    unsigned int i=0;
    for (; i<cp.size()-1; ++i)
    {
      Control_point_vec aux1(cp.size()-i-1); //in every step of the subdivision we reduce by 1 the size of the auxilary vector
      unsigned int j=0;
      for (; j<cp.size()-i-1; ++j)
      {
        aux1[j] = Point_2(t1*(aux[j].x()) + t0*(aux[j+1].x()), t1*(aux[j].y()) + t0*(aux[j+1].y())); 
      }
      lcp.push_back(aux1.front());
      rcp.push_front(aux1.back());

      aux = aux1;
    }
  }
  */

  ////////////////////////////////////////////////////////////////
  //Intersection points
  ///////////////////////////////////////////////////////////////
  void CrvCrvInter(const Control_point_vec& cp1, const Control_point_vec& cp2, 
                  Bound_pair_lst& intersection_pairs)
  {
    //iddo: in order to support interval, we will need to convert Rational to
    //      interval do we need to do it outside the call to the function or
    //      here (needs conversion from Rat_Point...).
    
    std::set<NT> endpoint_intersections; // Auxilary variable needed in the recursion
    // to identify already existing points due to endpoint intersections.

    // Call to recursive function.
    _CrvCrvInterAux (cp1, 0, 1,
                     cp2, 0, 1,
                     true,
                     endpoint_intersections,
                     intersection_pairs);
  }

  void refine_intersection_point (const Bound_pair& intersection_point,
                                  Bound_pair& refined_point)
  {
    // Assumes we do not need to check span, intersection_point has all the
    //necessary data.
    CGAL_precondition(intersection_point.bound1.can_refine);
    CGAL_precondition(intersection_point.bound2.can_refine);
    // The following precondition makes sure the point is not rational already
    // (in which case it would not need refinement).
//    CGAL_precondition(intersection_point.bbox.min_x != intersection_point.bbox.max_x);

    const Control_point_vec& cv1 = intersection_point.bound1.bounding_polyline;
    const NT& l1 = intersection_point.bound1.t_min;
    const NT& r1 = intersection_point.bound1.t_max;
    const Control_point_vec& cv2 = intersection_point.bound2.bounding_polyline;
    const NT& l2 = intersection_point.bound2.t_min;
    const NT& r2 = intersection_point.bound2.t_max;

    /* Subdivide the two curves and recurse. */
    // This is probably not needed as it will be handled
    // in the recursion.
    bool can_refine1 = can_refine(cv1, l1, r1);
    bool can_refine2 = can_refine(cv2, l2, r2);
    if (!can_refine1 || !can_refine2) 
    {
      // Failed in the bounded approach - stop the subdivision and indicate
      // that the subdivision has failed.
      refined_point = intersection_point;
      refined_point.bound1.can_refine = false;
      refined_point.bound2.can_refine = false;

      return;
    }

    Control_point_vec cv1a, cv1b;
    
    bisect_control_polygon_2 (cv1.begin(), cv1.end(),
                              std::back_inserter(cv1a),
                              std::front_inserter(cv1b));
    //DeCasteljau(cv1, 0.5, cv1a, cv1b);
    NT t_mid1 = NT(0.5) * (l1 + r1);

    Control_point_vec cv2a, cv2b; 

    bisect_control_polygon_2 (cv2.begin(), cv2.end(),
                              std::back_inserter(cv2a),
                              std::front_inserter(cv2b));
    //DeCasteljau(cv2, 0.5, cv2a, cv2b);
    NT t_mid2 = NT(0.5) * (l2 + r2);

    // Iddo: Catch situations that we cannot refine, in roder to prevent
    //       having multiple points as the refined point.
    bool can_refine_pt_cv1a = can_refine (cv1a, l1, t_mid1);
    bool can_refine_pt_cv1b = can_refine (cv1b, t_mid1, r1);
    bool can_refine_pt_cv2a = can_refine (cv2a, l2, t_mid2);
    bool can_refine_pt_cv2b = can_refine (cv2b, t_mid2, r2);
    
    if (!can_refine_pt_cv1a || !can_refine_pt_cv1b ||
        !can_refine_pt_cv2a || !can_refine_pt_cv2b)
    {
      // Construct a Bez_point_bound that is inconclusive:
      // all parameters found so far + can_refine = false, with the bounding
      // box of cv1, and add to output.
      Bez_point_bbox bbox1;

      cp_bbox (cv1, bbox1);
      refined_point = 
        Bound_pair (Bez_point_bound (cv1, l1, r1,
                                     Bez_point_bound::INTERSECTION_PT, false),
                    Bez_point_bound (cv2, l2, r2,
                                     Bez_point_bound::INTERSECTION_PT, false),
                    bbox1);

      return;
    }

    //Recursion:
    Bound_pair_lst intersection_pairs;
    std::set<NT> endpoint_intersections;
    bool should_check_span = false; //this precondition is hard to check efficiently.

    _CrvCrvInterAux (cv1a, l1, t_mid1, cv2a, l2, t_mid2,
                     should_check_span, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1a, l1, t_mid1, cv2b, t_mid2, r2,
                     should_check_span, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1b, t_mid1, r1, cv2a, l2, t_mid2,
                     should_check_span, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1b, t_mid1, r1, cv2b, t_mid2, r2,
                     should_check_span, endpoint_intersections,
                     intersection_pairs);

    refined_point = intersection_pairs.front(); 

    CGAL_assertion (intersection_pairs.size() == 1);

    return;
  }

  ////////////////////////////////////////////////////////////////
  // Vertical tangency points
  ///////////////////////////////////////////////////////////////

  // Vertical tangency points (cp is a Bez curve between t0-t1).
  // Output are pairs of <bound,bbox> representing the vertical tangency points.
  void vertical_tangency_points
      (const Control_point_vec& cp,
       const NT& t0, const NT& t1,
       std::list<std::pair<Bez_point_bound, Bez_point_bbox> >& tangency_points)
  {
    //TODO - handle special case of degree two curves.

    bool can_refine_pt = can_refine(cp, t0, t1);
    if (!can_refine_pt)
    {
      // Construct a Bez_point_bound that is inconclusive:
      // all parameters found so far + can_refine = false,
      // and add to output.
      Bez_point_bbox bbox;
      cp_bbox(cp, bbox);
      tangency_points.push_back
        (std::make_pair(Bez_point_bound(cp, t0, t1,
                                        Bez_point_bound::VERTICAL_TANGENCY_PT,
                                        false),
                        bbox));
      return;
    }

    // Check for x-monoticity (if exists return)
    bool is_x_monotone = _is_x_monotone(cp);
    if (is_x_monotone)
      return;

    // Check for y-monoticity (else recurse)
    bool is_y_monotone = _is_y_monotone(cp);
    if (!is_y_monotone)
    {
      // Recurse
      Control_point_vec left, right;

      bisect_control_polygon_2 (cp.begin(), cp.end(),
                                std::back_inserter(left),
                                std::front_inserter(right));
      //DeCasteljau(cp, 0.5, left, right);
      NT t_mid = NT(0.5) * (t0 + t1);
      // TODO - handle the case where t_mid is a vertical (rational) tangency
      // point, by checking whether the first vector of right is vertical.
      // Ron: This happens in the Bezier.dat input file!
      CGAL_assertion(CGAL::compare(right[0].x(), right[1].x()) != EQUAL);

      vertical_tangency_points(left, t0, t_mid, tangency_points);
      vertical_tangency_points(right, t_mid, t1, tangency_points);
      return;
    }

    // Check for convexity (since it is y-monotone and not x-monotone)
    // if convex - add to list, else recurse.
    if (is_convex_2 (cp.begin(), cp.end(), kernel))
    {
      Bez_point_bbox bbox;
      typename Control_point_vec::const_iterator res;
      res = bottom_vertex_2 (cp.begin(), cp.end(), kernel);
      bbox.min_y = res -> y();
      res = top_vertex_2 (cp.begin(), cp.end(), kernel);
      bbox.max_y = res -> y();
      res = left_vertex_2 (cp.begin(), cp.end(), kernel);
      bbox.min_x = res -> x();
      res = right_vertex_2 (cp.begin(), cp.end(), kernel);
      bbox.max_x = res -> x();

      tangency_points.push_back(std::make_pair(
        Bez_point_bound(cp, t0, t1, Bez_point_bound::VERTICAL_TANGENCY_PT, true),
        bbox));
    }
    else
    {
      // Recurse
      Control_point_vec left, right;

      bisect_control_polygon_2 (cp.begin(), cp.end(),
                                std::back_inserter(left),
                                std::front_inserter(right));
      //DeCasteljau(cp, 0.5, left, right);
      NT t_mid = NT(0.5) * (t0 + t1);

      // TODO - handle the case where t_mid is a vertical (rational) tangency
      // point, by checking whether the first vector of right is vertical.
      CGAL_assertion(CGAL::compare(right[0].x(), right[1].x()) != EQUAL);

      vertical_tangency_points(left, t0, t_mid, tangency_points);
      vertical_tangency_points(right, t_mid, t1, tangency_points);
      return;
    }

    // Debug print
    //std::cout << "number of tangency points is " << tangency_points.size() << std::endl;
  }

  // Refinement of vertical tangency point.
  // Assumes the tangency point is isolated (only a single point in the interval).
  void refine_tangency_point(const Control_point_vec& cp, const NT& t0, const NT& t1,
    std::pair<Bez_point_bound, Bez_point_bbox>& tangency_point)
  {
    // TODO - a better implmentation that just checks which of the two subdivided
    // curves is x-monotone (since they should both be y-monotone and convex).

    Control_point_vec left, right;
    /*
    // iddo: the can_refine check is probably not needed here (the check
    // in the recursion will handle it).
    bool can_refine_pt = can_refine(cp, left, right);
    if (!can_refine)
    {
      CGAL_assertion(false); 
    }
    */

    bisect_control_polygon_2 (cp.begin(), cp.end(),
                              std::back_inserter(left),
                              std::front_inserter(right));
    //DeCasteljau(cp, 0.5, left, right);
    NT t_mid = NT(0.5) * (t0 + t1);

    // TODO - handle the case where t_mid is a vertical (rational) tangency point,
    // by checking whether the first vector of right is vertical.
    CGAL_assertion(CGAL::compare(right[0].x(), right[1].x()) != EQUAL);

    bool can_refine_pt_left = can_refine(left, t0, t_mid);
    bool can_refine_pt_right = can_refine(right, t_mid, t1);

    if (!can_refine_pt_left || !can_refine_pt_right)
    {
      // Construct a Bez_point_bound that is inconclusive:
      // all parameters found so far + can_refine = false,
      // and add to output.
      Bez_point_bbox bbox;
      cp_bbox(cp, bbox);
      tangency_point = 
        std::make_pair(Bez_point_bound(cp, t0, t1,
                                       Bez_point_bound::VERTICAL_TANGENCY_PT,
                                       false),
                       bbox);
      return;
    }

    std::list<std::pair<Bez_point_bound, Bez_point_bbox> > aux1, aux2;
    vertical_tangency_points(left, t0, t_mid, aux1);
    vertical_tangency_points(right, t_mid, t1, aux2);

    CGAL_assertion(aux1.size() + aux2.size() == 1);
    if (aux1.size() == 1)
    {
      tangency_point = *(aux1.begin());
    }
    else 
    {
      tangency_point = *(aux2.begin());
    }
  }

  ////////////////////////////////////////////////////////////////
  // Slope comparison.
  ///////////////////////////////////////////////////////////////
  Comparison_result compare_slopes_of_bounded_intersection_point(
    const Bez_point_bound& bound1, const Bez_point_bound& bound2)
  {

    // Since we assume the spans do not overlap, we can just compare any vector
    // of the hodograph spans.
    const Control_point_vec& cv1 = bound1.bounding_polyline;
    const Control_point_vec& cv2 = bound2.bounding_polyline;

    // An (expensive) check for the precondition that the spans do not overlap. 
    // Therefor it is commented out.
    /*
    Vector_2 dir1, dir2, leftMost1, leftMost2, rightMost1, rightMost2;
    bool spansOverlap;
    bool cv1SpanOK = _CrvTanAngularSpan(cv1, dir1, leftMost1, rightMost1);
	  bool cv2SpanOK = _CrvTanAngularSpan(cv2, dir2, leftMost2, rightMost2);
    if (cv1SpanOK && cv2SpanOK)
    {
  	  spansOverlap = _AngularSpansOverlap(leftMost1, rightMost1, leftMost2, rightMost2);
    }
    else
    {
      spansOverlap = true;
    }
    CGAL_precondition(!spansOverlap);
    */

    Vector_2 dir1 = cv1.back() - cv1.front();
    Vector_2 dir2 = cv2.back() - cv2.front();

    // RI: use compare_angle_with_x_axis (after doing: if (dir.x()<0) dir=-dir)
    // instead of division.
    NT slope1 = dir1.y()/dir1.x();
    NT slope2 = dir2.y()/dir2.x();
    return CGAL::compare(slope1, slope2);
  }


  // Can be made a template function <Control_point_vec, Kernel>
  void cp_bbox(const Control_point_vec& cp, Bez_point_bbox& bez_bbox) const
  {
   CGAL_assertion (! cp.empty());

    typename Kernel::Compare_x_2       comp_x = kernel.compare_x_2_object();
    typename Kernel::Compare_y_2       comp_y = kernel.compare_y_2_object();

    typename Control_point_vec::const_iterator it = cp.begin();
    typename Control_point_vec::const_iterator min_x = it;
    typename Control_point_vec::const_iterator max_x = it;
    typename Control_point_vec::const_iterator min_y = it;
    typename Control_point_vec::const_iterator max_y = it;

    while (it != cp.end())
    {
      if (comp_x (*it, *min_x) == SMALLER)
        min_x = it;
      else if (comp_x (*it, *max_x) == LARGER)
        max_x = it;

      if (comp_y (*it, *min_y) == SMALLER)
        min_y = it;
      else if (comp_y (*it, *max_y) == LARGER)
        max_y = it;

      ++it;
    }

    bez_bbox.min_x = min_x->x(); 
    bez_bbox.max_x = max_x->x(); 
    bez_bbox.min_y = min_y->y(); 
    bez_bbox.max_y = max_y->y(); 
  }

private:

  /*
  void _construct_bboxes (const Control_point_vec& cv1,
                          const Control_point_vec& cv2,
                          Bez_point_bbox& bbox1, Bez_point_bbox& bbox2) const
  {
    const double   diff_bound = 0.1875; // radians, a bit more than 10 degrees.
    const double   _pi = 3.14159265;

    // Compute the angles theta1 and theta2 the lines that connect the first
    // and last control points in each curve form with the x-axis.
    const double   dx1 = CGAL::to_double (cv1.back().x()) -
                         CGAL::to_double (cv1.front().x());
    const double   dy1 = CGAL::to_double (cv1.back().y()) -
                         CGAL::to_double (cv1.front().y());
    const double   theta1 = atan2 (dy1, dx1);

    const double   dx2 = CGAL::to_double (cv2.back().x()) -
                         CGAL::to_double (cv2.front().x());
    const double   dy2 = CGAL::to_double (cv2.back().y()) -
                         CGAL::to_double (cv2.front().y());
    const double   theta2 = atan2 (dy2, dx2);

    // Check if the two angles are close, and compute half their average.
    const double   diff_theta = fabs (theta2 - theta1);

    double         half_theta;

    if (diff_theta < diff_bound)
    {
      // The two angles are close:
      half_theta = (theta1 + theta2) / 4;
    }
    else if (fabs (diff_theta - _pi) < diff_bound)
    {
      // We have theta2 ~= theta1 +/- PI, thus we take:
      if (theta2 > 0)
        half_theta = (theta1 + theta2 - _pi) / 4;
      else
        half_theta = (theta1 + theta2 + _pi) / 4;
    }
    else
    {
      // In this case the angles are not close, so we compute their
      // axis-alligned bounding boxes.
      cp_bbox(cv1, bbox1);
      cp_bbox(cv2, bbox2);
      return;
    }

    // Check if theta/2 is close to 0 or to +/- PI/2.
    if (fabs(half_theta) < 0.0625 || _pi/2 - fabs(half_theta) < 0.0625)
    {
      // In this case, using the axis-alligned bounding box is fine.
      cp_bbox(cv1, bbox1);
      cp_bbox(cv2, bbox2);
      return;
    }

    // Compute a rational tau ~= tan(theta/2) and use it to find a rotation
    // angle close to theta with rational sine and cosine.
    const double   tan_half_theta = tan (half_theta);
    const int      denom = 1000;
    const int      numer = static_cast<int> (denom * tan_half_theta + 0.5);
    const NT       tau (numer, denom);
    const NT       sqr_tau = tau*tau;
    const NT       sin_theta = 2*tau / (1 + sqr_tau);
    const NT       cos_theta = (1 - sqr_tau) / (1 + sqr_tau);

    // Rotate both control polygons by the minus the approximated angle. 
    Control_point_vec                          rot_cv1;
    Control_point_vec                          rot_cv2;
    typename Control_point_vec::const_iterator it;

    for (it = cv1.begin(); it != cv1.end(); ++it)
    {
      rot_cv1.push_back (Point_2 (cos_theta * it->x() + sin_theta * it->y(),
                                  -sin_theta * it->x() + cos_theta * it->y()));
    }

    for (it = cv2.begin(); it != cv2.end(); ++it)
    {
      rot_cv2.push_back (Point_2 (cos_theta * it->x() + sin_theta * it->y(),
                                  -sin_theta * it->x() + cos_theta * it->y()));
    }

    // Compute the axis-alligned bounding boxes of the rotated curves, which
    // are equivalent of the rotated skewed bounding boxes of the curves.
    cp_bbox (rot_cv1, bbox1);
    cp_bbox (rot_cv2, bbox2);

    return;
  }
  */

  // Can be made a template function <Control_point_vec, Kernel>
  bool _is_monotone(const Control_point_vec& cp, bool check_x) const
  {
    CGAL_precondition(cp.size() > 1);
    Comparison_result                          curr_res;
    Comparison_result                          res = EQUAL;

    typename Control_point_vec::const_iterator cp_iter = cp.begin();
    typename Control_point_vec::const_iterator cp_iter_end = cp.end();
    typename Control_point_vec::const_iterator cp_iter_next = cp_iter; 
    ++cp_iter_next;

    while (cp_iter_next != cp_iter_end)
    {
      curr_res = check_x ? compare_x(*cp_iter, *cp_iter_next) :
                           compare_y(*cp_iter, *cp_iter_next);
      if (curr_res != EQUAL)
      {
        res = curr_res;
        ++cp_iter; 
        ++cp_iter_next;
        break;
      }
      ++cp_iter;
      ++cp_iter_next;
    }
    if (cp_iter_next == cp_iter_end)
      // iso-parallel segment.
      // In our context, we consider a horizontal segment as non-y-monotone
      // but a vertical segment is considered x-monotone (is this right?).
    {
      if (check_x)
        return true; // Vertical segment when checking x-monoticity
      return false;  // Horizontal segment when checking y-monoticity
    }

    while (cp_iter_next != cp_iter_end)
    {
      curr_res = check_x ? compare_x(*cp_iter, *cp_iter_next) :
                           compare_y(*cp_iter, *cp_iter_next);
      if (curr_res == EQUAL)
      {
        ++cp_iter;
        ++cp_iter_next;
        continue;
      }
      if (res != curr_res)
        break;

      ++cp_iter;
      ++cp_iter_next;
    }
    if (cp_iter_next != cp_iter_end)
      return false;
    return true;
  }

  inline bool _is_x_monotone(const Control_point_vec& cp) const
  {
    return _is_monotone(cp, true);
  }

  inline bool _is_y_monotone(const Control_point_vec& cp) const
  {
    return _is_monotone(cp, false);
  }


  // Can be made a template function <Control_point_vec, Kernel>
  // Returns true if the span is less than 90 degrees (we don't want to work with larger spans).
  bool _CrvTanAngularSpan(const Control_point_vec& cp,
				Vector_2& dir,
				Vector_2& leftMost,
				Vector_2& rightMost)
  {
    Vector_2 v;
    const Point_2& O(ORIGIN);

    dir = cp.back() - cp.front();

	  // Initialize LeftMost and RightMost
	  leftMost = dir;
	  rightMost = dir;

    typename Control_point_vec::const_iterator iterCP0 = cp.begin();
    typename Control_point_vec::const_iterator iterCP1 = iterCP0;
    ++iterCP1;
    for (; iterCP1!=cp.end(); ++iterCP1, ++iterCP0) 
    {
      v = (*iterCP1) - (*iterCP0);

      // Test for 90 degree angle, we can loosen it to test for 180 degree between
	    // left and right but we don't really want to work with larger spans.
      // TODO - try get rid of these ORIGINs
      if (angle(ORIGIN+v, O, ORIGIN+dir) != ACUTE)
	      return false;

	    if (left_turn(ORIGIN+v, O, ORIGIN+dir)) 
      {
		    if (left_turn(ORIGIN+v, O, ORIGIN+leftMost)) 
        {
			    leftMost = v;
		    }
	    }
	    else 
      {
		    if (right_turn(ORIGIN+v, O, ORIGIN+rightMost)) 
        {
			    rightMost = v;		
		    }
	    }
    } //for

    return true;
  }

  // TODO - there is already a kernel function that does this ccw?
  // Can be made a template function <Kernel>
  inline bool _VecIsBetween(const Vector_2& vec,
							      const Vector_2& leftVec,
							      const Vector_2& rightVec)
  {
    const Point_2& O(ORIGIN);
	  return 	(right_turn(ORIGIN+rightVec, O, ORIGIN+vec) && left_turn(ORIGIN+leftVec, O, ORIGIN+vec));
  }

  // Can be made a template function <Kernel>
  bool _AngularSpansOverlap(
				  const Vector_2& leftMost1,
				  const Vector_2& rightMost1,
				  const Vector_2& leftMost2,
				  const Vector_2& rightMost2
				  )
  {
	  Vector_2 negVec;

	  // Check if any of span2 vectors is between span1.
	  if (_VecIsBetween(leftMost2, leftMost1, rightMost1))
		  return true;

	  negVec = -leftMost2;
	  if (_VecIsBetween(negVec, leftMost1, rightMost1))
		  return true;

	  if (_VecIsBetween(rightMost2, leftMost1, rightMost1))
		  return true;

	  negVec = -rightMost2;
	  if (_VecIsBetween(negVec, leftMost1, rightMost1))
		  return true;


	  // Note: If we got here, the only possibility for intersection
	  // is that span1 is totally between span2

	  // Check if span1 is between span2
	  if (_VecIsBetween(leftMost1, leftMost2, rightMost2))
		  return true;

	  negVec = -leftMost1;
	  if (_VecIsBetween(negVec, leftMost2, rightMost2))
		  return true;


	  // The rest of the symmetric conditions below are redundant (see note above)
	  return false;

    /*
	  if (_VecIsBetween(rightMost1, leftMost2, rightMost2))
		  return true;

	  negVec = -rightMost1;
	  if (_VecIsBetween(negVec, leftMost2, rightMost2))
		  return true;

    return false;
    */			
  }

  // Constructing the curve skewed bbox - actually returning
  // two parallel lines (parallel to the first-last line) that bound
  // the curve from above and below.
  // Using less_signed_distance_to_line_2 enables not to construct
  // new points (only take maximum and minimum signed dist points
  // and add the first-last vector).
  // Can be made a template function <Control_point_vec, Kernel>
  void _cvSkewBbox(const Control_point_vec& cp,
				Line_2& la,
				Line_2& lb)
  {
    // Functor for comparing two points signed distance from a given line (used for min/max_elem).
    // Binder taken from <CGAL/functional.h>, if this doesn't compile, construct your own functor.
    typedef CGALi::Binder<typename Kernel::Less_signed_distance_to_line_2, Arity_tag< 3 >, Line_2, 1 > Less_signed_dist_to_line;
    // Explanation: Take the 3-ary kernel functor Less_signed_distance_to_line_2 and bind its
    // first parameter, the line, to be our given line l.

    Line_2 l(cp.front(), cp.back());
    Vector_2 v(cp.front(), cp.back());
    typename Kernel::Less_signed_distance_to_line_2 lsdtl2 = kernel.less_signed_distance_to_line_2_object();

    typename Control_point_vec::const_iterator aux;
    aux = std::min_element(cp.begin(), cp.end(), Less_signed_dist_to_line(lsdtl2, l));
    la = Line_2(*aux, v);

    aux = std::max_element(cp.begin(), cp.end(), Less_signed_dist_to_line(lsdtl2, l));
    lb = Line_2(*aux, v);
  }


  // Returns true if endpoints coincide, and adds the Bound_pair
  // to the list if the endpoint hasn't already been inserted there.
  // Precondition !spansOverlap
  bool _endpoints_coincide(const Control_point_vec& cv1, const NT& l1, const NT& r1,
			        const Control_point_vec& cv2, const NT& l2, const NT& r2,
              std::set<NT>& endpoint_intersections,
              Bound_pair_lst& intersection_pairs) const
  {
    const Point_2& s1 = cv1.front();
    const Point_2& t1 = cv1.back();
    const Point_2& s2 = cv2.front();
    const Point_2& t2 = cv2.back();

    // What happens when an endpoint of one is an inner point (e.g., 1/3) of another?
    // We will not be able to discover this endcase and therefore get to can_refine()==false.
    int endpoint_intersection_num = 0; // For sanity check and return value.

    NT x, y, t_param1, t_param2;
    if (s1==s2) // Using operator== of kernel Point_2
    {
      ++endpoint_intersection_num;
      x = s1.x();
      y = s1.y();
      t_param1 = l1;
      t_param2 = l2;
    }
    if (s1==t2)
    {
      ++endpoint_intersection_num;
      x = s1.x();
      y = s1.y();
      t_param1 = l1;
      t_param2 = r2;
    }
    if (t1==s2)
    {
      ++endpoint_intersection_num;
      x = t1.x();
      y = t1.y();
      t_param1 = r1;
      t_param2 = l2;
    }
    if (t1==t2)
    {
      ++endpoint_intersection_num;
      x = t1.x();
      y = t1.y();
      t_param1 = r1;
      t_param2 = r2;
    }

    CGAL_assertion(endpoint_intersection_num <= 1);
    if (endpoint_intersection_num == 0)
      return false;

    // Construct a degenerate rational Bound_pair
    if (endpoint_intersections.find(t_param1) == endpoint_intersections.end())
    {
      // Found a new intersection point that wasn't inserted yet.
      endpoint_intersections.insert(t_param1); 
      intersection_pairs.push_back(Bound_pair(
                                    Bez_point_bound(cv1, t_param1, t_param1, Bez_point_bound::RATIONAL_PT, true),
                                    Bez_point_bound(cv2, t_param2, t_param2, Bez_point_bound::RATIONAL_PT, true),
                                    Bez_point_bbox(x,x,y,y)
                                  )
                                );
    }
    return true;
  }


  // Auxilary recursive function for intersecting two curves.
  void _CrvCrvInterAux
      (const Control_point_vec& cv1, const NT& l1, const NT& r1,
       const Control_point_vec& cv2, const NT& l2, const NT& r2,
       bool should_check_span,
       std::set<NT>& endpoint_intersections,
       Bound_pair_lst& intersection_pairs)
  {
    //iddo debug
    //static unsigned int dbg = 0;
    //std::cout << ++dbg << ", " << recursion_level << "; ";

    // Check if we got to subdivision termination criteria. 
    bool can_refine1 = can_refine(cv1, l1, r1);
    bool can_refine2 = can_refine(cv2, l2, r2);;
    if (!can_refine1 || !can_refine2) 
    {
      //Failed in the bounded approach - stop the subdivision.
      //Construct output with can_refine=false...  We don't really need a
      // meaningful bbox since we can't refine, so we take cv1 bbox (which
      // certainly contains the point).
      Bez_point_bbox bbox1;

      cp_bbox(cv1, bbox1);
      intersection_pairs.push_back
        (Bound_pair (Bez_point_bound (cv1, l1, r1,
                                      Bez_point_bound::INTERSECTION_PT, false),
                     Bez_point_bound (cv2, l2, r2,
                                      Bez_point_bound::INTERSECTION_PT, false),
                     bbox1));

      return;
    }

    Bez_point_bbox bbox1, bbox2;

    //_construct_bboxes (cv1, cv2,
    //                   bbox1, bbox2);
    cp_bbox(cv1, bbox1);
    cp_bbox(cv2, bbox2);

    // The bounding boxes do not overlap - return.
    if (!bbox1.Overlaps(bbox2))
    {
      //iddo debug
      //std::cout << "boxes do not overlap - curve do not intersect"  << std::endl;
      return;
    }
    else // Debug print
    {
      //std::cout << "boxes overlap - continuing..." << std::endl;
    }

    Vector_2 dir1, dir2, leftMost1, rightMost1, leftMost2, rightMost2;

    bool cv1SpanOK = true;
	  bool cv2SpanOK = true;
	  bool spansOverlap = false;

    if (should_check_span)
    {
      cv1SpanOK = _CrvTanAngularSpan(cv1, dir1, leftMost1, rightMost1);
	    cv2SpanOK = _CrvTanAngularSpan(cv2, dir2, leftMost2, rightMost2);
      if (cv1SpanOK && cv2SpanOK)
      {
  	    spansOverlap = _AngularSpansOverlap(leftMost1, rightMost1, leftMost2, rightMost2);
      }
      else
      {
        spansOverlap = true;
      }
    }

    if (!spansOverlap)
    {
      // Debug print.
      //std::cout << "Spans do NOT Overlap." << std::endl;

      // Debug print
      /*
      std::cout << "CV1:\n";
      Control_point_vec::const_iterator debug_iter = cv1.begin();
      while (debug_iter != cv1.end())
      {
        std::cout << "[" << *debug_iter << "]";
        ++debug_iter;
      }
      std::cout << std::endl;
      std::cout << "CV2:\n";
      debug_iter = cv2.begin();
      while (debug_iter != cv2.end())
      {
        std::cout << "[" << *debug_iter << "]";
        ++debug_iter;
      }
      std::cout << std::endl;
      */

      // Checking for endpoint intersections.
      bool endpoint_coincidence = _endpoints_coincide(cv1, l1, r1,
			        cv2, l2, r2,
              endpoint_intersections,
              intersection_pairs);

      if (endpoint_coincidence)
      {
        // Debug print
        //std::cout << "endpoint coincidence intersection." << std::endl; 
        //std::cout << "curves intersect! intersectionIntervals.size() == " << intersection_pairs.size() << std::endl;

        //Since spans do not overlap and there is an endpoint coincidence
        //there is no need for further checking (the endpoint is the intersection).
        return;
      }

      const Point_2& s1 = cv1.front();
      const Point_2& t1 = cv1.back();
      const Point_2& s2 = cv2.front();
      const Point_2& t2 = cv2.back();

      Line_2 skew1a, skew1b;
      _cvSkewBbox(cv1, skew1a, skew1b);

      Line_2 skew2a, skew2b;
      _cvSkewBbox(cv2, skew2a, skew2b);

      // Check for opposite orientations between: 
      // (skew2a, s1)(skew2a, t1)
      // (Skew2b, s1)(Skew2b, t1)
      // and vice versa
      //RI: TODO not use < 0 (count on numeric value) on Oriented_side
      Oriented_side or_2a_s1 = skew2a.oriented_side(s1);
      Oriented_side or_2a_t1 = skew2a.oriented_side(t1);

      Oriented_side or_2b_s1 = skew2b.oriented_side(s1);
      Oriented_side or_2b_t1 = skew2b.oriented_side(t1);

      bool s1_t1_are_opposite = ((or_2a_s1*or_2a_t1 < 0) &&
								      (or_2b_s1*or_2b_t1 < 0));


      Oriented_side or_1a_s2 = skew1a.oriented_side(s2);
      Oriented_side or_1a_t2 = skew1a.oriented_side(t2);

      Oriented_side or_1b_s2 = skew1b.oriented_side(s2);
      Oriented_side or_1b_t2 = skew1b.oriented_side(t2);

      bool s2_t2_are_opposite = ((or_1a_s2*or_1a_t2 < 0) &&
								      (or_1b_s2*or_1b_t2 < 0));

      if (s1_t1_are_opposite && s2_t2_are_opposite) 
      {
        // Construct the bbox from intersection of skew lines.
        Bez_point_bbox inter_bbox;
        Control_point_vec aux_vec;
        Object  res;
        Point_2 p;

        res = intersection (skew1a, skew2a);
        if (!assign(p, res)) 
        {
          CGAL_assertion(false);
        }
        aux_vec.push_back(p);

        res = intersection (skew1a, skew2b);
        if (!assign(p, res)) 
        {
          CGAL_assertion(false);
        }
        aux_vec.push_back(p);

        res = intersection (skew1b, skew2a);
        if (!assign(p, res)) 
        {
          CGAL_assertion(false);
        }
        aux_vec.push_back(p);

        res = intersection (skew1b, skew2b);
        if (!assign(p, res)) 
        {
          CGAL_assertion(false);
        }
        aux_vec.push_back(p);

        cp_bbox(aux_vec, inter_bbox);

        intersection_pairs.push_back(Bound_pair(
                                        Bez_point_bound(cv1, l1, r1, Bez_point_bound::INTERSECTION_PT, true),
                                        Bez_point_bound(cv2, l2, r2, Bez_point_bound::INTERSECTION_PT, true),
                                        inter_bbox
                                      )
                                    );

        // Debug print
        //std::cout << "curves intersect! intersectionIntervals.size() == " << intersection_pairs.size() << std::endl;
        // gets(aux);

	      return;		
      }
    }

    // Subdivide the two curves and recurse. 
    Control_point_vec cv1a, cv1b; 
    
    bisect_control_polygon_2 (cv1.begin(), cv1.end(),
                              std::back_inserter(cv1a),
                              std::front_inserter(cv1b));
    //DeCasteljau(cv1, 0.5, cv1a, cv1b);
    NT t_mid1 = NT(0.5) * (l1 + r1);

    Control_point_vec cv2a, cv2b; 
    
    bisect_control_polygon_2 (cv2.begin(), cv2.end(),
                              std::back_inserter(cv2a),
                              std::front_inserter(cv2b));
    //DeCasteljau(cv2, 0.5, cv2a, cv2b);
    NT t_mid2 = NT(0.5) * (l2 + r2);

    // Recursion:
    _CrvCrvInterAux (cv1a, l1, t_mid1, cv2a, l2, t_mid2, 
                     spansOverlap, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1a, l1, t_mid1, cv2b, t_mid2, r2, 
                     spansOverlap, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1b, t_mid1, r1, cv2a, l2, t_mid2, 
                     spansOverlap, endpoint_intersections,
                     intersection_pairs);
    _CrvCrvInterAux (cv1b, t_mid1, r1, cv2b, t_mid2, r2,
                     spansOverlap, endpoint_intersections,
                     intersection_pairs);
  }

};

//TODO - These classes can both go outside to another file.
//if we decide to take some of the functionality to a template algo file.
template <class _NT> 
struct _Bez_point_bound {
  enum Point_type {RATIONAL_PT, VERTICAL_TANGENCY_PT, INTERSECTION_PT, UNDEFINED};

  typedef _NT                    NT;
  typedef typename Cartesian<NT>::Point_2 Point_2;
  typedef std::deque<Point_2>                   Control_point_vec;

  Control_point_vec bounding_polyline;

  NT t_min; // iddo: in the future, maybe a better rep for these numbers
  NT t_max; // possibly a mantissa+exp representation.
  Point_type point_type;
  bool can_refine;

  _Bez_point_bound() : bounding_polyline(),
    t_min(), t_max(), point_type(UNDEFINED), can_refine(false)
  {}

  _Bez_point_bound(const _Bez_point_bound& o) : bounding_polyline(o.bounding_polyline),
    t_min(o.t_min), t_max(o.t_max), point_type(o.point_type), can_refine(o.can_refine)
  {}

  _Bez_point_bound(const Control_point_vec& bp,
    const NT& min_t, const NT& max_t, 
    Point_type pt, bool cr) : bounding_polyline(bp),
    t_min(min_t), t_max(max_t), point_type(pt), can_refine(cr)
  {}

};

template <class NT> 
struct _Bez_point_bbox {
  NT min_x, max_x, min_y, max_y;

  _Bez_point_bbox() : min_x(), max_x(), min_y(), max_y() {}
  _Bez_point_bbox(const _Bez_point_bbox& other) : 
    min_x(other.min_x), max_x(other.max_x), 
    min_y(other.min_y), max_y(other.max_y) 
  {}

  _Bez_point_bbox(const NT& _min_x, const NT& _max_x,
    const NT& _min_y, const NT& _max_y) : 
  min_x(_min_x), max_x(_max_x), min_y(_min_y), max_y(_max_y) {}

  void AddS(const _Bez_point_bbox& other) {
    min_x = std::min(min_x, other.min_x);
    min_y = std::min(min_y, other.min_y);
    max_x = std::max(max_x, other.max_x);
    max_y = std::max(max_y, other.max_y);
  }

  _Bez_point_bbox Add(const _Bez_point_bbox& other) const {
    _Bez_point_bbox res(other);
    res.AddS(*this);
    return res;
  }

  bool Overlaps(const _Bez_point_bbox& other) const
  {
    if (max_x < other.min_x ||
        max_y < other.min_y ||
        other.max_x < min_x ||
        other.max_y < min_y)
        return false;
    return true;
  }

  bool Overlaps_x(const _Bez_point_bbox& other) const
  {
    if (max_x < other.min_x ||
        other.max_x < min_x)
        return false;
    return true;
  }

};

CGAL_END_NAMESPACE

#endif //CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H
