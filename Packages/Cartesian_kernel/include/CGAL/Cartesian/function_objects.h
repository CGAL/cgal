// Copyright (c) 1999-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra, Sylvain Pion, Michael Hoffmann

#ifndef CGAL_CARTESIAN_FUNCTION_OBJECTS_H
#define CGAL_CARTESIAN_FUNCTION_OBJECTS_H

#include <CGAL/Kernel/function_objects.h>
#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/predicates/kernel_ftC3.h>
#include <CGAL/constructions/kernel_ftC2.h>
#include <CGAL/constructions/kernel_ftC3.h>
#include <CGAL/Cartesian/solve_3.h>

CGAL_BEGIN_NAMESPACE

namespace CartesianKernelFunctors {

  using namespace CommonKernelFunctors;

  template <typename K>
  class Angle_2
  {
    typedef typename K::Point_2 Point_2;
  public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return angleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()); }
  };

  template <typename K>
  class Angle_3
  {
    typedef typename K::Point_3 Point_3;
  public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      return angleC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z(),
		     r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Are_parallel_2
  {
    typedef typename K::Line_2          Line_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename K::Ray_2           Ray_2;
  
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()(const Line_2& l1, const Line_2& l2) const
    { return parallelC2(l1.a(), l1.b(), l2.a(), l2.b()); }

    bool
    operator()(const Segment_2& s1, const Segment_2& s2) const
    { return parallelC2(s1.source().x(), s1.source().y(),
                        s1.target().x(), s1.target().y(),
                        s2.source().x(), s2.source().y(),
                        s2.target().x(), s2.target().y());
    }

    bool
    operator()(const Ray_2& r1, const Ray_2& r2) const
    { return parallelC2(r1.source().x(), r1.source().y(),
                        r1.second_point().x(), r1.second_point().y(),
                        r2.source().x(), r2.source().y(),
                        r2.second_point().x(), r2.second_point().y());
    }
  };

  template <typename K>
  class Are_parallel_3
  {
    typedef typename K::Line_3          Line_3;
    typedef typename K::Segment_3       Segment_3;
    typedef typename K::Ray_3           Ray_3;
    typedef typename K::Plane_3         Plane_3;
  
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()(const Line_3& l1, const Line_3& l2) const
    { return parallelC3(
                l1.to_vector().x(), l1.to_vector().y(), l1.to_vector().z(),
                l2.to_vector().x(), l2.to_vector().y(), l2.to_vector().z());
    }

    bool
    operator()(const Plane_3& h1, const Plane_3& h2) const
    { return parallelC3(h1.a(), h1.b(), h1.c(),
                        h2.a(), h2.b(), h2.c());
    }

    bool
    operator()(const Segment_3& s1, const Segment_3& s2) const
    { return parallelC3(s1.source().x(), s1.source().y(), s1.source().z(),
                        s1.target().x(), s1.target().y(), s1.target().z(),
                        s2.source().x(), s2.source().y(), s2.source().z(),
                        s2.target().x(), s2.target().y(), s2.target().z());
    }

    bool
    operator()(const Ray_3& r1, const Ray_3& r2) const
    { return parallelC3(r1.source().x(), r1.source().y(), r1.source().z(),
	r1.second_point().x(), r1.second_point().y(), r1.second_point().z(),
                        r2.source().x(), r2.source().y(), r2.source().z(),
	r2.second_point().x(), r2.second_point().y(), r2.second_point().z());
    }
  };

  template <typename K>
  class Bounded_side_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Point_3         Point_3;
    typedef typename K::Sphere_3        Sphere_3;
    typedef typename K::Tetrahedron_3   Tetrahedron_3;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
  public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 2 >   Arity;

    Bounded_side
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.bounded_side(p); }

    Bounded_side
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    {
      FT alpha, beta, gamma;

      solve(t.vertex(1)-t.vertex(0),
            t.vertex(2)-t.vertex(0),
            t.vertex(3)-t.vertex(0),
            p - t.vertex(0), alpha, beta, gamma);
      if (   (alpha < 0) || (beta < 0) || (gamma < 0)
          || (alpha + beta + gamma > 1) )
          return ON_UNBOUNDED_SIDE;

      if (   (alpha == 0) || (beta == 0) || (gamma == 0)
          || (alpha+beta+gamma == 1) )
        return ON_BOUNDARY;

      return ON_BOUNDED_SIDE;
    }

    Bounded_side
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.bounded_side(p); }
  };

  template <typename K>
  class Collinear_are_ordered_along_line_2
  {
    typedef typename K::Point_2         Point_2;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Collinear_2 Collinear_2;
    Collinear_2 c;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Collinear_are_ordered_along_line_2() {}
    Collinear_are_ordered_along_line_2(const Collinear_2& c_) : c(c_) {}
#endif // CGAL_kernel_exactness_preconditions 

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      CGAL_kernel_exactness_precondition( c(p, q, r) );
      return collinear_are_ordered_along_lineC2
	(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Collinear_are_ordered_along_line_3
  {
    typedef typename K::Point_3         Point_3;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Collinear_3 Collinear_3;
    Collinear_3 c;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Collinear_are_ordered_along_line_3() {}
    Collinear_are_ordered_along_line_3(const Collinear_3& c_) : c(c_) {}
#endif // CGAL_kernel_exactness_preconditions 

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      CGAL_kernel_exactness_precondition( c(p, q, r) );
      return collinear_are_ordered_along_lineC3(p.x(), p.y(), p.z(),
						q.x(), q.y(), q.z(),
						r.x(), r.y(), r.z());
    }  
  };

  template <typename K>
  class Collinear_are_strictly_ordered_along_line_2
  {
    typedef typename K::Point_2         Point_2;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Collinear_2 Collinear_2;
    Collinear_2 c;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Collinear_are_strictly_ordered_along_line_2() {}
    Collinear_are_strictly_ordered_along_line_2(const Collinear_2& c_) : c(c_) 
    {}
#endif // CGAL_kernel_exactness_preconditions 

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      CGAL_kernel_exactness_precondition( c(p, q, r) );
      return collinear_are_strictly_ordered_along_lineC2
	(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Collinear_are_strictly_ordered_along_line_3
  {
    typedef typename K::Point_3         Point_3;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Collinear_3 Collinear_3;
    Collinear_3 c;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Collinear_are_strictly_ordered_along_line_3() {}
    Collinear_are_strictly_ordered_along_line_3(const Collinear_3& c_) : c(c_)
    {}
#endif // CGAL_kernel_exactness_preconditions 

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      CGAL_kernel_exactness_precondition( c(p, q, r) );
      return collinear_are_strictly_ordered_along_lineC3(p.x(), p.y(), p.z(),
							 q.x(), q.y(), q.z(),
							 r.x(), r.y(), r.z());
    }  
  };

  template <typename K>
  class Collinear_has_on_2
  {
    typedef typename K::Point_2               Point_2;
    typedef typename K::Ray_2                 Ray_2;
    typedef typename K::Segment_2             Segment_2;
    typedef typename K::Construct_point_on_2  Construct_point_on_2;
    typedef typename K::Compare_x_2           Compare_x_2;
    typedef typename K::Compare_y_2           Compare_y_2;
    typedef typename K::Collinear_are_ordered_along_line_2  
    Collinear_are_ordered_along_line_2;
    Construct_point_on_2 cp;
    Compare_x_2 cx;
    Compare_y_2 cy;
    Collinear_are_ordered_along_line_2 co;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    Collinear_has_on_2() {}
    Collinear_has_on_2(const Construct_point_on_2& cp_,
		       const Compare_x_2& cx_,
		       const Compare_y_2& cy_,
		       const Collinear_are_ordered_along_line_2& co_) 
      : cp(cp_), cx(cx_), cy(cy_), co(co_)
    {}

    bool
    operator()( const Ray_2& r, const Point_2& p) const
    {
      Point_2 source = cp(r,0);      
      Point_2 second = cp(r,1);
      switch(cx(source, second)) {
      case SMALLER:
        return cx(source, p) != LARGER;
      case LARGER:
        return cx(p, source) != LARGER;
      default:
        switch(cy(source, second)){
        case SMALLER:
	  return cy(source, p) != LARGER;
        case LARGER:
	  return cy(p, source) != LARGER;
        default:
	  return true; // p == source
        }
      } // switch
    }
  
    bool
    operator()( const Segment_2& s, const Point_2& p) const
    { 
      return co(cp(s,0), p, cp(s,1));
    }
  };

  template <typename K>
  class Collinear_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    Orientation_2 o;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    Collinear_2() {}
    Collinear_2(const Orientation_2 o_) : o(o_) {}

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return o(p, q, r) == COLLINEAR; }
  };

  template <typename K>
  class Collinear_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return collinearC3(p.x(), p.y(), p.z(),
			 q.x(), q.y(), q.z(),
			 r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Compare_angle_with_x_axis_2
  {
    typedef typename K::Direction_2  Direction_2;
  public:
    typedef Comparison_result        result_type;
    typedef Arity_tag< 2 >           Arity;

    Comparison_result
    operator()(const Direction_2& d1, const Direction_2& d2) const
    {
      return compare_angle_with_x_axisC2(d1.dx(), d1.dy(), d2.dx(), d2.dy());
    }
  };

  template <typename K>
  class Compare_distance_2
  {
    typedef typename K::Point_2   Point_2;
  public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 3 >        Arity;

    Comparison_result
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return cmp_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Compare_distance_3
  {
    typedef typename K::Point_3   Point_3;
  public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 3 >        Arity;

    Comparison_result
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      return cmp_dist_to_pointC3(p.x(), p.y(), p.z(),
				 q.x(), q.y(), q.z(),
				 r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Compare_slope_2
  {
    typedef typename K::Line_2     Line_2;
    typedef typename K::Segment_2  Segment_2;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()(const Line_2& l1, const Line_2& l2) const
    { 
      return compare_slopesC2(l1.a(), l1.b(), l2.a(), l2.b());
    }

    Comparison_result
    operator()(const Segment_2& s1, const Segment_2& s2) const
    { 
      return compare_slopesC2(s1.source().x(), s1.source().y(),
			      s1.target().x(), s1.target().y(),
			      s2.source().x(), s2.source().y(),
			      s2.target().x(), s2.target().y());
    }
  };

  template <typename K>
  class Compare_x_at_y_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 3 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_y_at_xC2(p.y(), p.x(), h.b(), h.a(), h.c()); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { 
      return compare_y_at_xC2(p.y(), h1.b(), h1.a(), h1.c(),
			      h2.b(), h2.a(), h2.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { 
      return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
			      h.b(), h.a(), h.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { 
      return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
			      h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
    }
  };

  template <typename K>
  class Compare_xyz_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { 
      return compare_lexicographically_xyzC3(p.x(), p.y(), p.z(),
					     q.x(), q.y(), q.z());
    }
  };

  template <typename K>
  class Compare_xy_2
  {
    typedef typename K::Point_2    Point_2;
  public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y()); }
  };

  template <typename K>
  class Compare_xy_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef Comparison_result  result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y()); }
  };

  template <typename K>
  class Compare_x_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >     Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL_NTS compare(p.x(), q.x()); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l, const Line_2& h) const
    { return compare_xC2(p.x(), l.a(), l.b(), l.c(), h.a(), h.b(), h.c()); }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l.a(), l.b(), l.c(), h1.a(), h1.b(), h1.c(),
			 h2.a(), h2.b(), h2.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { 
      return compare_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
			 h1.a(), h1.b(), h1.c(), h2.a(), h2.b(), h2.c());
    }
  };

  template <typename K>
  class Compare_x_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL_NTS compare(p.x(), q.x()); }
  };

  template <typename K>
  class Compare_y_at_x_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
    typedef typename K::Segment_2  Segment_2;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 3 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_y_at_xC2(p.x(), p.y(), h.a(), h.b(), h.c()); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(p.x(), h1.a(), h1.b(), h1.c(),
			      h2.a(), h2.b(), h2.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    {
      return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
			      h.a(), h.b(), h.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
			      h1.a(), h1.b(), h1.c(), h2.a(), h2.b(), h2.c());
    }

    Comparison_result
    operator()( const Point_2& p, const Segment_2& s) const
    {
      return compare_y_at_xC2(p.x(), p.y(),
			      s.source().x(), s.source().y(),
			      s.target().x(), s.target().y());
    }

    Comparison_result
    operator()( const Point_2& p,
	        const Segment_2& s1, const Segment_2& s2) const
    {
      return compare_y_at_x_segment_C2(p.x(),
				       s1.source().x(), s1.source().y(),
				       s1.target().x(), s1.target().y(),
				       s2.source().x(), s2.source().y(),
				       s2.target().x(), s2.target().y());
    }
  };

  template <typename K>
  class Compare_y_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2    Line_2;
  public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL_NTS compare(p.y(), q.y()); }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { 
      return compare_xC2(p.y(), 
			 l1.b(), l1.a(), l1.c(), 
			 l2.b(), l2.a(), l2.c());
    }

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l.b(), l.a(), l.c(), h1.b(), h1.a(), h1.c(),
			 l.b(), l.a(), l.c(), h2.b(), h2.a(), h2.c());
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
			 h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
    }
  };

  template <typename K>
  class Compare_y_3
  {
    typedef typename K::Point_3   Point_3;
  public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >        Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL_NTS compare(p.y(), q.y()); }
  };

  template <typename K>
  class Compare_z_3
  {
    typedef typename K::Point_3   Point_3;
  public:
    typedef Comparison_result     result_type;
    typedef Arity_tag< 2 >        Arity;

    Comparison_result
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL_NTS compare(p.z(), q.z()); }
  };

  template <typename K>
  class Compute_area_2
  {
    typedef typename K::FT                FT;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Triangle_2        Triangle_2;
    typedef typename K::Point_2           Point_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r ) const
    {
      FT v1x = q.x() - p.x();
      FT v1y = q.y() - p.y();
      FT v2x = r.x() - p.x();
      FT v2y = r.y() - p.y();
      return det2x2_by_formula(v1x, v1y, v2x, v2y)/2;
    }

    FT
    operator()( const Iso_rectangle_2& r ) const
    { return r.area(); }

    FT
    operator()( const Triangle_2& t ) const
    { return t.area(); }
  };

  template <typename K>
  class Compute_scalar_product_2
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_2          Vector_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()(const Vector_2& v, const Vector_2& w) const
    {
	return v.x() * w.x() + v.y() * w.y();
    }
  };

  template <typename K>
  class Compute_scalar_product_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_3          Vector_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()(const Vector_3& v, const Vector_3& w) const
    {
	return v.x() * w.x() + v.y() * w.y() + v.z() * w.z();
    }
  };

  template <typename K>
  class Compute_squared_area_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Point_3           Point_3;
    typedef typename K::Triangle_3        Triangle_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Triangle_3& t ) const
    {
	return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
	return squared_areaC3(p.x(), p.y(), p.z(),
                              q.x(), q.y(), q.z(),
                              r.x(), r.y(), r.z());
    }
  };

  // FIXME
  template <typename K>
  class Compute_squared_distance_Point_Point_2 {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()( const Point_2& p, const Point_2& q) const
    { 
      return squared_distanceC2(p.x(), p.y(), q.x(), q.y());
    }
  };

  // TODO ...
  template <typename K>
  class Compute_squared_radius_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Circle_2& c) const
    { return c.squared_radius(); }

    FT
    operator()( const Point_2& p, const Point_2& q) const
    { return squared_radiusC2(p.x(), p.y(), q.x(), q.y()); }

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return squared_radiusC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()); }
  };

  template <typename K>
  class Compute_squared_radius_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Sphere_3    Sphere_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Sphere_3& s) const
    { return s.squared_radius(); }

    FT
    operator()( const Point_3& p, const Point_3& q) const
    {
      return squared_radiusC3(p.x(), p.y(), p.z(),
			      q.x(), q.y(), q.z());
    }

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return squared_radiusC3(p.x(), p.y(), p.z(),
			      q.x(), q.y(), q.z(),
			      r.x(), r.y(), r.z());
    }

    FT
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      return squared_radiusC3(p.x(), p.y(), p.z(),
			      q.x(), q.y(), q.z(),
			      r.x(), r.y(), r.z(),
			      s.x(), s.y(), s.z());
    }
  };

  template <typename K>
  class Compute_volume_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Iso_cuboid_3   Iso_cuboid_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()(const Point_3& p0, const Point_3& p1,
	       const Point_3& p2, const Point_3& p3) const
    {
      return det3x3_by_formula(p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z(),
                               p2.x()-p0.x(), p2.y()-p0.y(), p2.z()-p0.z(),
                               p3.x()-p0.x(), p3.y()-p0.y(), p3.z()-p0.z())/6;
    }

    FT
    operator()( const Tetrahedron_3& t ) const
    {
      return this->operator()(t.vertex(0), t.vertex(1),
		              t.vertex(2), t.vertex(3));
    }

    FT
    operator()( const Iso_cuboid_3& c ) const
    { return c.volume(); }
  };

  template <typename K>
  class Construct_base_vector_3
  {
    typedef typename K::Vector_3   Vector_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::FT         FT;
    typedef typename K::Construct_cross_product_vector_3
    Construct_cross_product_vector_3;
    typedef typename K::Construct_orthogonal_vector_3 
    Construct_orthogonal_vector_3;
    Construct_cross_product_vector_3 cp;
    Construct_orthogonal_vector_3 co;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Construct_base_vector_3() {}
    Construct_base_vector_3(const Construct_cross_product_vector_3& cp_,
			    const Construct_orthogonal_vector_3& co_)
      : cp(cp_), co(co_)
    {}
  
    Vector_3
    operator()( const Plane_3& h, int index ) const
    {
      if (index == 1) {
	if ( CGAL_NTS is_zero(h.a()) )  // parallel to x-axis
	  return Vector_3(FT(1), FT(0), FT(0));
	 
	if ( CGAL_NTS is_zero(h.b()) )  // parallel to y-axis
	  return Vector_3(FT(0), FT(1), FT(0));
	 
	if ( CGAL_NTS is_zero(h.c()) )  // parallel to z-axis
	  return Vector_3(FT(0), FT(0), FT(1));
	 
	return Vector_3(-h.b(), h.a(), FT(0));
      } else {
	return cp(co(h), this->operator()(h,1));
      }
    }
  };

  template <typename K>
  class Construct_bisector_2
  {
    typedef typename K::FT      FT;
    typedef typename K::Point_2 Point_2;
    typedef typename K::Line_2  Line_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      FT a, b, c;
      bisector_of_pointsC2(p.x(), p.y(), q.x(), q.y(), a, b, c);
      return Line_2(a, b, c);
    }

    Line_2
    operator()(const Line_2& p, const Line_2& q) const
    {
      FT a, b, c;
      bisector_of_linesC2(p.a(), p.b(), p.c(),
                          q.a(), q.b(), q.c(),
                          a, b, c);
      return Line_2(a, b, c);
    }
  };

  template <typename K>
  class Construct_bisector_3
  {
    typedef typename K::FT      FT;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Plane_3
    operator()(const Point_3& p, const Point_3& q) const
    {
      FT a, b, c, d;
      bisector_of_pointsC3(p.x(), p.y(), p.z(),
	                   q.x(), q.y(), q.z(),
			   a, b, c, d);
      return Plane_3(a, b, c, d);
    }

    Plane_3
    operator()(const Plane_3& p, const Plane_3& q) const
    {
      FT a, b, c, d;
      bisector_of_planesC3(p.a(), p.b(), p.c(), p.d(),
                           q.a(), q.b(), q.c(), q.d(),
                           a, b, c, d);
      return Plane_3(a, b, c, d);
    }
  };

  template <typename K>
  class Construct_centroid_2
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
      return construct_point_2(x, y);
    }

    Point_2
    operator()(const Point_2& p, const Point_2& q, 
               const Point_2& r, const Point_2& s) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), s.x(), s.y(), x, y);
      return construct_point_2(x, y);
    }
  };

  template <typename K>
  class Construct_centroid_3
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_3  Point_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      centroidC3(p.x(), p.y(), p.z(),
		 q.x(), q.y(), q.z(),
		 r.x(), r.y(), r.z(),
		 x, y, z);
      return construct_point_3(x, y, z);
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q, 
               const Point_3& r, const Point_3& s) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      centroidC3(p.x(), p.y(), p.z(),
		 q.x(), q.y(), q.z(),
		 r.x(), r.y(), r.z(),
		 s.x(), s.y(), s.z(),
		 x, y, z);
      return construct_point_3(x, y, z);
    }
  };

  template <typename K>
  class Construct_circumcenter_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { 
      typename K::Construct_point_2 construct_point_2;
      typedef typename K::FT        FT;
      FT x, y;
      circumcenterC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
      return construct_point_2(x, y);
    }

    Point_2
    operator()(const Triangle_2& t) const
    { 
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }
  };

  template <typename K>
  class Construct_circumcenter_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Triangle_3     Triangle_3;
    typedef typename K::Point_3        Point_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 4 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      circumcenterC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z(),
		     r.x(), r.y(), r.z(),
		     x, y, z);
      return construct_point_3(x, y, z);
    }

    Point_3
    operator()(const Triangle_3& t) const
    { 
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q,
	       const Point_3& r, const Point_3& s) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      circumcenterC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z(),
		     r.x(), r.y(), r.z(),
		     s.x(), s.y(), s.z(),
		     x, y, z);
      return construct_point_3(x, y, z);
    }

    Point_3
    operator()(const Tetrahedron_3& t) const
    { 
      return this->operator()(t.vertex(0), t.vertex(1),
                              t.vertex(2), t.vertex(3));
    }
  };

  template <typename K>
  class Construct_cross_product_vector_3
  {
    typedef typename K::Vector_3  Vector_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()(const Vector_3& v, const Vector_3& w) const
    {
      return Vector_3(v.y() * w.z() - v.z() * w.y(),
		      v.z() * w.x() - v.x() * w.z(),
		      v.x() * w.y() - v.y() * w.x());
    }
  };

  template <typename K>
  class Construct_lifted_point_3
  {
    typedef typename K::Point_2                    Point_2;
    typedef typename K::Point_3                    Point_3;
    typedef typename K::Plane_3                    Plane_3;
    typedef typename K::Construct_base_vector_3    Construct_base_vector_3;
    typedef typename K::Construct_point_on_3       Construct_point_on_3;
    typedef typename K::Construct_scaled_vector_3  Construct_scaled_vector_3;
    typedef typename K::Construct_translated_point_3  
    Construct_translated_point_3;
    Construct_base_vector_3 cb;
    Construct_point_on_3 cp;
    Construct_scaled_vector_3 cs;
    Construct_translated_point_3 ct;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Construct_lifted_point_3() {}
    Construct_lifted_point_3(const Construct_base_vector_3& cb_,
			     const Construct_point_on_3& cp_,
			     const Construct_scaled_vector_3& cs_,
			     const Construct_translated_point_3& ct_)
      : cb(cb_), cp(cp_), cs(cs_), ct(ct_)
    {}

    Point_3
    operator()(const Plane_3& h, const Point_2& p) const
    {  
      return ct(ct(cp(h), cs(cb(h,1), p.x())), cs(cb(h,2), p.y()));
    }
  };

  template <typename K>
  class Construct_line_2
  {
    typedef typename K::RT                        RT;
    typedef typename K::FT                        FT;
    typedef typename K::Point_2                   Point_2;
    typedef typename K::Direction_2               Direction_2;
    typedef typename K::Vector_2                  Vector_2;
    typedef typename K::Segment_2                 Segment_2;
    typedef typename K::Ray_2                     Ray_2;
    typedef typename K::Line_2                    Line_2;
    typedef typename K::Construct_point_on_2      Construct_point_on_2;
    Construct_point_on_2 c;
  public:
    typedef Line_2            result_type;
    typedef Arity_tag< 2 >    Arity;

    Construct_line_2() {}
    Construct_line_2(const Construct_point_on_2& c_) : c(c_) {}

    Line_2
    operator()() const
    { return Line_2(); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Line_2
    operator()(const RT& a, const RT& b, const RT& cc) const
    { return Line_2(a, b, cc); }
#endif // CGAL_NO_DEPRECATED_CODE

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { 
      FT a, b, cc;
      line_from_pointsC2(p.x(), p.y(), q.x(), q.y(), a, b, cc);
      return Line_2(a, b, cc);
    }

    Line_2
    operator()(const Point_2& p, const Direction_2& d) const
    { 
      FT a, b, cc;
      line_from_point_directionC2(p.x(), p.y(), d.dx(), d.dy(), a, b, cc);
      return Line_2(a, b, cc);
    }

    Line_2
    operator()(const Point_2& p, const Vector_2& v) const
    { 
      FT a, b, cc;
      line_from_point_directionC2(p.x(), p.y(), v.x(), v.y(), a, b, cc);
      return Line_2(a, b, cc);
    }

    Line_2
    operator()(const Segment_2& s) const
    { return this->operator()(c(s, 0), c(s, 1)); }

    Line_2
    operator()(const Ray_2& r) const
    { return this->operator()(c(r, 0), c(r, 1)); }
  };

  template <typename K>
  class Construct_line_3
  {
    typedef typename K::Point_3                   Point_3;
    typedef typename K::Direction_3               Direction_3;
    typedef typename K::Segment_3                 Segment_3;
    typedef typename K::Ray_3                     Ray_3;
    typedef typename K::Line_3                    Line_3;
    typedef typename K::Vector_3                  Vector_3;
    typedef typename K::Construct_vector_3        Construct_vector_3;
    typedef typename K::Construct_direction_3     Construct_direction_3;
    typedef typename K::Construct_point_on_3      Construct_point_on_3;
    Construct_vector_3 cv;
    Construct_point_on_3 cp;
  public:
    typedef Line_3            result_type;
    typedef Arity_tag< 2 >    Arity;

    Construct_line_3() {}
    Construct_line_3(const Construct_vector_3& cv_,
		     const Construct_point_on_3& cp_) 
      : cv(cv_), cp(cp_) 
    {}

    Line_3
    operator()() const
    { return Line_3(); }

    Line_3
    operator()(const Point_3& p, const Point_3& q) const
    { return Line_3(p, cv(p, q)); }

    Line_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return operator()(p, cv(d.dx(), d.dy(), d.dz())); }

    Line_3
    operator()(const Point_3& p, const Vector_3& v) const
    { return Line_3(p, v); }

    Line_3
    operator()(const Segment_3& s) const
    { return Line_3(cp(s,0), cv(cp(s,0), cp(s,1))); }

    Line_3
    operator()(const Ray_3& r) const
    { return Line_3(cp(r,0), cv(cp(r,0), cp(r,1))); }
  };

  template <typename K>
  class Construct_midpoint_2
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_2   Point_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q) const
    { 
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      midpointC2(p.x(), p.y(), q.x(), q.y(), x, y);
      return construct_point_2(x, y);
    }
  };

  template <typename K>
  class Construct_midpoint_3
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_3   Point_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q) const
    { 
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      midpointC3(p.x(), p.y(), p.z(), q.x(), q.y(), q.z(), x, y, z);
      return construct_point_3(x, y, z);
    }
  };

  template <typename K>
  class Construct_opposite_vector_2
  {
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_2
    operator()( const Vector_2& v) const
    { return Vector_2(-v.x(), -v.y()); }
  };

  template <typename K>
  class Construct_opposite_vector_3
  {
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_3
    operator()( const Vector_3& v) const
    { return Vector_3(-v.x(), -v.y(), -v.z()); }
  };

  template <typename K>
  class Construct_orthogonal_vector_3
  {
    typedef typename K::FT FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Vector_3    Vector_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_3
    operator()( const Plane_3& p ) const
    { return Vector_3(p.a(), p.b(), p.c()); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    { 
      FT rpx = p.x()-r.x();
      FT rpy = p.y()-r.y();
      FT rpz = p.z()-r.z();
      FT rqx = q.x()-r.x();
      FT rqy = q.y()-r.y();
      FT rqz = q.z()-r.z();
      // Cross product rp * rq
      FT vx = rpy*rqz - rqy*rpz;
      FT vy = rpz*rqx - rqz*rpx;
      FT vz = rpx*rqy - rqx*rpy;
      typename K::Construct_vector_3 construct_vector;
      
      return construct_vector(vx, vy, vz); 
    }
  };

  template <typename K>
  class Construct_projected_point_3
  {
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::FT         FT;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Line_3& l, const Point_3& p ) const
    {
      // projects p on the line l
      FT lpx = l.point().x();
      FT lpy = l.point().y();
      FT lpz = l.point().z();
      FT ldx = l.direction().dx();
      FT ldy = l.direction().dy();
      FT ldz = l.direction().dz();
      FT dpx = p.x()-lpx;
      FT dpy = p.y()-lpy;
      FT dpz = p.z()-lpz;
      FT lambda = (ldx*dpx+ldy*dpy+ldz*dpz) / (ldx*ldx+ldy*ldy+ldz*ldz);
      return Point_3(lpx + lambda * ldx,
                     lpy + lambda * ldy,
                     lpz + lambda * ldz);
    }

    Point_3
    operator()( const Plane_3& h, const Point_3& p ) const
    { return h.projection(p); }
  };

  template <typename K>
  class Construct_scaled_vector_2
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()( const Vector_2& v, const FT& c) const
    {  
      return Vector_2(c * v.x(), c * v.y());
    }
  };

  template <typename K>
  class Construct_scaled_vector_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_3   Vector_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()( const Vector_3& w, const FT& c) const
    {  
      return Vector_3(c * w.x(), c * w.y(), c * w.z());
    }
  };

  template <typename K>
  class Construct_translated_point_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Vector_2  Vector_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Point_2& p, const Vector_2& v) const
    {  
      typename K::Construct_point_2 construct_point_2;
      return construct_point_2(p.x() + v.x(), p.y() + v.y());
    }
  };

  template <typename K>
  class Construct_translated_point_3
  {
    typedef typename K::Point_3   Point_3;
    typedef typename K::Vector_3  Vector_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Point_3& p, const Vector_3& v) const
    { 
      typename K::Construct_point_3 construct_point_3;
      return construct_point_3(p.x() + v.x(), p.y() + v.y(), p.z() + v.z());
    }
  };

  template <typename K>
  class Construct_vector_2
  {
    typedef typename K::RT           RT;
    typedef typename K::FT           FT;
    typedef typename K::Segment_2    Segment_2;
    typedef typename K::Ray_2        Ray_2;
    typedef typename K::Line_2       Line_2;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Point_2      Point_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()() const
    { return Vector_2(); }

    Vector_2
    operator()( const Point_2& p, const Point_2& q) const
    { return Vector_2(q.x() - p.x(), q.y() - p.y()); }

    Vector_2
    operator()( const Segment_2& s) const
    { return s.to_vector(); }

    Vector_2
    operator()( const Ray_2& r) const
    { return r.to_vector(); }

    Vector_2
    operator()( const Line_2& l) const
    { return l.to_vector(); }

    Vector_2
    operator()( Null_vector) const
    { return Vector_2(FT(0), FT(0)); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Vector_2
    operator()( const RT& x, const RT& y) const
    { return Vector_2(x, y); }

    Vector_2
    operator()( const RT& x, const RT& y, const RT& w) const
    { return Vector_2(x, y, w); }
#endif // CGAL_NO_DEPRECATED_CODE
  };

  template <typename K>
  class Construct_vector_3
  {
    typedef typename K::RT           RT;
    typedef typename K::FT           FT;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Point_3      Point_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()() const
    { return Vector_3(); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q) const
    { 
      return Vector_3(q.x() - p.x(), q.y() - p.y(), q.z() - p.z());
    }

    Vector_3
    operator()( const Segment_3& s) const
    { return s.to_vector(); }

    Vector_3
    operator()( const Ray_3& r) const
    { return r.to_vector(); }

    Vector_3
    operator()( const Line_3& l) const
    { return l.to_vector(); }

    Vector_3
    operator()( const Null_vector&) const
    { return Vector_3(FT(0), FT(0), FT(0)); }

#ifndef CGAL_NO_DEPRECATED_CODE
    Vector_3
    operator()( const RT& x, const RT& y, const RT& z) const
    { return Vector_3(x, y, z); }

    Vector_3
    operator()( const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Vector_3(x, y, z, w); }
#endif // CGAL_NO_DEPRECATED_CODE
  };

  template <typename K>
  class Coplanar_orientation_3
  {
    typedef typename K::Point_3      Point_3;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Coplanar_3   Coplanar_3;
    typedef typename K::Collinear_3  Collinear_3;
    Coplanar_3 cp;
    Collinear_3 cl;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef Orientation  result_type;
    typedef Arity_tag< 4 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Coplanar_orientation_3() {}
    Coplanar_orientation_3(const Coplanar_3& cp_, const Collinear_3& cl_) 
      : cp(cp_), cl(cl_)
    {}
#endif // CGAL_kernel_exactness_preconditions 

    Orientation
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      return coplanar_orientationC3(p.x(), p.y(), p.z(),
				    q.x(), q.y(), q.z(),
				    r.x(), r.y(), r.z());
    }

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      // p,q,r,s supposed to be coplanar
      // p,q,r supposed to be non collinear
      // tests whether s is on the same side of p,q as r
      // returns :
      // COLLINEAR if pqr collinear
      // POSITIVE if qrp and qrs have the same orientation
      // NEGATIVE if qrp and qrs have opposite orientations
      CGAL_kernel_exactness_precondition( ! cl(p, q, r) );
      CGAL_kernel_exactness_precondition( cp(p, q, r, s) );
      return coplanar_orientationC3(p.x(), p.y(), p.z(),
				    q.x(), q.y(), q.z(),
				    r.x(), r.y(), r.z(),
				    s.x(), s.y(), s.z());
    }
  };

  template <typename K>
  class Coplanar_side_of_bounded_circle_3
  {
    typedef typename K::Point_3   Point_3;
#ifdef CGAL_kernel_exactness_preconditions 
    typedef typename K::Coplanar_3   Coplanar_3;
    typedef typename K::Collinear_3  Collinear_3;
    Coplanar_3 cp;
    Collinear_3 cl;
#endif // CGAL_kernel_exactness_preconditions 
  public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 4 >   Arity;

#ifdef CGAL_kernel_exactness_preconditions 
    Coplanar_side_of_bounded_circle_3() {}
    Coplanar_side_of_bounded_circle_3(const Coplanar_3& cp_, 
				      const Collinear_3& cl_) 
      : cp(cp_), cl(cl_)
    {}
#endif // CGAL_kernel_exactness_preconditions 

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    { 
      // p,q,r,t are supposed to be coplanar.
      // p,q,r determine an orientation of this plane (not collinear).
      // returns the equivalent of side_of_bounded_circle(p,q,r,t) 
      // in this plane
      CGAL_kernel_exactness_precondition( cp(p,q,r,t) );
      CGAL_kernel_exactness_precondition( !cl(p,q,r) );
      return coplanar_side_of_bounded_circleC3(p.x(), p.y(), p.z(),
					       q.x(), q.y(), q.z(),
					       r.x(), r.y(), r.z(),
					       t.x(), t.y(), t.z());
    }
  };

  template <typename K>
  class Equal_xy_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { 
      return p.x() == q.x() && p.y() == q.y();
    }
  };

  template <typename K>
  class Equal_x_2
  {
    typedef typename K::Point_2    Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return p.x() == q.x(); }
  };

  template <typename K>
  class Equal_x_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.x() == q.x(); }
  };

  template <typename K>
  class Equal_y_2
  {
    typedef typename K::Point_2    Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return p.y() == q.y(); }
  };

  template <typename K>
  class Equal_y_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.y() == q.y(); }
  };

  template <typename K>
  class Equal_z_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.z() == q.z(); }
  };

  template <typename K>
  class Has_on_3
  {
    typedef typename K::FT               FT;
    typedef typename K::Point_3          Point_3;
    typedef typename K::Vector_3         Vector_3;
    typedef typename K::Line_3           Line_3;
    typedef typename K::Ray_3            Ray_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Triangle_3       Triangle_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Line_3& l, const Point_3& p) const
    { return l.has_on(p); }

    bool
    operator()( const Ray_3& r, const Point_3& p) const
    { return r.has_on(p); }

    bool
    operator()( const Segment_3& s, const Point_3& p) const
    { return s.has_on(p); }

    bool
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.has_on(p); }

    bool
    operator()( const Triangle_3& t, const Point_3& p) const
    {
      Point_3  o  = t.vertex(0) + t.supporting_plane().orthogonal_vector();
      Vector_3 v0 = t.vertex(0)-o,
               v1 = t.vertex(1)-o,
               v2 = t.vertex(2)-o;

      FT alpha, beta, gamma;
      solve(v0, v1, v2, p-o, alpha, beta, gamma);
      return (alpha >= FT(0)) && (beta >= FT(0)) && (gamma >= FT(0))
          && ((alpha+beta+gamma == FT(1)));
    }
  };

  template <typename K>
  class Less_distance_to_point_2
  {
    typedef typename K::Point_2   Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { 
      return has_smaller_dist_to_pointC2(p.x(), p.y(), 
					 q.x(), q.y(), 
					 r.x(), r.y());
    }
  };

  template <typename K>
  class Less_distance_to_point_3
  {
    typedef typename K::Point_3   Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      return has_smaller_dist_to_pointC3(p.x(), p.y(), p.z(),
					 q.x(), q.y(), q.z(),
					 r.x(), r.y(), r.z());
    }
  };

  // TODO ...
  template <typename K>
  class Less_signed_distance_to_line_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2   Line_2;
    typedef typename K::Equal_2 Equal_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 4 >   Arity;

    bool
    operator()(const Point_2& a, const Point_2& b,
               const Point_2& c, const Point_2& d) const
    {
      CGAL_kernel_precondition_code(Equal_2 equal;)
      CGAL_kernel_precondition(! equal(a,b));
      Comparison_result res = cmp_signed_dist_to_lineC2(a.x(), a.y(), 
							b.x(), b.y(),
							c.x(), c.y(),
							d.x(), d.y());

      if ( res == SMALLER ) return true;

      return false;
    }

    bool
    operator()(const Line_2& l, const Point_2& p, const Point_2& q) const
    {
      return has_smaller_signed_dist_to_directionC2(l.a(), l.b(), 
						    p.x(), p.y(),
						    q.x(), q.y());
    }
  };

  template <typename K>
  class Less_signed_distance_to_plane_3
  {
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
    typedef typename K::Collinear_3 Collinear_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()( const Plane_3& h, const Point_3& p, const Point_3& q) const
    { 
      return has_smaller_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
						    p.x(), p.y(), p.z(),
						    q.x(), q.y(), q.z());
    }

    bool
    operator()( const Point_3& hp, const Point_3& hq,  const Point_3& hr,
		const Point_3& p, const Point_3& q) const
    { 
      CGAL_kernel_precondition_code(Collinear_3 collinear_3;)
      CGAL_kernel_precondition(! collinear_3(hp, hq, hr));
      return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
						hq.x(), hq.y(), hq.z(),
						hr.x(), hr.y(), hr.z(),
						p.x(),  p.y(),  p.z(),
						q.x(),  q.y(),  q.z());;
    }
  };

  template <typename K>
  class Less_xyz_3
  {
    typedef typename K::Point_3 Point_3;
    typedef typename K::Compare_xyz_3 Compare_xyz_3;
    Compare_xyz_3 c;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    Less_xyz_3() {}
    Less_xyz_3(const Compare_xyz_3& c_) : c(c_) {}

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_xy_2
  {
    typedef typename K::Point_2 Point_2;
    typedef typename K::Compare_xy_2 Compare_xy_2;
    Compare_xy_2 c;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    Less_xy_2() {}
    Less_xy_2(const Compare_xy_2& c_) : c(c_) {}

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_xy_3
  {
    typedef typename K::Point_3 Point_3;
    typedef typename K::Compare_xy_3 Compare_xy_3;
    Compare_xy_3 c;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    Less_xy_3() {}
    Less_xy_3(const Compare_xy_3& c_) : c(c_) {}

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_x_2
  {
    typedef typename K::Point_2 Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return p.x() < q.x(); }
  };

  template <typename K>
  class Less_x_3
  {
    typedef typename K::Point_3 Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.x() < q.x(); }
  };

  template <typename K>
  class Less_yx_2
  {
    typedef typename K::Point_2       Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { 
      return compare_lexicographically_xyC2(p.y(), p.x(), 
					    q.y(), q.x()) == SMALLER; 
    }
  };

  template <typename K>
  class Less_y_2
  {
    typedef typename K::Point_2 Point_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_2& p, const Point_2& q) const
    { return p.y() < q.y(); }
  };

  template <typename K>
  class Less_y_3
  {
    typedef typename K::Point_3 Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.y() < q.y(); }
  };

  template <typename K>
  class Less_z_3
  {
    typedef typename K::Point_3 Point_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    bool
    operator()( const Point_3& p, const Point_3& q) const
    { return p.z() < q.z(); }
  };

  template <typename K>
  class Orientation_2
  {
    typedef typename K::Point_2 Point_2;
  public:
    typedef Orientation      result_type;
    typedef Arity_tag< 3 >   Arity;

    Orientation
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { 
      return orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Orientation_3
  {
    typedef typename K::Point_3 Point_3;
  public:
    typedef Orientation      result_type;
    typedef Arity_tag< 4 >   Arity;

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      return orientationC3(p.x(), p.y(), p.z(),
			   q.x(), q.y(), q.z(),
			   r.x(), r.y(), r.z(),
			   s.x(), s.y(), s.z());
    }
  };

  template <typename K>
  class Side_of_bounded_circle_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 4 >   Arity;

    Bounded_side
    operator()( const Point_2& p, const Point_2& q, const Point_2& t) const
    { 
      return side_of_bounded_circleC2(p.x(), p.y(), 
				      q.x(), q.y(), 
				      t.x(), t.y());
    }

    Bounded_side
    operator()( const Point_2& p, const Point_2& q,
	        const Point_2& r, const Point_2& t) const
    { 
      return side_of_bounded_circleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(),
				      t.x(), t.y());
    }
  };

  template <typename K>
  class Side_of_bounded_sphere_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef Bounded_side   result_type;
    typedef Arity_tag< 5 >   Arity;

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& test) const
    { 
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
				      q.x(), q.y(), q.z(),
				      test.x(), test.y(), test.z());
    }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& test) const
    {
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
				      q.x(), q.y(), q.z(),
				      r.x(), r.y(), r.z(),
				      test.x(), test.y(), test.z());
    }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& test) const
    {
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
				      q.x(), q.y(), q.z(),
				      r.x(), r.y(), r.z(),
				      s.x(), s.y(), s.z(),
				      test.x(), test.y(), test.z());
    }
  };

  template <typename K>
  class Side_of_oriented_circle_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 4 >   Arity;

    Oriented_side
    operator()( const Point_2& p, const Point_2& q,
	        const Point_2& r, const Point_2& t) const
    {
      return side_of_oriented_circleC2(p.x(), p.y(), 
				       q.x(), q.y(), 
				       r.x(), r.y(),
				       t.x(), t.y());
    }
  };

  template <typename K>
  class Side_of_oriented_sphere_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 5 >   Arity;

    Oriented_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& test) const
    { 
      return side_of_oriented_sphereC3(p.x(), p.y(), p.z(),
				       q.x(), q.y(), q.z(),
				       r.x(), r.y(), r.z(),
				       s.x(), s.y(), s.z(),
				       test.x(), test.y(), test.z());
    }
  };

} // namespace CartesianKernelFunctors

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_FUNCTION_OBJECTS_H
