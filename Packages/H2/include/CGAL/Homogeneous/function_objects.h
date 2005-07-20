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

#ifndef CGAL_HOMOGENEOUS_FUNCTION_OBJECTS_H
#define CGAL_HOMOGENEOUS_FUNCTION_OBJECTS_H

#include <CGAL/Kernel/function_objects.h>
#include <CGAL/Cartesian/function_objects.h>

CGAL_BEGIN_NAMESPACE

namespace HomogeneousKernelFunctors {

  using namespace CommonKernelFunctors;

  // For lazyness...
  using CartesianKernelFunctors::Are_parallel_2;
  using CartesianKernelFunctors::Are_parallel_3;
  using CartesianKernelFunctors::Compute_squared_area_3;
  using CartesianKernelFunctors::Collinear_3;
  using CartesianKernelFunctors::Construct_line_3;

  template <typename K>
  class Angle_2
  {
    typedef typename K::Point_2             Point_2;
    typedef typename K::Construct_vector_2  Construct_vector_2;
    Construct_vector_2 c;
  public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;
  
    Angle_2() {}
    Angle_2(const Construct_vector_2& c_) : c(c_) {}

    Angle
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return (Angle) CGAL_NTS sign(c(q,p) * c(q,r)); } 
    // FIXME: scalar product
  };

  template <typename K>
  class Angle_3
  {
    typedef typename K::Point_3             Point_3;
    typedef typename K::Construct_vector_3  Construct_vector_3;
    Construct_vector_3 c;
  public:
    typedef Angle            result_type;
    typedef Arity_tag< 3 >   Arity;

    Angle_3() {}
    Angle_3(const Construct_vector_3& c_) : c(c_) {}

    Angle
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { return (Angle) CGAL_NTS sign(c(q,p) * c(q,r)); } 
    // FIXME: scalar product
  };


  template <typename K>
  class Bounded_side_2
  {
    typedef typename K::Point_2         Point_2;
    typedef typename K::Circle_2        Circle_2;
    typedef typename K::Triangle_2      Triangle_2;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
  public:
    typedef Bounded_side     result_type;
    typedef Arity_tag< 2 >   Arity;

    Bounded_side
    operator()( const Circle_2& c, const Point_2& p) const
    { 
      typename K::Compute_squared_distance_2 squared_distance;
      return Bounded_side(CGAL_NTS compare(c.squared_radius(),
					   squared_distance(c.center(),p)));
    }

    Bounded_side
    operator()( const Triangle_2& t, const Point_2& p) const
    { 
      typename K::Collinear_are_ordered_along_line_2 
	collinear_are_ordered_along_line;
      typename K::Orientation_2 orientation;
      Orientation o1 = orientation(t.vertex(0), t.vertex(1), p),
	o2 = orientation(t.vertex(1), t.vertex(2), p),
	o3 = orientation(t.vertex(2), t.vertex(3), p);
    
      if (o2 == o1 && o3 == o1)
	return ON_BOUNDED_SIDE;
      return
	(o1 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(0), p, t.vertex(1))) ||
	(o2 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(1), p, t.vertex(2))) ||
	(o3 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(2), p, t.vertex(3)))
	? ON_BOUNDARY
	: ON_UNBOUNDED_SIDE;
    }

    Bounded_side
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    { 
      return  r.rep().bounded_side(p); 
    }

  };

  template <typename K>
  class Bounded_side_3
  {
    typedef typename K::RT              RT;
    typedef typename K::Point_3         Point_3;
    typedef typename K::Vector_3        Vector_3;
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
      RT alpha;
      RT beta ;
      RT gamma;

      Vector_3 v1 = t.vertex(1)-t.vertex(0);
      Vector_3 v2 = t.vertex(2)-t.vertex(0);
      Vector_3 v3 = t.vertex(3)-t.vertex(0);

      Vector_3 vp =   p     - t.vertex(0);

      // want to solve  alpha*v1 + beta*v2 + gamma*v3 == vp
      // let vi' == vi*vi.hw()
      // we solve alpha'*v1' + beta'*v2' + gamma'*v3' == vp' / vp.hw()
      //          muliplied by vp.hw()
      // then we have  alpha = alpha'*v1.hw() / vp.hw()
      // and           beta  = beta' *v2.hw() / vp.hw()
      // and           gamma = gamma'*v3.hw() / vp.hw()

      RT v1x = v1.hx();
      RT v1y = v1.hy();
      RT v1z = v1.hz();
      RT v2x = v2.hx();
      RT v2y = v2.hy();
      RT v2z = v2.hz();
      RT v3x = v3.hx();
      RT v3y = v3.hy();
      RT v3z = v3.hz();
      RT vpx = vp.hx();
      RT vpy = vp.hy();
      RT vpz = vp.hz();

      alpha = det3x3_by_formula( vpx, v2x, v3x,
                                 vpy, v2y, v3y,
                                 vpz, v2z, v3z );
      beta  = det3x3_by_formula( v1x, vpx, v3x,
                                 v1y, vpy, v3y,
                                 v1z, vpz, v3z );
      gamma = det3x3_by_formula( v1x, v2x, vpx,
                                 v1y, v2y, vpy,
                                 v1z, v2z, vpz );
      RT det= det3x3_by_formula( v1x, v2x, v3x,
                                 v1y, v2y, v3y,
                                 v1z, v2z, v3z );

      CGAL_kernel_assertion( det != 0 );
      if (det < 0 )
      {
          alpha = - alpha;
          beta  = - beta ;
          gamma = - gamma;
          det   = - det  ;
      }

      bool t1 = ( alpha < 0 );
      bool t2 = ( beta  < 0 );
      bool t3 = ( gamma < 0 );
            // t1 || t2 || t3 == not contained in cone

      RT lhs = alpha*v1.hw() + beta*v2.hw() + gamma*v3.hw();
      RT rhs = det * vp.hw();

      bool t4 = ( lhs > rhs );
            // alpha + beta + gamma > 1 ?
      bool t5 = ( lhs < rhs );
            // alpha + beta + gamma < 1 ?
      bool t6 = ( (alpha > 0) && (beta > 0) && (gamma > 0) );

      if ( t1 || t2 || t3 || t4 )
      {
          return ON_UNBOUNDED_SIDE;
      }
      return (t5 && t6) ? ON_BOUNDED_SIDE : ON_BOUNDARY;
    }

    Bounded_side
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    { return c.bounded_side(p); }
  };

  template <typename K>
  class Collinear_are_ordered_along_line_2
  {
    typedef typename K::RT              RT;
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

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();

      if ( !(phx * rhw == rhx * phw ) )          // non-vertical ?
	{
	  return !( (  ( phx * qhw < qhx * phw)
		       &&( rhx * qhw < qhx * rhw))
		    ||(  ( qhx * phw < phx * qhw)
			 &&( qhx * rhw < rhx * qhw)) );
	}
      else if ( !(phy * rhw == rhy * phw ) )
	{
	  return !( (  ( phy * qhw < qhy * phw)
		       &&( rhy * qhw < qhy * rhw))
		    ||(  ( qhy * phw < phy * qhw)
			 &&( qhy * rhw < rhy * qhw)) );
	}
      else
	return (( phx*qhw == qhx*phw) && ( phy*qhw == qhy*phw));
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
      typedef typename K::RT RT;
      const RT phx = p.hx();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhw = q.hw();
      const RT rhx = r.hx();
      const RT rhw = r.hw();

      const RT pqx = phx*qhw;
      const RT qpx = qhx*phw;
      const RT prx = phx*rhw;
      const RT qrx = qhx*rhw;
      const RT rqx = rhx*qhw;
      const RT rpx = rhx*phw;

      if ( prx != rpx )   // px != rx
	{
	  //    (px <= qx)&&(qx <= rx) || (px >= qx)&&(qx >= rx)
	  // !(((qx <  px)||(rx <  qx))&&((px <  qx)||(qx <  rx)))
	  return ! (   ((qpx < pqx) || (rqx < qrx))
		       && ((pqx < qpx) || (qrx < rqx))  );
	}

      const RT phy = p.hy();
      const RT qhy = q.hy();
      const RT rhy = r.hy();

      const RT pqy = phy*qhw;
      const RT qpy = qhy*phw;
      const RT pry = phy*rhw;
      const RT qry = qhy*rhw;
      const RT rqy = rhy*qhw;
      const RT rpy = rhy*phw;

      if ( pry != rpy )
	{
	  return ! (   ((qpy < pqy) || (rqy < qry))
		       && ((pqy < qpy) || (qry < rqy))  );
	}

      const RT phz = p.hz();
      const RT qhz = q.hz();
      const RT rhz = r.hz();

      const RT pqz = phz*qhw;
      const RT qpz = qhz*phw;
      const RT prz = phz*rhw;
      const RT qrz = qhz*rhw;
      const RT rqz = rhz*qhw;
      const RT rpz = rhz*phw;

      if ( prz != rpz )
	{
	  return ! (   ((qpz < pqz) || (rqz < qrz))
		       && ((pqz < qpz) || (qrz < rqz))  );
	}
      // p == r
      return  ((rqx == qrx) && (rqy == qry) && (rqz == qrz));
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
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();

      if ( !(phx * rhw == rhx * phw ) )
	{
	  return (   ( phx * qhw < qhx * phw)
		     &&( qhx * rhw < rhx * qhw))
	    ||(   ( qhx * phw < phx * qhw)    // ( phx * qhw > qhx * phw)
		  &&( rhx * qhw < qhx * rhw));  // ( qhx * rhw > rhx * qhw)
	}
      else
	{
	  return (   ( phy * qhw < qhy * phw)
		     &&( qhy * rhw < rhy * qhw))
	    ||(   ( qhy * phw < phy * qhw)    // ( phy * qhw > qhy * phw)
		  &&( rhy * qhw < qhy * rhw));  // ( qhy * rhw > rhy * qhw)
	}
    }
  };

  template <typename K>
  class Collinear_are_strictly_ordered_along_line_3
  {
    typedef typename K::Point_3         Point_3;
    typedef typename K::Direction_3     Direction_3;
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
      if ( p == r) return false;
      Direction_3 dir_pq = (p - q).direction();
      Direction_3 dir_rq = (r - q).direction();
      return (dir_pq == -dir_rq);
    }  // FIXME
  };

  template <typename K>
  class Collinear_has_on_2
  {
    typedef typename K::Point_2               Point_2;
    typedef typename K::Direction_2           Direction_2;
    typedef typename K::Ray_2                 Ray_2;
    typedef typename K::Segment_2             Segment_2;
    typedef typename K::Construct_point_on_2  Construct_point_on_2;
    typedef typename K::Compare_xy_2          Compare_xy_2;
    typedef typename K::Collinear_are_ordered_along_line_2
                                           Collinear_are_ordered_along_line_2;
    Collinear_are_ordered_along_line_2 co;
    Construct_point_on_2 cp;
    Compare_xy_2 cxy;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 2 >   Arity;

    Collinear_has_on_2() {}
    Collinear_has_on_2(const Construct_point_on_2& cp_,
		       const Compare_xy_2& cxy_)
      : cp(cp_), cxy(cxy_)
    {}

    bool
    operator()( const Ray_2& r, const Point_2& p) const
    {
      const Point_2 & source = cp(r,0);      
      return p == source || Direction_2(p - source) == r.direction();
    } // FIXME
  
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
    {
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();
      const RT  RT0 = RT(0);

      // | A B |
      // | C D |

      RT  A = phx*rhw - phw*rhx;
      RT  B = phy*rhw - phw*rhy;
      RT  C = qhx*rhw - qhw*rhx;
      RT  D = qhy*rhw - qhw*rhy;

      RT  det =  A*D - B*C;

      /*
	RT det_old =   p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
	+ p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
	+ p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );
	
	if ( !(CGAL_NTS sign(det) == CGAL_NTS sign(det_old)) )
	{
	std::cerr << "det: " << det << " det_old: " << det_old << flush;
	}
      */
      
      return ( det == RT0 );
    }
  };

  template <typename K>
  class Compare_angle_with_x_axis_2
  {
    typedef typename K::Point_2      Point_2;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Direction_2  Direction_2;
  public:
    typedef Comparison_result        result_type;
    typedef Arity_tag< 2 >           Arity;

    Comparison_result
    operator()(const Direction_2& d1, const Direction_2& d2) const
    {
      typedef typename K::RT  RT;
      CGAL_kernel_precondition(
          static_cast<int>(COUNTERCLOCKWISE) == static_cast<int>(LARGER)
       && static_cast<int>(COLLINEAR)        == static_cast<int>(EQUAL)
       && static_cast<int>(CLOCKWISE)        == static_cast<int>(SMALLER) );

      const RT RT0(0);

      Vector_2 dirvec1(d1.x(), d1.y());      // Added
      Point_2   p1 = CGAL::ORIGIN + dirvec1; // Added
      Vector_2 dirvec2(d2.x(), d2.y());      // Added
      Point_2   p2 = ORIGIN + dirvec2;       // Added
      //  Point_2   p1 = ORIGIN + d1.vector(); // Commented out
      //  Point_2   p2 = ORIGIN + d2.vector(); // Commented out

      CGAL_kernel_precondition( RT0 < p1.hw() );
      CGAL_kernel_precondition( RT0 < p2.hw() );

      int       x_sign1 = static_cast<int>(CGAL_NTS sign( p1.hx() ));
      int       x_sign2 = static_cast<int>(CGAL_NTS sign( p2.hx() ));
      int       y_sign1 = static_cast<int>(CGAL_NTS sign( p1.hy() ));
      int       y_sign2 = static_cast<int>(CGAL_NTS sign( p2.hy() ));

      if ( y_sign1 * y_sign2 < 0)
	{
	  return (0 < y_sign1 ) ? SMALLER : LARGER;
	}

      Point_2   origin( RT0  , RT0   );

      if ( 0 < y_sign1 * y_sign2 )
	{
	  return static_cast<Comparison_result>(static_cast<int>(
			 orientation(origin, p2, p1)));

	  // Precondition on the enums:
	  // COUNTERCLOCKWISE == LARGER   ( ==  1 )
	  // COLLINEAR        == EQUAL    ( ==  0 )
	  // CLOCKWISE        == SMALLER  ( == -1 )
	}
      
      // ( y_sign1 * y_sign2 == 0 )
      
      bool b1 = (y_sign1 == 0) && (x_sign1 >= 0);
      bool b2 = (y_sign2 == 0) && (x_sign2 >= 0);
      
      if ( b1 ) { return  b2 ? EQUAL : SMALLER; }
      if ( b2 ) { return  b1 ? EQUAL : LARGER; }
      if ( y_sign1 == y_sign2 )  // == 0
	{
	  return EQUAL;
	}
      else
	{
	  return (orientation(origin, p1, p2) == COUNTERCLOCKWISE) ?
	    static_cast<Comparison_result>(SMALLER) :
	    static_cast<Comparison_result>(LARGER);
	}
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
      typedef typename K::RT RT;

      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhw = q.hw();
      const RT rhx = r.hx();
      const RT rhy = r.hy();
      const RT rhw = r.hw();
      const RT RT0 = RT(0);
      const RT RT2 = RT(2);

      RT dosd =   // difference of squared distances

	//            phx * phx   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phx * qhx   *   phw * qhw * rhw * rhw
	//   +        qhx * qhx   *   phw * phw * rhw * rhw
	//
	//   +        phy * phy   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phy * qhy   *   phw * qhw * rhw * rhw
	//   +        qhy * qhy   *   phw * phw * rhw * rhw
	//
	// - (        phx * phx   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phx * rhx   *   phw * qhw * qhw * rhw
	//   +        rhx * rhx   *   phw * phw * qhw * qhw
	//
	//   +        phy * phy   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phy * rhy   *   phw * qhw * qhw * rhw
	//   +        rhy * rhy   *   phw * phw * qhw * qhw
	
	rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy )
			    - RT2 * qhw * ( phx*qhx + phy*qhy )
			    )
	- qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy )
			      - RT2 * rhw * ( phx*rhx + phy*rhy )
			      );


      if ( RT0 < dosd )
	{
	  return LARGER;
	}
      else
	{
	  return (dosd < RT0) ? SMALLER : EQUAL;
	}
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
      typedef typename K::RT RT;

      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phz = p.hz();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhz = q.hz();
      const RT qhw = q.hw();
      const RT rhx = r.hx();
      const RT rhy = r.hy();
      const RT rhz = r.hz();
      const RT rhw = r.hw();
      const RT RT0 = RT(0);
      const RT RT2 = RT(2);

      RT dosd =   // difference of squared distances

	rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy + qhz*qhz )
			    - RT2 * qhw * ( phx*qhx + phy*qhy + phz*qhz )
			    )
	- qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy + rhz*rhz )
			      - RT2 * rhw * ( phx*rhx + phy*rhy + phz*rhz )
			      );

      if ( RT0 < dosd )
	{ return LARGER; }
      else
	{ return (dosd < RT0) ? SMALLER : EQUAL; }
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
      typedef typename K::RT RT;
      if (l1.is_horizontal())
	return l2.is_vertical() ? 
	  SMALLER : Comparison_result(CGAL_NTS sign<RT>(l2.a() * l2.b()));
      if (l2.is_horizontal()) 
	return l1.is_vertical() ? 
	  LARGER : Comparison_result(- CGAL_NTS sign<RT>(l1.a() * l1.b()));
      if (l1.is_vertical()) return l2.is_vertical() ? EQUAL : LARGER;
      if (l2.is_vertical()) return SMALLER;
      int l1_sign = CGAL_NTS sign<RT>(-l1.a() * l1.b());
      int l2_sign = CGAL_NTS sign<RT>(-l2.a() * l2.b());

      if (l1_sign < l2_sign) return SMALLER;
      if (l1_sign > l2_sign) return LARGER;

      if (l1_sign > 0)
	return CGAL_NTS compare( CGAL_NTS abs<RT>(l1.a() * l2.b()),
				 CGAL_NTS abs<RT>(l2.a() * l1.b()) );

      return CGAL_NTS compare( CGAL_NTS abs<RT>(l2.a() * l1.b()),
			       CGAL_NTS abs<RT>(l1.a() * l2.b()) );
    } // FIXME

    Comparison_result
    operator()(const Segment_2& s1, const Segment_2& s2) const
    {  
      typedef typename K::FT        FT;

      Comparison_result cmp_y1 = compare_y(s1.source(), s1.target());
      if (cmp_y1 == EQUAL) // horizontal
	{
	  Comparison_result cmp_x2 = compare_x(s2.source(), s2.target());

	  if (cmp_x2 == EQUAL) return SMALLER;
	  FT s_hw = s2.source().hw();
	  FT t_hw = s2.target().hw();
	  return Comparison_result (
	    - CGAL_NTS sign((s2.source().hy()*t_hw - s2.target().hy()*s_hw) *
			    (s2.source().hx()*t_hw - s2.target().hx()*s_hw)) );
	}

      Comparison_result cmp_y2 = compare_y(s2.source(), s2.target());
      if (cmp_y2 == EQUAL)
	{
	  Comparison_result cmp_x1 = compare_x(s1.source(), s1.target());

	  if (cmp_x1 == EQUAL) return LARGER;
	  FT s_hw = s1.source().hw();
	  FT t_hw = s1.target().hw();
	  return Comparison_result (
	    CGAL_NTS sign((s1.source().hy()*t_hw - s1.target().hy()*s_hw) *
			  (s1.source().hx()*t_hw - s1.target().hx()*s_hw)) );
	}

      Comparison_result cmp_x1 = compare_x(s1.source(), s1.target());
      Comparison_result cmp_x2 = compare_x(s2.source(), s2.target());
      if (cmp_x1 == EQUAL)
	return cmp_x2 == EQUAL ? EQUAL : LARGER;

      if (cmp_x2 == EQUAL) return SMALLER;

      FT s1_s_hw = s1.source().hw();
      FT s1_t_hw = s1.target().hw();
      FT s2_s_hw = s2.source().hw();
      FT s2_t_hw = s2.target().hw();
      FT s1_xdiff = s1.source().hx()*s1_t_hw - s1.target().hx()*s1_s_hw;
      FT s1_ydiff = s1.source().hy()*s1_t_hw - s1.target().hy()*s1_s_hw;
      FT s2_xdiff = s2.source().hx()*s2_t_hw - s2.target().hx()*s2_s_hw;
      FT s2_ydiff = s2.source().hy()*s2_t_hw - s2.target().hy()*s2_s_hw;
      Sign s1_sign = CGAL_NTS sign(s1_ydiff * s1_xdiff);
      Sign s2_sign = CGAL_NTS sign(s2_ydiff * s2_xdiff);
 
      if (s1_sign < s2_sign) return SMALLER;
      if (s1_sign > s2_sign) return LARGER;

      if (s1_sign > 0)
	return Comparison_result(
	 CGAL_NTS sign ( CGAL_NTS abs(s1_ydiff * s2_xdiff) -
			 CGAL_NTS abs(s2_ydiff * s1_xdiff)) );

      return Comparison_result(
       CGAL_NTS sign ( CGAL_NTS abs(s2_ydiff * s1_xdiff) -
		       CGAL_NTS abs(s1_ydiff * s2_xdiff)) );
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
    {
      typedef typename K::RT RT;
      CGAL_kernel_precondition( ! h.is_horizontal() );
      Oriented_side ors = h.oriented_side( p );
      if ( h.a() < RT(0) )
	  ors = opposite( ors );
      if ( ors == ON_POSITIVE_SIDE )
	  return LARGER;
      return ( ors == ON_NEGATIVE_SIDE ) ? SMALLER : EQUAL;
    } // FIXME

    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return CGAL_NTS compare(h1.x_at_y( p.y() ), h2.x_at_y( p.y() )); }
    // FIXME

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return compare_x_at_y( gp_linear_intersection( l1, l2 ), h); }
    // FIXME

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_x_at_y( gp_linear_intersection( l1, l2 ), h1, h2 ); }
    // FIXME
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
      typedef typename K::RT RT;
      RT pV = p.hx()*q.hw();
      RT qV = q.hx()*p.hw();
      if ( pV < qV )
	{
	  return SMALLER;
	}
      if ( qV < pV )    //   ( pV > qV )
	{
	  return LARGER;
	}
      // same x
      pV = p.hy()*q.hw();
      qV = q.hy()*p.hw();
      if ( pV < qV )
	{
	  return SMALLER;
	}
      if ( qV < pV )    //   ( pV > qV )
	{
	  return LARGER;
	}
      // same x and y
      pV = p.hz()*q.hw();
      qV = q.hz()*p.hw();
      if ( pV < qV )
	{
	  return SMALLER;
	}
      else
	{
	  return (qV < pV) ? LARGER : EQUAL;
	}
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
    { 
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();

      RT pV = phx*qhw;
      RT qV = qhx*phw;
      if ( pV == qV )
	{
	  pV = phy*qhw;
	  qV = qhy*phw;
	}
      if ( pV < qV )
	{
	  return SMALLER;
	}
      else
	{
	  return (qV < pV) ? LARGER : EQUAL;
	}
    }
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
    { 
      typedef typename K::RT RT;
      RT pV = p.hx()*q.hw();
      RT qV = q.hx()*p.hw();
      if ( pV < qV )
	{
	  return SMALLER;
	}
      if ( qV < pV )    //   ( pV > qV )
	{
	  return LARGER;
	}
      // same x
      pV = p.hy()*q.hw();
      qV = q.hy()*p.hw();
      if ( pV < qV )
	{
	  return SMALLER;
	}
      if ( qV < pV )    //   ( pV > qV )
	{
	  return LARGER;
	}
      // same x and y
      return EQUAL;
    }
  };

  template <typename K>
  class Compare_x_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
  public:
    typedef Comparison_result      result_type;
    typedef Arity_tag< 2 >         Arity;

    Comparison_result
    operator()( const Point_2& p, const Point_2& q) const
    {
      return CGAL_NTS compare(p.hx()*q.hw(), q.hx()*p.hw());
    }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    {
      Point_2 ip = gp_linear_intersection( l1, l2 );
      return this->operator()(p, ip);
    } // FIXME

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return this->operator()(l, h1, l, h2);
    } // FIXME

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { 
      Point_2 lip = gp_linear_intersection( l1, l2 );
      Point_2 hip = gp_linear_intersection( h1, h2 );
      return this->operator()(lip, hip);
    } // FIXME
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
    { return CGAL_NTS compare(p.hx() * q.hw(), q.hx() * p.hw() ); }
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
    { 
      typedef typename K::RT RT;
      CGAL_kernel_precondition( ! h.is_vertical() );
      Oriented_side ors = h.oriented_side( p );
      if ( h.b() < RT(0) )
	{
	  ors = opposite( ors );
	}
      if ( ors == ON_POSITIVE_SIDE )
	{
	  return LARGER;
	}
      return ( ors == ON_NEGATIVE_SIDE ) ? SMALLER : EQUAL;
    } // FIXME
  
    Comparison_result
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    { return CGAL_NTS compare(h1.y_at_x( p.x() ), h2.y_at_x( p.x() )); }
    // FIXME

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    { return compare_y_at_x( gp_linear_intersection( l1, l2 ), h); }
    // FIXME

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    { return compare_y_at_x( gp_linear_intersection( l1, l2 ), h1, h2 ); }
    // FIXME

    Comparison_result
    operator()( const Point_2& p, const Segment_2& s) const
    {
      // compares the y-coordinates of p and the vertical projection of p on s.
      // Precondition : p is in the x-range of s.

      if (compare_x(s.source(), s.target()) == SMALLER) {
        CGAL_kernel_precondition(compare_x(s.source(), p) != LARGER
				 && compare_x(p, s.target()) != LARGER);
        return (Comparison_result) orientation(p, s.source(), s.target());
      }
      else if (compare_x(s.source(), s.target()) == LARGER) {
        CGAL_kernel_precondition(compare_x(s.target(), p) != LARGER
				 && compare_x(p, s.source()) != LARGER);
        return (Comparison_result) orientation(p, s.target(), s.source());
      }
      else {
        CGAL_kernel_precondition(compare_x(s.target(), p) == EQUAL);
	if (compare_y(p, s.source()) == SMALLER &&
	    compare_y(p, s.target()) == SMALLER)
	  return SMALLER;
	if (compare_y(p, s.source()) == LARGER &&
	    compare_y(p, s.target()) == LARGER)
	  return LARGER;
	return EQUAL;
      }
    } // FIXME

    Comparison_result
    operator()( const Point_2& p,
	        const Segment_2& s1, const Segment_2& s2) const
    {
      // compares the y-coordinates of the vertical projections 
      //   of p on s1 and s2
      // Precondition : p is in the x-range of s1 and s2.
      // - if one or two segments are vertical :
      //   - if the segments intersect, return EQUAL
      //   - if not, return the obvious SMALLER/LARGER.

      typedef typename K::FT FT;
      FT px = p.x();
      FT s1sx = s1.source().x();
      FT s1sy = s1.source().y();
      FT s1tx = s1.target().x();
      FT s1ty = s1.target().y();
      FT s2sx = s2.source().x();
      FT s2sy = s2.source().y();
      FT s2tx = s2.target().x();
      FT s2ty = s2.target().y();

      CGAL_kernel_precondition(px >= CGAL_NTS min(s1sx, s1tx) &&
	                       px <= CGAL_NTS max(s1sx, s1tx));
      CGAL_kernel_precondition(px >= CGAL_NTS min(s2sx, s2tx) &&
	                       px <= CGAL_NTS max(s2sx, s2tx));

      if (s1sx != s1tx && s2sx != s2tx) {
	FT s1stx = s1sx-s1tx;
	FT s2stx = s2sx-s2tx;

	return Comparison_result(
	 CGAL_NTS compare(s1sx, s1tx) *
	 CGAL_NTS compare(s2sx, s2tx) *
	 CGAL_NTS compare(-(s1sx-px)*(s1sy-s1ty)*s2stx,
			  (s2sy-s1sy)*s2stx*s1stx
			  -(s2sx-px)*(s2sy-s2ty)*s1stx ));
      }
      else {
	if (s1sx == s1tx) { // s1 is vertical
	  Comparison_result c1, c2;
	  c1 = compare_y_at_x(s1.source(), s2);
	  c2 = compare_y_at_x(s1.target(), s2);
	  if (c1 == c2)
	    return c1;
	  return EQUAL;
	}
	// s2 is vertical
	Comparison_result c3, c4;
	c3 = compare_y_at_x(s2.source(), s1);
	c4 = compare_y_at_x(s2.target(), s1);
	if (c3 == c4)
	  return opposite(c3);
	return EQUAL;
      }
    } // FIXME
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
    { 
      typedef typename K::RT RT;
      
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT  RT0 = RT(0);
      RT com = phy * qhw - qhy * phw;
      if ( com < RT0 )
	  return SMALLER;
      else if ( RT0 < com )
	  return LARGER;
      return EQUAL;
    }

    Comparison_result
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    { 
      Point_2 ip = gp_linear_intersection( l1, l2 );
      return compare_y( p, ip );
    } // FIXME

    Comparison_result
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return this->operator()(l, h1, l, h2);
    }

    Comparison_result
    operator()( const Line_2& l1, const Line_2& l2,
	        const Line_2& h1, const Line_2& h2) const
    {
      Point_2 lip = gp_linear_intersection( l1, l2 );
      Point_2 hip = gp_linear_intersection( h1, h2 );
      return this->operator()( lip, hip );
    } // FIXME
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
    { return CGAL_NTS compare(p.hy() * q.hw(), q.hy() * p.hw() ); }
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
    { return CGAL_NTS compare(p.hz() * q.hw(), q.hz() * p.hw() ); }
  };

  template <typename K>
  class Compute_area_2
  {
    typedef typename K::RT                RT;
    typedef typename K::FT                FT;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Triangle_2        Triangle_2;
    typedef typename K::Point_2           Point_2;
    typedef typename K::Vector_2          Vector_2;
    typedef typename K::Construct_vector_2 Construct_vector_2;
    Construct_vector_2 co;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r ) const
    {
      Vector_2 v1 = co(p, q);
      Vector_2 v2 = co(p, r);

      RT num = v1.hx()*v2.hy() - v2.hx()*v1.hy();
      RT den = RT(2) * v1.hw() * v2.hw();
      return FT(num)/FT(den);
    }

    FT
    operator()( const Iso_rectangle_2& r ) const
    { return (r.xmax()-r.xmin()) * (r.ymax()-r.ymin()); }

    FT
    operator()( const Triangle_2& t ) const
    { return t.area(); }
  };

  template <typename K>
  class Compute_scalar_product_2
  {
    typedef typename K::RT                RT;
    typedef typename K::FT                FT;
    typedef typename K::Vector_2          Vector_2;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()(const Vector_2& v, const Vector_2& w) const
    {
        return FT( RT(v.hx()*w.hx() + v.hy()*w.hy()) ) /
               FT( RT(v.hw()*w.hw() ) );
    }
  };

  template <typename K>
  class Compute_scalar_product_3
  {
    typedef typename K::RT                RT;
    typedef typename K::FT                FT;
    typedef typename K::Vector_3          Vector_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 2 >   Arity;

    FT
    operator()(const Vector_3& v, const Vector_3& w) const
    {
        return FT( RT(v.hx()*w.hx() + v.hy()*w.hy()) + v.hz()*w.hz() ) /
               FT( RT(v.hw()*w.hw() ) );
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
    { return c.rep().squared_radius(); }

    FT
    operator()( const Point_2& p, const Point_2& q) const
    { 
      typedef typename K::FT FT;
      return squared_distance(p, q)/FT(4);
    }  // FIXME

    FT
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return squared_distance(p, circumcenter(p, q, r)); }
    // FIXME
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
      typedef typename K::FT FT;
      return squared_distance(p, q) / FT(4);
    } // FIXME

    FT
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return squared_distance(p, circumcenter(p, q, r));
    } // FIXME

    FT
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      return squared_distance(p, circumcenter(p, q, r, s));
    } // FIXME
  };

  template <typename K>
  class Compute_volume_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Iso_cuboid_3   Iso_cuboid_3;
  public:
    typedef FT               result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()(const Point_3& p0, const Point_3& p1,
	       const Point_3& p2, const Point_3& p3) const
    {
      Vector_3 vec1 = p1 - p0;
      Vector_3 vec2 = p2 - p0;
      Vector_3 vec3 = p3 - p0;

      // first compute (vec1.hw * vec2.hw * vec3.hw * det(vec1, vec2, vec3))
      // then divide by (6 * vec1.hw * vec2.hw * vec3.hw)
      const FT w123 (vec1.hw() * vec2.hw() * vec3.hw());
      const FT& hx1 =  vec1.hx();
      const FT& hy1 =  vec1.hy();
      const FT& hz1 =  vec1.hz();
      const FT& hx2 =  vec2.hx();
      const FT& hy2 =  vec2.hy();
      const FT& hz2 =  vec2.hz();
      const FT& hx3 =  vec3.hx();
      const FT& hy3 =  vec3.hy();
      const FT& hz3 =  vec3.hz();

      return (  (hx1 * (hy2 * hz3 - hy3 * hz2))
              - (hy1 * (hx2 * hz3 - hx3 * hz2))
              + (hz1 * (hx2 * hy3 - hx3 * hy2)))/ (6 * w123);
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
  class Compute_x_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2        Vector_2;


  public:
    typedef  FT result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()(const Point_2& p) const
    {
      return p.rep().x();
    }

    FT
    operator()(const Vector_2& v) const
    {
      return v.rep().x();
    }
  };

  template <typename K>
  class Compute_y_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2        Vector_2;


  public:
    typedef  FT result_type;
    typedef Arity_tag< 1 >   Arity;

    FT
    operator()(const Point_2& p) const
    {
      return p.rep().y();
    }

    FT
    operator()(const Vector_2& v) const
    {
      return v.rep().y();
    }
  };

  template <typename K>
  class Construct_base_vector_3
  {
    typedef typename K::Vector_3   Vector_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::RT         RT;
    typedef typename K::Construct_orthogonal_vector_3 
    Construct_orthogonal_vector_3;
    Construct_orthogonal_vector_3 co;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Construct_base_vector_3() {}
    Construct_base_vector_3(const Construct_orthogonal_vector_3& co_)
      : co(co_)
    {}
  
    Vector_3
    operator()( const Plane_3& h, int index ) const
    {
      if (index == 1) {
	// point():
	// a() != RT0 : Point_3( -d(), RT0, RT0, a() );
	// b() != RT0 : Point_3( RT0, -d(), RT0, b() );
	//            : Point_3( RT0, RT0, -d(), c() );
	// point1():
	// a() != RT0 : Point_3( -b()-d(), a(), RT0, a() );
	// b() != RT0 : Point_3( RT0, -c()-d(), b(), b() );
	//            : Point_3( c(), RT0, -a()-d(), c() );
	 
	const RT RT0(0);
	if ( h.a() != RT0 )
	  {
	    return Vector_3( -h.b(), h.a(), RT0, h.a() );
	  }
	if ( h.b() != RT0 )
	  {
	    return Vector_3( RT0, -h.c(), h.b(), h.b() );
	  }
	CGAL_kernel_assertion ( h.c() != RT(0) );
	return Vector_3( h.c(), RT0, -h.a(), h.c() );
      } else {
	Vector_3 a = co(h);
	Vector_3 b = this->operator()(h, 1);
	return Vector_3(a.hy()*b.hz() - a.hz()*b.hy(),
			a.hz()*b.hx() - a.hx()*b.hz(),
			a.hx()*b.hy() - a.hy()*b.hx(),
			a.hw()*b.hw() );
      }
    }
  };

template <typename K>
  class Construct_bbox_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
    typedef typename K::Circle_2         Circle_2;
  public:
    typedef Bbox_2          result_type;
    typedef Arity_tag< 1 >   Arity;
    
    Bbox_2
    operator()( const Point_2& p) const
    { 
      Interval_nt<> ihx = CGAL_NTS to_interval(p.hx());
      Interval_nt<> ihy = CGAL_NTS to_interval(p.hy());
      Interval_nt<> ihw = CGAL_NTS to_interval(p.hw());
      
      Interval_nt<> ix = ihx/ihw;
      Interval_nt<> iy = ihy/ihw;
      
      return Bbox_2(ix.inf(), iy.inf(), ix.sup(), iy.sup());
    }
    
    Bbox_2
    operator()( const Segment_2& s) const
    { return s.source().bbox() + s.target().bbox(); }

    
    Bbox_2
    operator()( const Triangle_2& t) const
    { 
      typename K::Construct_bbox_2 construct_bbox_2;
      return construct_bbox_2(t.vertex(0)) 
	+ construct_bbox_2(t.vertex(1)) 
	+ construct_bbox_2(t.vertex(2));
    }

    Bbox_2
    operator()( const Iso_rectangle_2& r) const
    { 
      typename K::Construct_bbox_2 construct_bbox_2;
      return construct_bbox_2(r.min()) + construct_bbox_2(r.max()); 
    }
  

    Bbox_2
    operator()( const Circle_2& c) const
    { 
      typename K::Construct_bbox_2 construct_bbox_2;
      Bbox_2 b = construct_bbox_2(c.center());

      Interval_nt<> x (b.xmin(), b.xmax());
      Interval_nt<> y (b.ymin(), b.ymax());
      
      Interval_nt<> sqr = CGAL_NTS to_interval(c.squared_radius());
      Interval_nt<> r = CGAL::sqrt(sqr);
      Interval_nt<> minx = x-r;
      Interval_nt<> maxx = x+r;
      Interval_nt<> miny = y-r;
      Interval_nt<> maxy = y+r;
      
      return Bbox_2(minx.inf(), miny.inf(), maxx.sup(), maxy.sup()); }
  };

  template <typename K>
  class Construct_bisector_2
  {
    typedef typename K::RT      RT;
    typedef typename K::FT      FT;
    typedef typename K::Point_2 Point_2;
    typedef typename K::Line_2  Line_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      // Bisector equation is based on equation
      // ( X - p.x())^2 + (Y - p.y())^2 == ( X - q.x())^2 + (Y - q.y())
      // and x() = hx()/hw() ...

      const RT &phx = p.hx();
      const RT &phy = p.hy();
      const RT &phw = p.hw();
      const RT &qhx = q.hx();
      const RT &qhy = q.hy();
      const RT &qhw = q.hw();

      RT a = RT(2) * ( phx*phw*qhw*qhw - qhx*qhw*phw*phw );
      RT b = RT(2) * ( phy*phw*qhw*qhw - qhy*qhw*phw*phw );
      RT c = qhx*qhx*phw*phw + qhy*qhy*phw*phw 
	   - phx*phx*qhw*qhw - phy*phy*qhw*qhw;

      return Line_2( a, b, c );
    }

    Line_2
    operator()(const Line_2& p, const Line_2& q) const
    {
      RT a, b, c;
      bisector_of_linesC2(p.a(), p.b(), p.c(),
                          q.a(), q.b(), q.c(),
                          a, b, c);
      return Line_2(a, b, c);
    }
  };

  template <typename K>
  class Construct_bisector_3
  {
    typedef typename K::RT      RT;
    typedef typename K::FT      FT;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Plane_3 Plane_3;
  public:
    typedef Plane_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Plane_3
    operator()(const Point_3& p, const Point_3& q) const
    {
      // Bisector equation is based on equation
      // ( X - p.x())^2 + (Y - p.y())^2 == ( X - q.x())^2 + (Y - q.y())
      // and x() = hx()/hw() ...

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phz = p.hz();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhz = q.hz();
      const RT& qhw = q.hw();

      RT a = RT(2) * ( phx*phw*qhw*qhw - qhx*qhw*phw*phw );
      RT b = RT(2) * ( phy*phw*qhw*qhw - qhy*qhw*phw*phw );
      RT c = RT(2) * ( phz*phw*qhw*qhw - qhz*qhw*phw*phw );
      RT d = qhx*qhx*phw*phw + qhy*qhy*phw*phw + qhz*qhz*phw*phw
	   - phx*phx*qhw*qhw - phy*phy*qhw*qhw - phz*phz*qhw*qhw;

      return Plane_3( a, b, c, d );
    }

    Plane_3
    operator()(const Plane_3& p, const Plane_3& q) const
    {
      RT a, b, c, d;
      bisector_of_planesC3(p.a(), p.b(), p.c(), p.d(),
	                   q.a(), q.b(), q.c(), q.d(),
			   a, b, c, d);
      return Plane_3(a, b, c, d);
    }
  };


  template <typename K>
  class Construct_centroid_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      typedef typename K::RT  RT;
      const RT phw(p.hw());
      const RT qhw(q.hw());
      const RT rhw(r.hw());
      RT hx(p.hx()*qhw*rhw + q.hx()*phw*rhw + r.hx()*phw*qhw);
      RT hy(p.hy()*qhw*rhw + q.hy()*phw*rhw + r.hy()*phw*qhw);
      RT hw( phw*qhw*rhw * 3);
      return Point_2(hx, hy, hw);
    }

    Point_2
    operator()(const Triangle_2& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    Point_2
    operator()(const Point_2& p, const Point_2& q, 
               const Point_2& r, const Point_2& s) const
    {
      typedef typename K::RT  RT;
      const RT phw(p.hw());
      const RT qhw(q.hw());
      const RT rhw(r.hw());
      const RT shw(s.hw());
      RT hx(p.hx()*qhw*rhw*shw + q.hx()*phw*rhw*shw + r.hx()*phw*qhw*shw 
	    + s.hx()*phw*qhw*rhw);
      RT hy(p.hy()*qhw*rhw*shw + q.hy()*phw*rhw*shw + r.hy()*phw*qhw*shw
	    + s.hy()*phw*qhw*rhw);
      RT hw( phw*qhw*rhw*shw * 4);
      return Point_2(hx, hy, hw);
    }
  };

  template <typename K>
  class Construct_centroid_3
  {
    typedef typename K::RT             RT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Triangle_3     Triangle_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      const RT& phw = p.hw();
      const RT& qhw = q.hw();
      const RT& rhw = r.hw();
      RT hx(p.hx()*qhw*rhw + q.hx()*phw*rhw + r.hx()*phw*qhw);
      RT hy(p.hy()*qhw*rhw + q.hy()*phw*rhw + r.hy()*phw*qhw);
      RT hz(p.hz()*qhw*rhw + q.hz()*phw*rhw + r.hz()*phw*qhw);
      RT hw( phw*qhw*rhw * RT(3));
      return Point_3(hx, hy, hz, hw);
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q, 
               const Point_3& r, const Point_3& s) const
    {
      const RT& phw = p.hw();
      const RT& qhw = q.hw();
      const RT& rhw = r.hw();
      const RT& shw = s.hw();
      RT hx(p.hx()*qhw*rhw*shw + q.hx()*phw*rhw*shw + r.hx()*phw*qhw*shw
	    + s.hx()*phw*qhw*rhw);
      RT hy(p.hy()*qhw*rhw*shw + q.hy()*phw*rhw*shw + r.hy()*phw*qhw*shw
	    + s.hy()*phw*qhw*rhw);
      RT hz(p.hz()*qhw*rhw*shw + q.hz()*phw*rhw*shw + r.hz()*phw*qhw*shw
	    + s.hz()*phw*qhw*rhw);
      RT hw( phw*qhw*rhw*shw * RT(4));
      return Point_3(hx, hy, hz, hw);
    }

    Point_3
    operator()(const Triangle_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    Point_3
    operator()(const Tetrahedron_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1),
                              t.vertex(2), t.vertex(3));
    }
  };

  template <typename K>
  class Construct_circumcenter_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 3 >   Arity;

    Point_2
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { 
      typedef typename K::RT RT;
      RT phx = p.hx();
      RT phy = p.hy();
      RT phw = p.hw();
      RT qhx = q.hx();
      RT qhy = q.hy();
      RT qhw = q.hw();
      RT rhx = r.hx();
      RT rhy = r.hy();
      RT rhw = r.hw();

#ifdef CGAL_EXPANDED_CIRCUMCENTER_COMPUTATION
      RT vvx =
	( qhy*qhw*phw*phw - phy*phw*qhw*qhw )
	*( phx*phx*rhw*rhw + phy*phy*rhw*rhw - 
	   rhx*rhx*phw*phw - rhy*rhy*phw*phw )
	-  ( rhy*rhw*phw*phw - phy*phw*rhw*rhw )
	*( phx*phx*qhw*qhw + phy*phy*qhw*qhw - 
	   qhx*qhx*phw*phw - qhy*qhy*phw*phw );

      RT vvy =
	-  ( qhx*qhw*phw*phw - phx*phw*qhw*qhw )
	*( phx*phx*rhw*rhw + phy*phy*rhw*rhw - 
	   rhx*rhx*phw*phw - rhy*rhy*phw*phw )
	+  ( rhx*rhw*phw*phw - phx*phw*rhw*rhw )
	*( phx*phx*qhw*qhw + phy*phy*qhw*qhw - 
	   qhx*qhx*phw*phw - qhy*qhy*phw*phw );

      RT vvw = RT(2) *
	(  ( qhx*qhw*phw*phw - phx*phw*qhw*qhw )
	   *( rhy*rhw*phw*phw - phy*phw*rhw*rhw )
	   -  ( rhx*rhw*phw*phw - phx*phw*rhw*rhw )
	   *( qhy*qhw*phw*phw - phy*phw*qhw*qhw ) );
#endif // CGAL_EXPANDED_CIRCUMCENTER_COMPUTATION

      RT qy_py = ( qhy*qhw*phw*phw - phy*phw*qhw*qhw );
      RT qx_px = ( qhx*qhw*phw*phw - phx*phw*qhw*qhw );
      RT rx_px = ( rhx*rhw*phw*phw - phx*phw*rhw*rhw );
      RT ry_py = ( rhy*rhw*phw*phw - phy*phw*rhw*rhw );

      RT px2_py2_rx2_ry_2 =
	phx*phx*rhw*rhw + phy*phy*rhw*rhw - 
	rhx*rhx*phw*phw - rhy*rhy*phw*phw ;
      RT px2_py2_qx2_qy_2 =
	phx*phx*qhw*qhw + phy*phy*qhw*qhw - 
	qhx*qhx*phw*phw - qhy*qhy*phw*phw ;

      RT vvx = qy_py * px2_py2_rx2_ry_2 - ry_py * px2_py2_qx2_qy_2;
      RT vvy = rx_px * px2_py2_qx2_qy_2 - qx_px * px2_py2_rx2_ry_2;
      RT vvw = RT(2) * ( qx_px * ry_py - rx_px * qy_py );

      return Point_2( vvx, vvy, vvw );
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
    typedef typename K::Point_3        Point_3;
    typedef typename K::Triangle_3     Triangle_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Plane_3        Plane_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 4 >   Arity;

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    { 
      return gp_linear_intersection( Plane_3(p,q,r),
				     bisector(p,q),
				     bisector(p,r));
    } // FIXME

    Point_3
    operator()(const Triangle_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q,
	       const Point_3& r, const Point_3& s) const
    {
      typedef typename K::RT RT;

      RT phw( p.hw() );
      RT qhw( q.hw() );
      RT rhw( r.hw() );
      RT shw( s.hw() );

      RT phx( p.hx() );
      RT phy( p.hy() );
      RT phz( p.hz() );
      RT qhx( q.hx() );
      RT qhy( q.hy() );
      RT qhz( q.hz() );
      RT rhx( r.hx() );
      RT rhy( r.hy() );
      RT rhz( r.hz() );
      RT shx( s.hx() );
      RT shy( s.hy() );
      RT shz( s.hz() );

      RT pssq( phx*phx + phy*phy + phz*phz );
      RT qssq( qhx*qhx + qhy*qhy + qhz*qhz );
      RT rssq( rhx*rhx + rhy*rhy + rhz*rhz );
      RT sssq( shx*shx + shy*shy + shz*shz );

      phx *= phw;
      phy *= phw;
      phz *= phw;
      phw *= phw;
      qhx *= qhw;
      qhy *= qhw;
      qhz *= qhw;
      qhw *= qhw;
      rhx *= rhw;
      rhy *= rhw;
      rhz *= rhw;
      rhw *= rhw;
      shx *= shw;
      shy *= shw;
      shz *= shw;
      shw *= shw;

      RT chx =  det4x4_by_formula(phy, phz, pssq, phw,
				  qhy, qhz, qssq, qhw,
				  rhy, rhz, rssq, rhw,
				  shy, shz, sssq, shw );
      RT chy =  det4x4_by_formula(phx, phz, pssq, phw,
				  qhx, qhz, qssq, qhw,
				  rhx, rhz, rssq, rhw,
				  shx, shz, sssq, shw );
      RT chz =  det4x4_by_formula(phx, phy, pssq, phw,
				  qhx, qhy, qssq, qhw,
				  rhx, rhy, rssq, rhw,
				  shx, shy, sssq, shw );
      RT chw =  det4x4_by_formula(phx, phy, phz, phw,
				  qhx, qhy, qhz, qhw,
				  rhx, rhy, rhz, rhw,
				  shx, shy, shz, shw );

      return Point_3( chx, -chy, chz, RT(2)*chw);
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
    operator()(const Vector_3& a, const Vector_3& b) const
    {
      return Vector_3(a.hy()*b.hz() - a.hz()*b.hy(),
		      a.hz()*b.hx() - a.hx()*b.hz(),
		      a.hx()*b.hy() - a.hy()*b.hx(),
		      a.hw()*b.hw() );
    }
  };

  template <typename K>
  class Construct_difference_of_vectors_2
  {
    typedef typename K::Vector_2  Vector_2; 
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()(const Vector_2& v, const Vector_2& w) const
    {
      return Vector_2( v.hx()*w.hw() - w.hx()*v.hw(),
                      v.hy()*w.hw() - w.hy()*v.hw(),
                      v.hw()*w.hw() );
    }
  };

  template <typename K>
  class Construct_direction_2
  {
    typedef typename K::Direction_2     Direction_2;
    typedef typename Direction_2::Rep   Rep;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2        Vector_2;
    typedef typename K::Line_2          Line_2;
    typedef typename K::Ray_2           Ray_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename K::RT              RT;

  public:
    typedef Direction_2       result_type;
    typedef Arity_tag< 1 >    Arity;

    Direction_2
    operator()() const
    { return Rep(); }

    Direction_2
    operator()(const RT& x, const RT& y) const
    { return Rep(x, y); }

    Direction_2
    operator()(const Vector_2& v) const
    { return Rep(v.hx(),v.hy()); }

    Direction_2
    operator()(const Line_2& l) const
    { return Rep(l.b(), -l.a()); }

    Direction_2
    operator()(const Ray_2& r) const
    { return K().construct_direction_2_object()(r.source(), r.second_point()); }

    Direction_2
    operator()(const Segment_2& s) const
    { return K().construct_direction_2_object()(s.source(), s.target()); }

    Direction_2
    operator()(const Point_2& q, const Point_2& p) const
    {
      return Rep( p.hx()*q.hw() - q.hx()*p.hw(),
		  p.hy()*q.hw() - q.hy()*p.hw() );
    }
  };

  template <typename K>
  class Construct_sum_of_vectors_2
  {
    typedef typename K::Vector_2  Vector_2; 
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()(const Vector_2& v, const Vector_2& w) const
    {
      return Vector_2(v.hx()*w.hw() + w.hx()*v.hw(),
                      v.hy()*w.hw() + w.hy()*v.hw(),
                      v.hw()*w.hw());
    }
  };

  template <typename K>
  class Construct_divided_vector_2
  {
    typedef typename K::FT FT;
    typedef typename K::RT RT;
    typedef typename K::Vector_2  Vector_2; 
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()(const Vector_2& v, const FT& f ) const
    {
      return Vector_2( v.hx()*f.denominator(), v.hy()*f.denominator(),
		       v.hw()*f.numerator() ); 
    }

    Vector_2
    operator()(const Vector_2& v, const RT& f ) const
    {
      return Vector_2( v.hx(), v.hy(), v.hw()*f );
    }
  };


  template <typename K>
  class Construct_iso_rectangle_2
  {
    typedef typename K::RT               RT;
    typedef typename K::FT               FT;
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename Iso_rectangle_2::Rep     Rep;

  public:
    typedef Iso_rectangle_2   result_type;
    typedef Arity_tag< 2 >    Arity;

    Iso_rectangle_2
    operator()() const
    { return Rep(); }

    Iso_rectangle_2
    operator()(const Point_2& p, const Point_2& q) const
    {
  bool px_g_qx = ( p.hx()*q.hw() > q.hx()*p.hw() );
  bool py_g_qy = ( p.hy()*q.hw() > q.hy()*p.hw() );

  if ( px_g_qx || py_g_qy)
  {
      if ( px_g_qx && py_g_qy )
      {
          return Rep(q, p);
      }
      else
      {
         if ( px_g_qx )
         {
             return Rep(
             Point_2(q.hx()*p.hw(), p.hy()*q.hw(), q.hw()*p.hw() ),
             Point_2(p.hx()*q.hw(), q.hy()*p.hw(), q.hw()*p.hw() ));
         }
         if ( py_g_qy )
         {
             return Rep(
             Point_2(p.hx()*q.hw(), q.hy()*p.hw(), q.hw()*p.hw() ),
             Point_2(q.hx()*p.hw(), p.hy()*q.hw(), q.hw()*p.hw() ));
         }
      }
  }
  return Rep(p, q);
}

    Iso_rectangle_2
    operator()(const Point_2 &left,   const Point_2 &right,
               const Point_2 &bottom, const Point_2 &top) const
    { 
      CGAL_kernel_assertion_code(typename K::Less_x_2 less_x;)
	CGAL_kernel_assertion_code(typename K::Less_y_2 less_y;)
	CGAL_kernel_assertion(!less_x(right, left));
      CGAL_kernel_assertion(!less_y(top, bottom));
      return Rep(Point_2(left.hx()   * bottom.hw(),
			 bottom.hy() * left.hw(),
			 left.hw()   * bottom.hw()),
		 Point_2(right.hx()  * top.hw(),
			 top.hy()    * right.hw(),
			 right.hw()  * top.hw()));
    }

    Iso_rectangle_2
    operator()(const RT& min_hx, const RT& min_hy, 
	       const RT& max_hx, const RT& max_hy) const
    {
      CGAL_kernel_precondition(min_hx <= max_hx);
      CGAL_kernel_precondition(min_hy <= max_hy);
      return Rep(Point_2(min_hx, min_hy),
		 Point_2(max_hx, max_hy));
    }

    Iso_rectangle_2
    operator()(const RT& min_hx, const RT& min_hy, 
	       const RT& max_hx, const RT& max_hy, const RT& hw) const
    {

	return Rep(Point_2(min_hx, min_hy, hw),
		   Point_2(max_hx, max_hy, hw));
    }


  };


  template <typename K>
  class Construct_lifted_point_3
  {
    typedef typename K::RT                         RT;
    typedef typename K::Point_2                    Point_2;
    typedef typename K::Point_3                    Point_3;
    typedef typename K::Plane_3                    Plane_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()(const Plane_3& h, const Point_2& p) const
    {  
      Point_3 hp( p.hx(), p.hy(), RT(0.0), p.hw());
      return hp.transform( h.transform_to_2d().inverse() );
    }
  };

  template <typename K>
  class Construct_line_2
  {
    typedef typename K::RT                        RT;
    typedef typename K::FT                        FT;
    typedef typename K::Point_2                   Point_2;
    typedef typename K::Vector_2                  Vector_2;
    typedef typename K::Direction_2               Direction_2;
    typedef typename K::Segment_2                 Segment_2;
    typedef typename K::Ray_2                     Ray_2;
    typedef typename K::Line_2                    Line_2;
    typedef typename Line_2::Rep                  Rep;
    typedef typename K::Construct_point_on_2      Construct_point_on_2;
    Construct_point_on_2 cp;
  public:
    typedef Line_2            result_type;
    typedef Arity_tag< 2 >    Arity;

    Construct_line_2() {}
    Construct_line_2(const Construct_point_on_2& cp_) : cp(cp_) {}

    Line_2
    operator()() const
    { return Rep(); }

    Line_2
    operator()(const RT& a, const RT& b, const RT& c) const
    { return Rep(a, b, c); }

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { 
      return Rep(
		    //  a() * X + b() * Y + c() * W() == 0
		    //      |    X        Y       W     |
		    //      |  p.hx()   p.hy()  p.hw()  |
		    //      |  q.hx()   q.hy()  q.hw()  |
		    
		    p.hy()*q.hw() - p.hw()*q.hy(),
		    p.hw()*q.hx() - p.hx()*q.hw(),
		    p.hx()*q.hy() - p.hy()*q.hx() );
    }

    Line_2
    operator()(const Point_2& p, const Vector_2& v) const
    { 
      Point_2 q = p + v;
      return Rep( p.hy()*q.hw() - p.hw()*q.hy(),
		     p.hw()*q.hx() - p.hx()*q.hw(),
		     p.hx()*q.hy() - p.hy()*q.hx() );
    }

    Line_2
    operator()(const Point_2& p, const Direction_2& d) const
    { 
      Point_2 q = p + d.to_vector();
      return Rep( p.hy()*q.hw() - p.hw()*q.hy(),
		     p.hw()*q.hx() - p.hx()*q.hw(),
		     p.hx()*q.hy() - p.hy()*q.hx() );
    }

    Line_2
    operator()(const Segment_2& s) const
    { return this->operator()(cp(s, 0), cp(s, 1)); }

    Line_2
    operator()(const Ray_2& r) const
    { return this->operator()(cp(r, 0), cp(r, 1)); }
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
      typedef typename K::RT RT;
      const RT& phw = p.hw();
      const RT& qhw = q.hw();
      return Point_2( p.hx()*qhw + q.hx()*phw, 
		      p.hy()*qhw + q.hy()*phw,
		      phw * qhw * RT( 2));
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
      typedef typename K::RT RT;
      RT phw = p.hw();
      RT qhw = q.hw();
      return Point_3( p.hx()*qhw + q.hx()*phw,
		      p.hy()*qhw + q.hy()*phw,
		      p.hz()*qhw + q.hz()*phw,
		      RT(2) * phw * qhw );
    }
  };

  // TODO ...
  template <typename K>
  class Construct_opposite_vector_2
  {
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_2
    operator()( const Vector_2& v) const
    { return Vector_2(-v.hx(), -v.hy(), v.hw()); }
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
    { return Vector_3(-v.hx(), -v.hy(), -v.hz(), v.hw()); }
  };

  template <typename K>
  class Construct_orthogonal_vector_3
  {
    typedef typename K::Point_3     Point_3;
    typedef typename K::Vector_3    Vector_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 1 >   Arity;

    Vector_3
    operator()( const Plane_3& p ) const
    { return p.orthogonal_vector(); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
      return operator()(Plane_3(p, q, r));
    }
  };

  template <typename K>
  class Construct_perpendicular_vector_2
  {
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()( const Vector_2& v, Orientation o) const
    { 
      CGAL_kernel_precondition( o != COLLINEAR );
      if (o == COUNTERCLOCKWISE)
	return K().construct_vector_2_object()(-v.hy(), v.hx(), v.hw());
      else
	return K().construct_vector_2_object()(v.hy(), -v.hx(), v.hw());
    }
  };

  template <typename K>
  class Construct_perpendicular_direction_2
  {
    typedef typename K::Direction_2   Direction_2;
  public:
    typedef Direction_2      result_type;
    typedef Arity_tag< 2 >   Arity;

    Direction_2
    operator()( const Direction_2& d, Orientation o) const
    {  
      CGAL_kernel_precondition(o != COLLINEAR);
      if (o == COUNTERCLOCKWISE) {
	  return Direction_2(-d.dy(), d.dx());
	} else {
	  return Direction_2(d.dy(), -d.dx());
	}
    }
  };


  template <typename K>
  class Construct_perpendicular_line_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Point_2   Point_2;
  public:
    typedef Line_2           result_type;
    typedef Arity_tag< 2 >   Arity;

    Line_2
    operator()( const Line_2& l, const Point_2& p) const
    { return typename K::Line_2( -l.b()*p.hw(), l.a()*p.hw(), l.b()*p.hx() - l.a()*p.hy()); }
  };


  template <typename K>
  class Construct_point_2
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_2    Point_2;
    typedef typename K::Vector_2   Vector_2;
    typedef typename K::Line_2     Line_2;
    typedef typename Point_2::Rep  Rep;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 1 >   Arity;

    Point_2
    operator()() const
    { return Rep(); }

    Point_2
    operator()(Origin o) const
    { return Rep(o); }

    Point_2
    operator()(const RT& x, const RT& y) const
    { return Rep(x, y); }

    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { return Rep(x, y, w); }

    Point_2
    operator()(const Line_2& l) const
    { 
      CGAL_kernel_precondition( ! l.is_degenerate() );
      if (l.is_vertical() )
	{
	  return Rep(-l.c(), RT(0)  , l.a() );
	} else {
	  return Rep(RT(0)  , -l.c(), l.b() );
	}
    }
    Point_2
    operator()(const Line_2& l, int i) const
    {
      Point_2 p = K().construct_point_2_object()(l);
      Vector_2 v = K().construct_vector_2_object()(l);
      return K().construct_translated_point_2_object()(p,
						       K().construct_scaled_vector_2_object()(v, RT(i)));    
    }
    
  };

  template <typename K>
  class Construct_projected_point_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Direction_2 Direction_2;
    typedef typename K::Line_2      Line_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_2
    operator()( const Line_2& l, const Point_2& p ) const
    { 
      CGAL_kernel_precondition( ! l.is_degenerate() );
      Line_2  l2( p, Direction_2( l.a(), l.b() ));
      return Point_2( l.b()*l2.c() - l2.b()*l.c(),
		      l2.a()*l.c() - l.a()*l2.c(),
		      l.a()*l2.b() - l2.a()*l.b() );
    }
  };

  template <typename K>
  class Construct_projected_point_3
  {
    typedef typename K::RT         RT;
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::Vector_3   Vector_3;
  public:
    typedef Point_3          result_type;
    typedef Arity_tag< 2 >   Arity;

    Point_3
    operator()( const Line_3& l, const Point_3& p ) const
    {
      if ( l.has_on(p) )
          return p;
      Vector_3  v = p - l.point();
      const RT&  vx = v.hx();
      const RT&  vy = v.hy();
      const RT&  vz = v.hz();
      const RT&  vw = v.hw();
      Vector_3 dir = l.to_vector();
      const RT&  dx = dir.hx();
      const RT&  dy = dir.hy();
      const RT&  dz = dir.hz();
      const RT&  dw = dir.hw();

      RT lambda_num = (vx*dx + vy*dy + vz*dz)*dw; // *dw
      RT lambda_den = (dx*dx + dy*dy + dz*dz)*vw; // *dw

      return l.point() + ( (lambda_num * dir)/lambda_den );
    }

    Point_3
    operator()( const Plane_3& h, const Point_3& p ) const
    { return h.projection(p); }
  };

  template <typename K>
  class Construct_scaled_vector_2
  {
    typedef typename K::RT         RT;
    typedef typename K::FT         FT;
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()( const Vector_2& v, const RT& c) const
    {  
      return Vector_2(c * v.hx(), c * v.hy(), v.hw());
    }

    Vector_2
    operator()( const Vector_2& v, const FT& c) const
    {
      return Vector_2( v.hx()*c.numerator(), v.hy()*c.numerator(),
	               v.hw()*c.denominator() ); }
  };

  template <typename K>
  class Construct_scaled_vector_3
  {
    typedef typename K::RT         RT;
    typedef typename K::FT         FT;
    typedef typename K::Vector_3   Vector_3;
  public:
    typedef Vector_3         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_3
    operator()( const Vector_3& v, const RT& c) const
    {  
      return Vector_3(c * v.hx(), c * v.hy(), c * v.hz(), v.hw());
    }

    Vector_3
    operator()( const Vector_3& v, const FT& c) const
    {
      return Vector_3( v.hx()*c.numerator(), v.hy()*c.numerator(),
                       v.hz()*c.numerator(), v.hw()*c.denominator() ); }
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
      return Point_2( p.hx()*v.hw() + v.hx()*p.hw(),
		      p.hy()*v.hw() + v.hy()*p.hw(),
		      p.hw()*v.hw() );
    }

    Point_2
    operator()( const Origin&, const Vector_2& v) const
    {  
      return Point_2( v.hx(),
		      v.hy(),
		      v.hw() );
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
      return Point_3(p.hx()*v.hw() + v.hx()*p.hw(),
		     p.hy()*v.hw() + v.hy()*p.hw(),
		     p.hz()*v.hw() + v.hz()*p.hw(),
		     p.hw()*v.hw() );
    }

    Point_3
    operator()( const Origin&, const Vector_3& v) const
    {  
      return Point_3( v.hx(),
		      v.hy(),
		      v.hz(),
		      v.hw() );
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
    typedef typename K::Direction_2  Direction_2;
    typedef typename Vector_2::Rep   Rep;
  public:
    typedef Vector_2         result_type;
    typedef Arity_tag< 2 >   Arity;

    Vector_2
    operator()() const
    { return Rep(); }

    Vector_2
    operator()( const Point_2& p, const Point_2& q) const
    { 
      return Rep( q.hx()*p.hw() - p.hx()*q.hw(),
		       q.hy()*p.hw() - p.hy()*q.hw(),
		       p.hw()*q.hw() );
    }

    Vector_2
    operator()( const Origin& , const Point_2& q) const
    { 
      return Rep( q.hx() ,
		       q.hy() ,
		       q.hw() );
    }
    Vector_2
    operator()( const Point_2& p, const Origin& q) const
    { 
      return Rep( - p.hx(),
		       - p.hy(),
		       p.hw() );
    }

    Vector_2
    operator()( const Direction_2& d ) const
    { return Rep(d.dx(), d.dy()); }

    Vector_2
    operator()( const Segment_2& s) const
    { return K().construct_vector_2_object()(s.source(), s.target()); }

    Vector_2
    operator()( const Ray_2& r) const
    { return K().construct_vector_2_object()(r.source(), r.point(1)); }

    Vector_2
    operator()( const Line_2& l) const
    { return K().construct_vector_2_object()( l.b(), -l.a()); }

    Vector_2
    operator()( Null_vector) const
    { return Rep(RT(0), RT(0), RT(1)); }

    Vector_2
    operator()( const RT& x, const RT& y) const
    { return Rep(x, y); }

    Vector_2
    operator()( const RT& x, const RT& y, const RT& w) const
    { return Rep(x, y, w); }
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
      return Vector_3(q.hx()*p.hw() - p.hx()*q.hw(),
		      q.hy()*p.hw() - p.hy()*q.hw(),
		      q.hz()*p.hw() - p.hz()*q.hw(),
		      q.hw()*p.hw() );
    }

    Vector_3
    operator()( const Origin& , const Point_3& q) const
    { 
      return Vector_3( q.hx() ,
		       q.hy() ,
		       q.hz() ,
		       q.hw() );
    }
    Vector_3
    operator()( const Point_3& p, const Origin& q) const
    { 
      return Vector_3( - p.hx(),
		       - p.hy(),
		       - p.hz(),
		       p.hw() );
    }

    Vector_3
    operator()( const Segment_3& s) const
    { return this->operator()(s.start(), s.end()); }

    Vector_3
    operator()( const Ray_3& r) const
    { return r.to_vector(); }

    Vector_3
    operator()( const Line_3& l) const
    { return l.to_vector(); }

    Vector_3
    operator()( const Null_vector&) const
    { return Vector_3(RT(0), RT(0), RT(0), RT(1)); }

// #ifndef CGAL_NO_DEPRECATED_CODE
    Vector_3
    operator()( const RT& x, const RT& y, const RT& z) const
    { return Vector_3(x, y, z); }

    Vector_3
    operator()( const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Vector_3(x, y, z, w); }
// #endif // CGAL_NO_DEPRECATED_CODE
  };

  template <typename K>
  class Construct_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    typedef Point_2          result_type;
    typedef Arity_tag< 2 >   Arity;

    const Point_2 &
    operator()( const Segment_2& s, int i) const
    { return s.vertex(i); }

    const Point_2 &
    operator()( const Triangle_2& t, int i) const
    { return t.rep().vertex(i); }

    const Point_2
    operator()( const Iso_rectangle_2& r, int i) const
    {
      switch (i%4) {
      case 0: return r.min();
      case 1: 
	return Point_2( r.max().hx()*r.min().hw(),
                           r.min().hy()*r.max().hw(),
                           r.min().hw()*r.max().hw() );
      case 2: return r.max();
      default: return Point_2( r.min().hx()*r.max().hw(),
                           r.max().hy()*r.min().hw(),
                           r.min().hw()*r.max().hw() );
      }
    }
      
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
      Orientation oxy_pqr = orientationH2(p.hx(), p.hy(), p.hw(),
					  q.hx(), q.hy(), q.hw(),
					  r.hx(), r.hy(), r.hw());
      if (oxy_pqr != COLLINEAR)
	return oxy_pqr;

      Orientation oyz_pqr = orientationH2(p.hy(), p.hz(), p.hw(),
					  q.hy(), q.hz(), q.hw(),
					  r.hy(), r.hz(), r.hw());
      if (oyz_pqr != COLLINEAR)
	return oyz_pqr;

      return orientationH2(p.hx(), p.hz(), p.hw(),
			   q.hx(), q.hz(), q.hw(),
			   r.hx(), r.hz(), r.hw());
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

      // compute orientation of p,q,s in this plane P:
      Orientation save;
      if ( (save = orientationH2( p.hy(), p.hz(), p.hw(),
				  q.hy(), q.hz(), q.hw(),
				  r.hy(), r.hz(), r.hw())) != COLLINEAR)
	{ return
	    static_cast<Orientation>(
	     static_cast<int>( save)
	     * static_cast<int>( orientationH2( p.hy(), p.hz(), p.hw(),
						q.hy(), q.hz(), q.hw(),
						s.hy(), s.hz(), s.hw())) );
	}
      if ( (save = orientationH2( p.hx(), p.hz(), p.hw(),
				  q.hx(), q.hz(), q.hw(),
				  r.hx(), r.hz(), r.hw())) != COLLINEAR)
	{ return
	    static_cast<Orientation>(
	     static_cast<int>( save)
	     * static_cast<int>( orientationH2( p.hx(), p.hz(), p.hw(),
						q.hx(), q.hz(), q.hw(),
						s.hx(), s.hz(), s.hw())) );
	}
      if ( (save = orientationH2( p.hx(), p.hy(), p.hw(),
				  q.hx(), q.hy(), q.hw(),
				  r.hx(), r.hy(), r.hw())) != COLLINEAR)
	{ return
	    static_cast<Orientation>(
	     static_cast<int>( save)
	     * static_cast<int>( orientationH2( p.hx(), p.hy(), p.hw(),
						q.hx(), q.hy(), q.hw(),
						s.hx(), s.hy(), s.hw())) );
	}
      CGAL_kernel_assertion( false);
      return COLLINEAR;
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
      return (Bounded_side)
	side_of_oriented_sphere(p, q, r, t+cross_product(q-p, r-p), t);
    } // FIXME
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
      return   (p.hx() * q.hw() == q.hx() * p.hw() )
        && (p.hy() * q.hw() == q.hy() * p.hw() );
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
    { return p.hx()*q.hw() == q.hx()*p.hw(); }
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
    { return p.hx()*q.hw() == q.hx()*p.hw(); }
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
    { return p.hy()*q.hw() == q.hy()*p.hw(); }
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
    { return p.hy()*q.hw() == q.hy()*p.hw(); }
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
    { return p.hz()*q.hw() == q.hz()*p.hw(); }
  };

  template <typename K>
  class Has_on_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Line_3           Line_3;
    typedef typename K::Ray_3            Ray_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
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
      if (!t.is_degenerate() )
      {
        Plane_3 sup_pl = t.supporting_plane();
        if ( !sup_pl.has_on(p) )
        {
            return false;
        }
        Tetrahedron_3 tetrapak( t.vertex(0),
                                t.vertex(1),
                                t.vertex(2),
                                t.vertex(0) + sup_pl.orthogonal_vector());
        return tetrapak.has_on_boundary(p);
      }
      Point_3 minp( t.vertex(0) );
      Point_3 maxp( t.vertex(1) );
      if (lexicographically_xyz_smaller(t.vertex(1),t.vertex(0)) )
      {
          minp = t.vertex(1);
          maxp = t.vertex(0);
      }
      if (lexicographically_xyz_smaller(t.vertex(2),minp ) )
      {
          minp = t.vertex(2);
      }
      if (lexicographically_xyz_smaller(maxp, t.vertex(2)) )
      {
          maxp = t.vertex(2);
      }
      if (minp == maxp)
      {
          return (p == maxp);
      }
      Segment_3 s(minp,maxp);
      return s.has_on(p);
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
      typedef typename K::RT RT;

      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhw = q.hw();
      const RT rhx = r.hx();
      const RT rhy = r.hy();
      const RT rhw = r.hw();
      const RT RT0 = RT(0);
      const RT RT2 = RT(2);

      RT dosd =   // difference of squared distances

	//            phx * phx   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phx * qhx   *   phw * qhw * rhw * rhw
	//   +        qhx * qhx   *   phw * phw * rhw * rhw
	//
	//   +        phy * phy   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phy * qhy   *   phw * qhw * rhw * rhw
	//   +        qhy * qhy   *   phw * phw * rhw * rhw
	//
	// - (        phx * phx   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phx * rhx   *   phw * qhw * qhw * rhw
	//   +        rhx * rhx   *   phw * phw * qhw * qhw
	//
	//   +        phy * phy   *   qhw * qhw * rhw * rhw
	//   -RT(2) * phy * rhy   *   phw * qhw * qhw * rhw
	//   +        rhy * rhy   *   phw * phw * qhw * qhw
	
	rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy )
			    - RT2 * qhw * ( phx*qhx + phy*qhy )
			    )
	- qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy )
			      - RT2 * rhw * ( phx*rhx + phy*rhy )
			      );
      
      
      return ( dosd < RT0 );
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
      typedef typename K::RT RT;

      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phz = p.hz();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhz = q.hz();
      const RT qhw = q.hw();
      const RT rhx = r.hx();
      const RT rhy = r.hy();
      const RT rhz = r.hz();
      const RT rhw = r.hw();
      const RT RT0 = RT(0);
      const RT RT2 = RT(2);

      RT dosd =   // difference of squared distances

	rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy + qhz*qhz )
			    - RT2 * qhw * ( phx*qhx + phy*qhy + phz*qhz )
			    )
	- qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy + rhz*rhz )
			      - RT2 * rhw * ( phx*rhx + phy*rhy + phz*rhz )
			      );

      return ( dosd < RT0 );
    }
  };

  template <typename K>
  class Less_signed_distance_to_line_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2   Line_2;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 4 >   Arity;

    bool
    operator()(const Point_2& p, const Point_2& q,
               const Point_2& r, const Point_2& s) const
    {
      typedef typename K::RT RT;

      const RT phx= p.hx();
      const RT phy= p.hy();
      const RT phw= p.hw();
      const RT qhx= q.hx();
      const RT qhy= q.hy();
      const RT qhw= q.hw();
      const RT rhx= r.hx();
      const RT rhy= r.hy();
      const RT rhw= r.hw();
      const RT shx= s.hx();
      const RT shy= s.hy();
      const RT shw= s.hw();
      const RT RT0  = RT(0);

      RT  scaled_dist_r_minus_scaled_dist_s =
	( rhx*shw - shx*rhw ) * (phy*qhw - qhy*phw)
	- ( rhy*shw - shy*rhw ) * (phx*qhw - qhx*phw);

      return ( scaled_dist_r_minus_scaled_dist_s < RT0 );
    }

    bool
    operator()(const Line_2& l, const Point_2& p,
               const Point_2& q) const
    {  
      typedef typename K::RT RT;

      const RT la = l.a();
      const RT lb = l.b();
      const RT phx= p.hx();
      const RT phy= p.hy();
      const RT phw= p.hw();
      const RT qhx= q.hx();
      const RT qhy= q.hy();
      const RT qhw= q.hw();
      const RT RT0 = RT(0);
      
      RT scaled_dist_p_minus_scaled_dist_q =
	la*( phx*qhw - qhx*phw )
	+ lb*( phy*qhw - qhy*phw );
      
      return ( scaled_dist_p_minus_scaled_dist_q < RT0 );

    }
  };

  template <typename K>
  class Less_signed_distance_to_plane_3
  {
    typedef typename K::RT       RT;
    typedef typename K::Point_3  Point_3;
    typedef typename K::Plane_3  Plane_3;
    typedef typename K::Construct_plane_3  Construct_plane_3;
  public:
    typedef bool             result_type;
    typedef Arity_tag< 3 >   Arity;

    bool
    operator()( const Plane_3& pl, const Point_3& p, const Point_3& q) const
    { 
      const RT pla = pl.a();
      const RT plb = pl.b();
      const RT plc = pl.c();
      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phz = p.hz();
      const RT phw = p.hw();
      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhz = q.hz();
      const RT qhw = q.hw();
      const RT RT0 = RT(0);

      RT scaled_dist_p_minus_scaled_dist_q =
	pla*( phx*qhw - qhx*phw )
	+ plb*( phy*qhw - qhy*phw )
	+ plc*( phz*qhw - qhz*phw );

      return ( scaled_dist_p_minus_scaled_dist_q < RT0 );
    }

    bool
    operator()(const Point_3& plp, const Point_3& plq,const Point_3& plr,
	       const Point_3& p, const Point_3& q) const
    { 
      Construct_plane_3 construct_plane_3;
      return operator()(construct_plane_3(plp, plq, plr), p, q);
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
    { return ( p.hx()*q.hw() < q.hx()*p.hw() ); }
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
    { return ( p.hx()*q.hw() < q.hx()*p.hw() ); }
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
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();

      RT pV = phy * qhw;
      RT qV = qhy * phw;
      if ( qV < pV )
	{
	  return false;
	}
      else if ( pV < qV )
	{
	  return true;
	}
      pV = phx * qhw;
      qV = qhx * phw;
      return ( pV < qV );
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
    { return ( p.hy()*q.hw() < q.hy()*p.hw() ); }
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
    { return ( p.hy()*q.hw() < q.hy()*p.hw() ); }
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
    { return   (p.hz() * q.hw() < q.hz() * p.hw() ); }
  };

  template <typename K>
  class Orientation_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Vector_2  Vector_2;
    typedef typename K::Circle_2  Circle_2;
  public:
    typedef Orientation      result_type;
    typedef Arity_tag< 3 >   Arity;

    Orientation
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { 
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();

      // | A B |
      // | C D |

      RT  A = phx*rhw - phw*rhx;
      RT  B = phy*rhw - phw*rhy;
      RT  C = qhx*rhw - qhw*rhx;
      RT  D = qhy*rhw - qhw*rhy;

      return static_cast<Orientation>((int) CGAL_NTS compare(A*D, B*C));
    }

    Orientation
    operator()(const Vector_2& u, const Vector_2& v) const
    { 
      return Orientation (sign_of_determinant2x2(u.hx(), u.hy(),
                                                 v.hx(), v.hy()));
    }

    Orientation
    operator()(const Circle_2& c) const
    { 
      return c.rep().orientation();
    }
  };

  template <typename K>
  class Orientation_3
  {
    typedef typename K::Point_3   Point_3;
    typedef typename K::Vector_3  Vector_3;
  public:
    typedef Orientation      result_type;
    typedef Arity_tag< 4 >   Arity;

    Orientation
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& s) const
    { 
      // Two rows are switched, because of the homogeneous column.
      return (Orientation) 
	sign_of_determinant4x4( p.hx(), p.hy(), p.hz(), p.hw(),
				r.hx(), r.hy(), r.hz(), r.hw(),
				q.hx(), q.hy(), q.hz(), q.hw(),
				s.hx(), s.hy(), s.hz(), s.hw());
    }

    Orientation
    operator()( const Vector_3& u, const Vector_3& v, const Vector_3& w) const
    { 
      return (Orientation) 
	sign_of_determinant3x3( u.hx(), u.hy(), u.hz(),
				v.hx(), v.hy(), v.hz(),
				w.hx(), w.hy(), w.hz());
    }
  };


  template <typename K>
  class Oriented_side_2
  {
    typedef typename K::RT          RT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
    typedef typename K::Line_2      Line_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 2 >   Arity;

    Oriented_side
    operator()( const Circle_2& c, const Point_2& p) const
    { return Oriented_side(c.bounded_side(p) * c.orientation()); }

    Oriented_side
    operator()( const Line_2& l, const Point_2& p) const
    { 
      CGAL_kernel_precondition( ! l.is_degenerate() );
      RT v = l.a()*p.hx() + l.b()*p.hy() + l.c()*p.hw();
      if (v > RT(0)   )
	{
	  return ON_POSITIVE_SIDE;
	}
      else
	{
	  return (v < RT(0)   ) ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
	}
    }

    Oriented_side
    operator()( const Triangle_2& t, const Point_2& p) const
    { 
      typename K::Collinear_are_ordered_along_line_2 
	collinear_are_ordered_along_line;
      typename K::Orientation_2 orientation;
      // depends on the orientation of the vertices
      Orientation o1 = orientation(t.vertex(0), t.vertex(1), p),
	o2 = orientation(t.vertex(1), t.vertex(2), p),
	o3 = orientation(t.vertex(2), t.vertex(3), p),
	ot = orientation(t.vertex(0), t.vertex(1), t.vertex(2));

      if (o1 == ot && o2 == ot && o3 == ot) // ot cannot be COLLINEAR
	return Oriented_side(ot);
      return
	(o1 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(0), p, t.vertex(1))) ||
	(o2 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(1), p, t.vertex(2))) ||
	(o3 == COLLINEAR
	 && collinear_are_ordered_along_line(t.vertex(2), p, t.vertex(3)))
	? ON_ORIENTED_BOUNDARY
	: Oriented_side(-ot);
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
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& thx = t.hx();
      const RT& thy = t.hy();
      const RT& thw = t.hw();

      return Bounded_side( 
	    CGAL_NTS compare((thx*phw-phx*thw)*(qhx*thw-thx*qhw),
			     (thy*phw-phy*thw)*(thy*qhw-qhy*thw)) );
    }

    Bounded_side
    operator()( const Point_2& q, const Point_2& r,
	        const Point_2& s, const Point_2& t) const
    { 
      typedef typename K::RT RT;

      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();
      const RT& shx = s.hx();
      const RT& shy = s.hy();
      const RT& shw = s.hw();
      const RT& thx = t.hx();
      const RT& thy = t.hy();
      const RT& thw = t.hw();
      const RT  RT0 = RT(0);

      CGAL_kernel_precondition( ! collinear(q,r,s) );

      // compute sign of      |qx  qy  qx^2+qy^2  1 |   | a b c d |
      //                      |      --  r  --      | = | e f g h |
      //     determinant      |      --  s  --      | = | i j k l |
      //                      |      --  t  --      |   | m n o p |
      //           where

      RT a = qhx*qhw;
      RT b = qhy*qhw;
      RT c = qhx*qhx + qhy*qhy;
      RT d = qhw*qhw;

      RT e = rhx*rhw;
      RT f = rhy*rhw;
      RT g = rhx*rhx + rhy*rhy;
      RT h = rhw*rhw;

      RT i = shx*shw;
      RT j = shy*shw;
      RT k = shx*shx + shy*shy;
      RT l = shw*shw;

      RT m = thx*thw;
      RT n = thy*thw;
      RT o = thx*thx + thy*thy;
      RT p = thw*thw;

      RT det =   a * ( f*(k*p - l*o) + j*(h*o - g*p) + n*(g*l - h*k) )
	- e * ( b*(k*p - l*o) + j*(d*o - c*p) + n*(c*l - d*k) )
	+ i * ( b*(g*p - h*o) + f*(d*o - c*p) + n*(c*h - d*g) )
	- m * ( b*(g*l - h*k) + f*(d*k - c*l) + j*(c*h - d*g) );

      if ( det == RT0 )
	  return ON_BOUNDARY;
      else
	{
	  if (orientation(q,r,s) == CLOCKWISE)
	      det = -det;
	  return (RT0 < det ) ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
	}
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
    operator()( const Point_3& p, const Point_3& q, const Point_3& t) const
    { 
      typedef typename K::RT RT;

      const RT& phx = p.hx();
      const RT& phy = p.hy();
      const RT& phz = p.hz();
      const RT& phw = p.hw();
      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhz = q.hz();
      const RT& qhw = q.hw();
      const RT& thx = t.hx();
      const RT& thy = t.hy();
      const RT& thz = t.hz();
      const RT& thw = t.hw();

      return Bounded_side( 
	  CGAL_NTS sign<RT>((thx*phw-phx*thw)*(qhx*thw-thx*qhw)
			    + (thy*phw-phy*thw)*(qhy*thw-thy*qhw)
			    + (thz*phw-phz*thw)*(qhz*thw-thz*qhw)));
    }

    Bounded_side
    operator()( const Point_3& p, const Point_3& q,
	        const Point_3& r, const Point_3& t) const
    {
      Point_3 center = circumcenter(p, q, r);
      return Bounded_side( compare_distance_to_point(center, p, t) );
    } // FIXME

    Bounded_side
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
	        const Point_3& s, const Point_3& test) const
    {
      Oriented_side  oside = side_of_oriented_sphere(p,q,r,s,test);
      if ( are_positive_oriented( p,q,r,s) )
	{
	  switch (oside)
	    {
	    case ON_POSITIVE_SIDE    :   return ON_BOUNDED_SIDE;
	    case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
	    case ON_NEGATIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
	    }
	}
      else
	{
	  switch (oside)
	    {
	    case ON_POSITIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
	    case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
	    case ON_NEGATIVE_SIDE    :   return ON_BOUNDED_SIDE;
	    }
	}
      return ON_BOUNDARY;  // Pls, no warnings anylonger
    } // FIXME
  };

  template <typename K>
  class Side_of_oriented_circle_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef Oriented_side    result_type;
    typedef Arity_tag< 4 >   Arity;

    Oriented_side
    operator()( const Point_2& q, const Point_2& r,
	        const Point_2& s, const Point_2& t) const
    {
      typedef typename K::RT RT;

      const RT& qhx = q.hx();
      const RT& qhy = q.hy();
      const RT& qhw = q.hw();
      const RT& rhx = r.hx();
      const RT& rhy = r.hy();
      const RT& rhw = r.hw();
      const RT& shx = s.hx();
      const RT& shy = s.hy();
      const RT& shw = s.hw();
      const RT& thx = t.hx();
      const RT& thy = t.hy();
      const RT& thw = t.hw();

      CGAL_kernel_precondition( ! collinear(q,r,s) );

      // compute sign of      |qx  qy  qx^2+qy^2  1 |   | a b c d |
      //                      |      --  r  --      | = | e f g h |
      //     determinant      |      --  s  --      | = | i j k l |
      //                      |      --  t  --      |   | m n o p |
      //           where

      RT a = qhx*qhw;
      RT b = qhy*qhw;
      RT c = qhx*qhx + qhy*qhy;
      RT d = qhw*qhw;

      RT e = rhx*rhw;
      RT f = rhy*rhw;
      RT g = rhx*rhx + rhy*rhy;
      RT h = rhw*rhw;

      RT i = shx*shw;
      RT j = shy*shw;
      RT k = shx*shx + shy*shy;
      RT l = shw*shw;

      RT m = thx*thw;
      RT n = thy*thw;
      RT o = thx*thx + thy*thy;
      RT p = thw*thw;

      RT det =   a * ( f*(k*p - l*o) + j*(h*o - g*p) + n*(g*l - h*k) )
	- e * ( b*(k*p - l*o) + j*(d*o - c*p) + n*(c*l - d*k) )
	+ i * ( b*(g*p - h*o) + f*(d*o - c*p) + n*(c*h - d*g) )
	- m * ( b*(g*l - h*k) + f*(d*k - c*l) + j*(c*h - d*g) );

      return static_cast<Oriented_side>((int) CGAL_NTS sign(det));
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
	        const Point_3& s, const Point_3& t) const
    { 
      typedef typename K::RT RT;
      CGAL_kernel_precondition( !coplanar(p,q,r,s) );
      const RT phx = p.hx();
      const RT phy = p.hy();
      const RT phz = p.hz();
      const RT phw = p.hw();
      const RT phw2 = phw*phw;

      const RT qhx = q.hx();
      const RT qhy = q.hy();
      const RT qhz = q.hz();
      const RT qhw = q.hw();
      const RT qhw2 = qhw*qhw;

      const RT rhx = r.hx();
      const RT rhy = r.hy();
      const RT rhz = r.hz();
      const RT rhw = r.hw();
      const RT rhw2 = rhw*rhw;

      const RT shx = s.hx();
      const RT shy = s.hy();
      const RT shz = s.hz();
      const RT shw = s.hw();
      const RT shw2 = shw*shw;

      const RT thx = t.hx();
      const RT thy = t.hy();
      const RT thz = t.hz();
      const RT thw = t.hw();
      const RT thw2 = thw*thw;

      const RT det = det5x5_by_formula<RT>(
	   phx*phw, phy*phw, phz*phw, phx*phx + phy*phy + phz*phz, phw2,
	   qhx*qhw, qhy*qhw, qhz*qhw, qhx*qhx + qhy*qhy + qhz*qhz, qhw2,
	   rhx*rhw, rhy*rhw, rhz*rhw, rhx*rhx + rhy*rhy + rhz*rhz, rhw2,
	   shx*shw, shy*shw, shz*shw, shx*shx + shy*shy + shz*shz, shw2,
	   thx*thw, thy*thw, thz*thw, thx*thx + thy*thy + thz*thz, thw2);
      if (det < RT(0))
	{
	  return ON_POSITIVE_SIDE;
	}
      else
	{
	  return (RT(0) < det) ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
	}
    }
  };

} // namespace HomogeneousKernelFunctors

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_FUNCTION_OBJECTS_H
