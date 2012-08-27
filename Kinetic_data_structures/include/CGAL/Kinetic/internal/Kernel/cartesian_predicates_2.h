// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_CARTESIAN_PREDICATES_2_H_
#define CGAL_KINETIC_CARTESIAN_PREDICATES_2_H_
#include <CGAL/Kinetic/basic.h>
#include <CGAL/determinant.h>

namespace CGAL { namespace Kinetic { namespace internal {

template <class KK>
struct Cartesian_orientation_2
{
  Cartesian_orientation_2(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_2 first_argument_type;
  typedef typename KK::Point_2 second_argument_type;
  typedef typename KK::Point_2 third_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c) const
  {
    typedef typename KK::Certificate_function FT;
    FT a00= a.x();
    FT a01= a.y();
    FT a10= b.x();
    FT a11= b.y();
    FT a20= c.x();
    FT a21= c.y();

    // First compute the det2x2
    const FT m01 = a00*a11 - a10*a01;
    const FT m02 = a00*a21 - a20*a01;
    const FT m12 = a10*a21 - a20*a11;
    // Now compute the minors of rank 3
    const FT m012 =  m01 - m02 + m12;
    //std::cout << "Orientation 2 is " << m012 << std::endl;
    return m012;
  }
};

template <class KK>
struct Cartesian_side_of_oriented_circle_2
{
  Cartesian_side_of_oriented_circle_2(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_2 first_argument_type;
  typedef typename KK::Point_2 second_argument_type;
  typedef typename KK::Point_2 third_argument_type;
  typedef typename KK::Point_2 fourth_argument_type;
  result_type operator()(const first_argument_type &ap,
			 const second_argument_type &bp,
			 const third_argument_type &cp,
			 const fourth_argument_type &dp) const
  {
    typedef typename KK::Certificate_function RT;

#if 0
    const RT qhx = ap.hx();
    const RT qhy = ap.hy();
    const RT rhx = bp.hx();
    const RT rhy = bp.hy();
    const RT shx = cp.hx();
    const RT shy = cp.hy();
    const RT thx = dp.hx();
    const RT thy = dp.hy();

    // compute sign of      |qx  qy  qx^2+qy^2  1 |   | a b c d |
    //                      |      --  r  --      | = | e f g h |
    //     determinant      |      --  s  --      | = | i j k l |
    //                      |      --  t  --      |   | m n o p |
    //           where

    RT a = qhx;
    RT b = qhy;
    RT c = qhx*qhx + qhy*qhy;

    RT e = rhx;
    RT f = rhy;
    RT g = rhx*rhx + rhy*rhy;

    RT i = shx;
    RT j = shy;
    RT k = shx*shx + shy*shy;

    RT m = thx;
    RT n = thy;
    RT o = thx*thx + thy*thy;

    RT det =   a * ( f*(k - o) + j*(o - g) + n*(g - k) )
      - e * ( b*(k - o) + j*(o - c) + n*(c - k) )
      + i * ( b*(g - o) + f*(o - c) + n*(c - g) )
      - m * ( b*(g - k) + f*(k - c) + j*(c - g) );
#endif
    RT qpx = bp.x()-ap.x();
    RT qpy = bp.y()-ap.y();
    RT rpx = cp.x()-ap.x();
    RT rpy = cp.y()-ap.y();
    RT tpx = dp.x()-ap.x();
    RT tpy = dp.y()-ap.y();
	  
    RT det= (qpx*tpy - qpy*tpx)*(rpx*(cp.x()-bp.x()) + rpy*(cp.y()-bp.y()))
      -(tpx*(dp.x()-bp.x()) + tpy*(dp.y()-bp.y()))*( qpx*rpy - qpy*rpx) ;
    //std::cout << "New is " << nret << std::endl;
    //std::cout << "Old is " << det << std::endl;

    return det;
  }
};


template <class KK>
struct Cartesian_compare_distance_2
{
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_2 first_argument_type;
  typedef typename KK::Point_2 second_argument_type;
  typedef typename KK::Point_2 third_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const second_argument_type &c) const
  {
    typename KK::Motion_function db= CGAL::square(b.x()-a.x()) + CGAL::square(b.y()-a.y());
    typename KK::Motion_function dc= CGAL::square(c.x()-a.x()) + CGAL::square(c.y()-a.y());
    return dc-db;
  }
};

template <class KK>
struct Cartesian_compare_distance_3
{
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  typedef typename KK::Point_3 third_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const second_argument_type &c) const
  {
    typename KK::Motion_function db= CGAL::square(b.x()-a.x()) + CGAL::square(b.y()-a.y()) 
      + CGAL::square(b.z()-a.z());
    typename KK::Motion_function dc= CGAL::square(c.x()-a.x()) + CGAL::square(c.y()-a.y())
      + CGAL::square(c.z()-a.z());
    return dc-db;
  }
};

template <class KK>
struct Cartesian_less_x_1
{
  Cartesian_less_x_1(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_1 first_argument_type;
  typedef typename KK::Point_1 second_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b) const
  {
    //std::cout << "Evaluating compare of " << a << " and " << b << std::endl;
    return a.x() - b.x();
  }
};

template <class KK>
struct Cartesian_less_x_2
{
  Cartesian_less_x_2(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_2 first_argument_type;
  typedef typename KK::Point_2 second_argument_type;
  typedef typename KK::Motion_function::NT NT;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b) const
  {
    return a.x() - b.x();
  }
  result_type operator()(const typename result_type::NT &c, const second_argument_type &a ) const {
    return result_type(c) - a.x();
  }
  result_type operator()( const first_argument_type &b, const typename result_type::NT &c) const {
    return b.x() - result_type(c);
  }
};

template <class KK>
struct Cartesian_less_y_2
{
  Cartesian_less_y_2(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_2 first_argument_type;
  typedef typename KK::Point_2 second_argument_type;
  typedef typename KK::Motion_function::NT NT;
 
  result_type operator()(const NT &c, const first_argument_type &a) const {
    return result_type(c) - a.y();
  }

  result_type operator()(const second_argument_type &b, const NT &c) const {
    return b.y() - result_type(c);
  }

  result_type operator()(const second_argument_type &a,
			 const first_argument_type &b) const
  {
    return a.y() - b.y();
  }

};

} } } //namespace CGAL::Kinetic::internal
#endif
