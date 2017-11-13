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

#ifndef CGAL_KINETIC_CARTESIAN_PREDICATES_3_H
#define CGAL_KINETIC_CARTESIAN_PREDICATES_3_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/determinant.h>


namespace CGAL { namespace Kinetic { namespace internal {



template <class KK>
struct Cartesian_side_of_oriented_sphere_3
{
  Cartesian_side_of_oriented_sphere_3(){}

  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  typedef typename KK::Point_3 third_argument_type;
  typedef typename KK::Point_3 fourth_argument_type;
  typedef typename KK::Point_3 fifth_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c,
			 const fourth_argument_type &d,
			 const fourth_argument_type &e)const
  {
    typedef typename KK::Motion_function FT;
    // We translate the points so that T becomes the origin.
    FT dpx = (a.x()) - (e.x());
    FT dpy = (a.y()) - (e.y());
    FT dpz = (a.z()) - (e.z());
    FT dpt = dpx*dpx + dpy*dpy + dpz*dpz;
    FT dqx = (b.x()) - (e.x());
    FT dqy = (b.y()) - (e.y());
    FT dqz = (b.z()) - (e.z());
    FT dqt = dqx*dqx + dqy*dqy + dqz*dqz;
    FT drx = (c.x()) - (e.x());
    FT dry = (c.y()) - (e.y());
    FT drz = (c.z()) - (e.z());
    FT drt = drx*drx + dry*dry + drz*drz;
    FT dsx = (d.x()) - (e.x());
    FT dsy = (d.y()) - (e.y());
    FT dsz = (d.z()) - (e.z());
    FT dst = dsx*dsx + dsy*dsy + dsz*dsz;
    FT ret= CGAL::determinant(dpx, dpy, dpz, dpt,
				    dqx, dqy, dqz, dqt,
				    drx, dry, drz, drt,
				    dsx, dsy, dsz, dst);
    /*  CGAL_KINETIC_MAPLE_LOG( << std::endl << std::endl;
	CGAL_KINETIC_MAPLE_LOG( << "pt3\n";
	CGAL_KINETIC_MAPLE_LOG( << "m:=matrix(5,5,[[";
	CGAL_KINETIC_MAPLE_LOG(
	<< a.x() << ", " << a.y() << ", " << a.z() << ", " << a.x()*a.x()+a.y()*a.y()+a.z()*a.z() << ", 1], ["
	<< a.x() << ", " << b.y() << ", " << b.z() << ", " << b.x()*b.x()+b.y()*b.y()+b.z()*b.z() << ", 1], ["
	<< c.x() << ", " << c.y() << ", " << c.z() << ", " << c.x()*c.x()+c.y()*c.y()+c.z()*c.z() << ", 1], ["
	<< d.x() << ", " << d.y() << ", " << d.z() << ", " << d.x()*d.x()+d.y()*d.y()+d.z()*a.z() << ", 1], ["
	<< e.x() << ", " << e.y() << ", " << e.z() << ", " << e.x()*e.x()+e.y()*e.y()+e.z()*e.z() << ", 1]]);\n";
	CGAL_KINETIC_MAPLE_LOG( << "det(m)-( " << ret << ");\n";
	::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m2:= matrix(4,4,[[";
	::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dpx << "," << dpy << "," << dpz << "," << dpt << "], [";
	::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dqx << "," << dqy << "," << dqz << "," << dqt << "], [";
	::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << drx << "," << dry << "," << drz << "," << drt << "], [";
	::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dsx << "," << dsy << "," << dsz << "," << dst << "]]);\n";
    */
    return ret;
  }
};

template <class KK>
struct Cartesian_power_side_of_oriented_power_sphere_3
{
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Weighted_point_3 first_argument_type;
  typedef typename KK::Weighted_point_3 second_argument_type;
  typedef typename KK::Weighted_point_3 third_argument_type;
  typedef typename KK::Weighted_point_3 fourth_argument_type;
  typedef typename KK::Weighted_point_3 fifth_argument_type;

  Cartesian_power_side_of_oriented_power_sphere_3(){}

  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c,
			 const fourth_argument_type &d,
			 const fourth_argument_type &e)const
  {
    typedef typename KK::Motion_function FT;
    // We translate the points so that T becomes the origin.
    FT dpx = (a.point().x()) - (e.point().x());
    FT dpy = (a.point().y()) - (e.point().y());
    FT dpz = (a.point().z()) - (e.point().z());
    FT dpt = dpx*dpx + dpy*dpy
      + dpz*dpz - (a.weight()) + (e.weight());
    FT dqx = (b.point().x()) - (e.point().x());
    FT dqy = (b.point().y()) - (e.point().y());
    FT dqz = (b.point().z()) - (e.point().z());
    FT dqt = dqx*dqx + dqy*dqy
      + dqz*dqz - (b.weight()) + (e.weight());
    FT drx = (c.point().x()) - (e.point().x());
    FT dry = (c.point().y()) - (e.point().y());
    FT drz = (c.point().z()) - (e.point().z());
    FT drt = drx*drx + dry*dry
      + drz*drz - (c.weight()) + (e.weight());
    FT dsx = (d.point().x()) - (e.point().x());
    FT dsy = (d.point().y()) - (e.point().y());
    FT dsz = (d.point().z()) - (e.point().z());
    FT dst = dsx*dsx + dsy*dsy
      + dsz*dsz - (d.weight()) + (e.weight());

    return ::CGAL::determinant(dpx, dpy, dpz, dpt,
				     dqx, dqy, dqz, dqt,
				     drx, dry, drz, drt,
				     dsx, dsy, dsz, dst);
  }

 result_type operator()(const first_argument_type &,
			const second_argument_type &,
			const third_argument_type &,
			const fourth_argument_type &)const
  {
    CGAL_error();
    return result_type();
  }
  result_type operator()(const first_argument_type &,
			 const second_argument_type &,
			 const third_argument_type &)const
  {
    CGAL_error();
    return result_type();
  }
  result_type operator()(const first_argument_type &,
			 const second_argument_type &)const
  {
    CGAL_error();
    return result_type();
  }
};

// Works for unweighted points.
template <class KK>
struct Cartesian_lifted_power_test_3
{
  Cartesian_lifted_power_test_3(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Weighted_point_3 first_argument_type;
  typedef typename KK::Weighted_point_3 second_argument_type;
  typedef typename KK::Weighted_point_3 third_argument_type;
  typedef typename KK::Weighted_point_3 fourth_argument_type;
  typedef typename KK::Weighted_point_3 fifth_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c,
			 const fourth_argument_type &d,
			 const fourth_argument_type &e)const
  {
    typedef typename KK::Motion_function FT;

    FT dpx = (a.point().x()) - (e.point().x());
    FT dpy = (a.point().y()) - (e.point().y());
    FT dpz = (a.point().z()) - (e.point().z());
    FT dpt = (a.lifted()) - (e.lifted());
    FT dqx = (b.point().x()) - (e.point().x());
    FT dqy = (b.point().y()) - (e.point().y());
    FT dqz = (b.point().z()) - (e.point().z());
    FT dqt = (b.lifted()) - (e.lifted());
    FT drx = (c.point().x()) - (e.point().x());
    FT dry = (c.point().y()) - (e.point().y());
    FT drz = (c.point().z()) - (e.point().z());
    FT drt = (c.lifted()) - (e.lifted());
    FT dsx = (d.point().x()) - (e.point().x());
    FT dsy = (d.point().y()) - (e.point().y());
    FT dsz = (d.point().z()) - (e.point().z());
    FT dst = (d.lifted()) - (e.lifted());

    return CGAL::determinant(dpx, dpy, dpz, dpt,
				   dqx, dqy, dqz, dqt,
				   drx, dry, drz, drt,
				   dsx, dsy, dsz, dst);
  }
};

template <class KK>
struct Cartesian_linear_lifted_power_test_3
{
  Cartesian_linear_lifted_power_test_3(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Weighted_point_3 first_argument_type;
  typedef typename KK::Weighted_point_3 second_argument_type;
  typedef typename KK::Weighted_point_3 third_argument_type;
  typedef typename KK::Weighted_point_3 fourth_argument_type;
  typedef typename KK::Weighted_point_3 fifth_argument_type;
  result_type operator()(first_argument_type a,
			 second_argument_type b,
			 third_argument_type c,
			 fourth_argument_type d,
			 fourth_argument_type e)const
  {
    typedef typename KK::Motion_function FT;
    typedef typename FT::NT NT;
    /*if (a.is_constant()){
      std::swap(a,e);
      std::swap(b,c);
      }*/

    //typedef typename RT::NT FT;
    // We translate the points so that T becomes the origin.
    FT safe_ret;

    {
      FT dpx = (a.point().x()) - (e.point().x());
      FT dpy = (a.point().y()) - (e.point().y());
      FT dpz = (a.point().z()) - (e.point().z());
      FT dpt = (a.lifted()) - (e.lifted());
      FT dqx = (b.point().x()) - (e.point().x());
      FT dqy = (b.point().y()) - (e.point().y());
      FT dqz = (b.point().z()) - (e.point().z());
      FT dqt = (b.lifted()) - (e.lifted());
      FT drx = (c.point().x()) - (e.point().x());
      FT dry = (c.point().y()) - (e.point().y());
      FT drz = (c.point().z()) - (e.point().z());
      FT drt = (c.lifted()) - (e.lifted());
      FT dsx = (d.point().x()) - (e.point().x());
      FT dsy = (d.point().y()) - (e.point().y());
      FT dsz = (d.point().z()) - (e.point().z());
      FT dst = (d.lifted()) - (e.lifted());

      safe_ret = CGAL::determinant(dpx, dpy, dpz, dpt,
					 dqx, dqy, dqz, dqt,
					 drx, dry, drz, drt,
					 dsx, dsy, dsz, dst);
    }
    FT ret;
    {
      if (!b.is_constant()) {
	std::swap(a,b);
	std::swap(c,d);
      }
      else if (! c.is_constant()) {
	std::swap(a,c);
	std::swap(b,d);
      }
      else if (! d.is_constant()) {
	std::swap(a,d);
	std::swap(b,c);
      }
      else if (! e.is_constant()) {
	std::swap(a,e);
	std::swap(b,c);
      }
      CGAL_precondition(b.is_constant());
      CGAL_precondition(c.is_constant());
      CGAL_precondition(d.is_constant());
      CGAL_precondition(e.is_constant());

      FT a00 = a.point().x() - e.point().x();
      FT a01 = a.point().y() - e.point().y();
      FT a02 = a.point().z() - e.point().z();
      FT a03 = a.lifted() - e.lifted();
      NT a10 = b.point().x()[0] - e.point().x()[0];
      NT a11 = b.point().y()[0] - e.point().y()[0];
      NT a12 = b.point().z()[0] - e.point().z()[0];
      NT a13 = b.lifted()[0] - e.lifted()[0];
      NT a20 = c.point().x()[0] - e.point().x()[0];
      NT a21 = c.point().y()[0] - e.point().y()[0];
      NT a22 = c.point().z()[0] - e.point().z()[0];
      NT a23 = c.lifted()[0] - e.lifted()[0];
      NT a30 = d.point().x()[0] - e.point().x()[0];
      NT a31 = d.point().y()[0] - e.point().y()[0];
      NT a32 = d.point().z()[0] - e.point().z()[0];
      NT a33 = d.lifted()[0] - e.lifted()[0];
      // First ompute the det2x2
      const FT m01 = a10*a01 - a00*a11;
      const FT m02 = a20*a01 - a00*a21;
      const FT m03 = a30*a01 - a00*a31;
      const NT m12 = a20*a11 - a10*a21;
      const NT m13 = a30*a11 - a10*a31;
      const NT m23 = a30*a21 - a20*a31;
      // Now compute the minors of rank 3
      const FT m012 = m12*a02 - m02*a12 + m01*a22;
      const FT m013 = m13*a02 - m03*a12 + m01*a32;
      const FT m023 = m23*a02 - m03*a22 + m02*a32;
      const NT m123 = m23*a12 - m13*a22 + m12*a32;
      // Now compute the minors of rank 4
      const FT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
      ret= m0123;
    }
    if (ret != safe_ret) {
      ret.print();
      safe_ret.print();
      std::cout << a << std::endl;
      std::cout << b << std::endl;
      std::cout << c << std::endl;
      std::cout << d << std::endl;
      std::cout << e << std::endl;
    }
    CGAL_exactness_postcondition(ret==safe_ret);
    return ret;
  }
};

template <class Pt>
typename Pt::Coordinate co3(const Pt &a, const Pt &b, const Pt &c, const Pt &d)
{
  //std::cout << "Computing orientation of matrix(4,4, [[" << a << ",1], [" << b << ",1], [" << c << ",1], [" << d << ",1]]);\n";
  typedef typename Pt::Coordinate RT;

  RT px= (a.x());
  RT py= (a.y());
  RT pz= (a.z());
  RT qx= (b.x());
  RT qy= (b.y());
  RT qz= (b.z());
  RT rx= (c.x());
  RT ry= (c.y());
  RT rz= (c.z());
  RT sx= (d.x());
  RT sy= (d.y());
  RT sz= (d.z());
  if (d.is_constant()) {
    std::swap(px, sx);
    std::swap(py, sy);
    std::swap(pz, sz);
    std::swap(rx, qx);
    std::swap(ry, qy);
    std::swap(rz, qz);
  }

  RT a00= qx-px;
  RT a01= rx-px;
  RT a02= sx-px;
  RT a10= qy-py;
  RT a11= ry-py;
  RT a12= sy-py;
  RT a20= qz-pz;
  RT a21= rz-pz;
  RT a22= sz-pz;
  RT ret= CGAL::determinant(a00, a01, a02,
				  a10, a11, a12,
				  a20, a21, a22);
  /*CGAL_LOG_MAPLE(std::endl << std::endl);
    CGAL_LOG_MAPLE( "co3\n");
    CGAL_LOG_MAPLE( "m:= matrix(4,4, [[1," << a.x()<<","<<a.y() <<","<<a.z() << "], [1," << b.x()<<","<<b.y()<<","<<b.z());
    CGAL_LOG_MAPLE( "], [1," << c.x() <<","<<c.y()<<","<<c.z() << "], [1,");
    CGAL_LOG_MAPLE( d.x() << ","<<d.y()<<","<<d.z() << "]]);"<< std::endl);
    CGAL_LOG_MAPLE( "det(m)- (" << ret << ");" << std::endl);
    CGAL_LOG_MAPLE( "m2:= matrix(3,3,[[,");
    CGAL_LOG_MAPLE( a00 << "," << a01 << "," << a02 << "], [");
    CGAL_LOG_MAPLE( a10 << "," << a11 << "," << a12 << "], [");
    CGAL_LOG_MAPLE( a20 << "," << a21 << "," << a22 << "]]);");*/

  //std::cout << "returning " << ret << std::endl;
  return ret;
}




template <class KK>
struct Cartesian_orientation_3
{
  Cartesian_orientation_3(){}

  typedef typename KK::Certificate_function result_type;

  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  typedef typename KK::Point_3 third_argument_type;
  typedef typename KK::Point_3 fourth_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c,
			 const fourth_argument_type &d) const
  {
    return co3(a, b, c, d);
  }
};

template <class KK>
struct Cartesian_weighted_orientation_3: public Cartesian_orientation_3<KK>
{
  Cartesian_weighted_orientation_3(){}
  typedef typename KK::Certificate_function result_type;

  typedef typename KK::Weighted_point_3 first_argument_type;
  typedef typename KK::Weighted_point_3 second_argument_type;
  typedef typename KK::Weighted_point_3 third_argument_type;
  typedef typename KK::Weighted_point_3 fourth_argument_type;
  template <class NWP>
    result_type operator()(const NWP &a,
			   const NWP &b,
			   const NWP &c,
			   const NWP &d) const {
    return Cartesian_orientation_3<KK>::operator()(a,b,c,d);
  }
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b,
			 const third_argument_type &c,
			 const fourth_argument_type &d) const
  {
    return co3(a.point(), b.point(), c.point(), d.point());
  }

};

template <class KK>
struct Cartesian_less_x_3
{
  Cartesian_less_x_3(){}
  typedef typename KK::Certificate_function result_type;
  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b) const
  {
    return a.x() - b.x();
  }
  typedef typename KK::Motion_function::NT NT;  
  result_type operator()( const NT &c, const first_argument_type &a) const {
    return result_type(c) - a.x();
  }

  result_type operator()(const second_argument_type &b, const NT &c ) const {
    return b.x() - result_type(c);
  }
};

template <class KK>
struct Cartesian_less_y_3
{
  Cartesian_less_y_3(){}
  typedef typename KK::Certificate_function result_type;

  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b) const
  {
    return a.y() - b.y();
  }
  typedef typename KK::Motion_function::NT NT;  
  result_type operator()(const NT &c, const first_argument_type &a) const {
    return result_type(c) - a.y();
  }

  result_type operator()( const second_argument_type &b, const NT &c) const {
    return b.y() - result_type(c);
  }
};

template <class KK>
struct Cartesian_less_z_3
{
  Cartesian_less_z_3(){}
  typedef typename KK::Certificate_function result_type;

  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  result_type operator()(const first_argument_type &a,
			 const second_argument_type &b) const
  {
    return a.z() - b.z();
  }
  typedef typename KK::Motion_function::NT NT;  
  result_type operator()(const NT &c, const first_argument_type &a) const {
    return result_type(c) - a.z();
  }

  result_type operator()( const second_argument_type &b, const NT &c) const {
    return b.z() - result_type(c);
  }

};


template <class KK>
struct Cartesian_equal_3
{
  Cartesian_equal_3(){}
  typedef typename KK::Certificate_function result_type;

  typedef typename KK::Point_3 first_argument_type;
  typedef typename KK::Point_3 second_argument_type;
  template <class AT>
  result_type operator()(const AT &a,
			 const AT &b) const
  {
    if (a==b) return result_type(1);
    else return result_type(-1);
  }
 

 
};

/*PREDICATE_2_BEGIN(Point_sphere_orientation_3){
  return ;
  }
  PREDICATE_2_END(Cartesian_less_z_3);*/

#if 0
template <class Pt, class CC>
typename CC::result_type co3(const Pt &a, const Pt &b, const Pt &c, const Pt &d, CC cc)
{
  //std::cout << "Computing orientation of matrix(4,4, [[" << a << ",1], [" << b << ",1], [" << c << ",1], [" << d << ",1]]);\n";
  typedef typename CC::result_type RT;

  RT px= cc(a.x());
  RT py= cc(a.y());
  RT pz= cc(a.z());
  RT qx= cc(b.x());
  RT qy= cc(b.y());
  RT qz= cc(b.z());
  RT rx= cc(c.x());
  RT ry= cc(c.y());
  RT rz= cc(c.z());
  RT sx= cc(d.x());
  RT sy= cc(d.y());
  RT sz= cc(d.z());
  if (d.is_constant()) {
    std::swap(px, sx);
    std::swap(py, sy);
    std::swap(pz, sz);
    std::swap(rx, qx);
    std::swap(ry, qy);
    std::swap(rz, qz);
  }

  RT a00= qx-px;
  RT a01= rx-px;
  RT a02= sx-px;
  RT a10= qy-py;
  RT a11= ry-py;
  RT a12= sy-py;
  RT a20= qz-pz;
  RT a21= rz-pz;
  RT a22= sz-pz;
  RT ret= CGAL::determinant(a00, a01, a02,
				  a10, a11, a12,
				  a20, a21, a22);
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << std::endl << std::endl;
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "co3\n";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m:= matrix(4,4, [[1," << a.x()<<","<<a.y() <<","<<a.z() << "], [1," << b.x()<<","<<b.y()<<","<<b.z();
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "], [1," << c.x() <<","<<c.y()<<","<<c.z() << "], [1,";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << d.x() << ","<<d.y()<<","<<d.z() << "]]);"<< std::endl;
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "det(m)- (" << ret << ");" << std::endl;
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m2:= matrix(3,3,[[,";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << a00 << "," << a01 << "," << a02 << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << a10 << "," << a11 << "," << a12 << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << a20 << "," << a21 << "," << a22 << "]]);";

  //std::cout << "returning " << ret << std::endl;
  return ret;
}


template <class Pt, class CC>
typename CC::result_type pt3(const Pt &a, const Pt &b, const Pt &c, const Pt &d, const Pt &e, CC cc)
{
  typedef typename CC::result_type FT;
  //typedef typename RT::NT FT;
  // We translate the points so that T becomes the origin.
  FT dpx = cc(a.point().x()) - cc(e.point().x());
  FT dpy = cc(a.point().y()) - cc(e.point().y());
  FT dpz = cc(a.point().z()) - cc(e.point().z());
  FT dpt = dpx*dpx + dpy*dpy
    + dpz*dpz - cc(a.weight()) + cc(e.weight());
  FT dqx = cc(b.point().x()) - cc(e.point().x());
  FT dqy = cc(b.point().y()) - cc(e.point().y());
  FT dqz = cc(b.point().z()) - cc(e.point().z());
  FT dqt = dqx*dqx + dqy*dqy
    + dqz*dqz - cc(b.weight()) + cc(e.weight());
  FT drx = cc(c.point().x()) - cc(e.point().x());
  FT dry = cc(c.point().y()) - cc(e.point().y());
  FT drz = cc(c.point().z()) - cc(e.point().z());
  FT drt = drx*drx + dry*dry
    + drz*drz - cc(c.weight()) + cc(e.weight());
  FT dsx = cc(d.point().x()) - cc(e.point().x());
  FT dsy = cc(d.point().y()) - cc(e.point().y());
  FT dsz = cc(d.point().z()) - cc(e.point().z());
  FT dst = dsx*dsx + dsy*dsy
    + dsz*dsz - cc(d.weight()) + cc(e.weight());

  return CGAL::determinant(dpx, dpy, dpz, dpt,
				 dqx, dqy, dqz, dqt,
				 drx, dry, drz, drt,
				 dsx, dsy, dsz, dst);
}


template <class C, class CC>
typename CC::result_type pt3(const Cartesian_moving_point_3<C>  &a,
			     const Cartesian_moving_point_3<C>  &b,
			     const Cartesian_moving_point_3<C>  &c,
			     const Cartesian_moving_point_3<C>  &d,
			     const Cartesian_moving_point_3<C>  &e, CC cc)
{
  typedef typename CC::result_type FT;
  //typedef typename RT::NT FT;
  // We translate the points so that T becomes the origin.
  FT dpx = cc(a.x()) - cc(e.x());
  FT dpy = cc(a.y()) - cc(e.y());
  FT dpz = cc(a.z()) - cc(e.z());
  FT dpt = dpx*dpx + dpy*dpy + dpz*dpz;
  FT dqx = cc(b.x()) - cc(e.x());
  FT dqy = cc(b.y()) - cc(e.y());
  FT dqz = cc(b.z()) - cc(e.z());
  FT dqt = dqx*dqx + dqy*dqy + dqz*dqz;
  FT drx = cc(c.x()) - cc(e.x());
  FT dry = cc(c.y()) - cc(e.y());
  FT drz = cc(c.z()) - cc(e.z());
  FT drt = drx*drx + dry*dry + drz*drz;
  FT dsx = cc(d.x()) - cc(e.x());
  FT dsy = cc(d.y()) - cc(e.y());
  FT dsz = cc(d.z()) - cc(e.z());
  FT dst = dsx*dsx + dsy*dsy + dsz*dsz;
  FT ret= CGAL::determinant(dpx, dpy, dpz, dpt,
				  dqx, dqy, dqz, dqt,
				  drx, dry, drz, drt,
				  dsx, dsy, dsz, dst);
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << std::endl << std::endl;
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "pt3\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "m:=matrix(5,5,[[";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE)
    << a.x() << ", " << a.y() << ", " << a.z() << ", " << a.x()*a.x()+a.y()*a.y()+a.z()*a.z() << ", 1], ["
    << a.x() << ", " << b.y() << ", " << b.z() << ", " << b.x()*b.x()+b.y()*b.y()+b.z()*b.z() << ", 1], ["
    << c.x() << ", " << c.y() << ", " << c.z() << ", " << c.x()*c.x()+c.y()*c.y()+c.z()*c.z() << ", 1], ["
    << d.x() << ", " << d.y() << ", " << d.z() << ", " << d.x()*d.x()+d.y()*d.y()+d.z()*a.z() << ", 1], ["
    << e.x() << ", " << e.y() << ", " << e.z() << ", " << e.x()*e.x()+e.y()*e.y()+e.z()*e.z() << ", 1]]);\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "det(m)-( " << ret << ");\n";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m2:= matrix(4,4,[[";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dpx << "," << dpy << "," << dpz << "," << dpt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dqx << "," << dqy << "," << dqz << "," << dqt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << drx << "," << dry << "," << drz << "," << drt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dsx << "," << dsy << "," << dsz << "," << dst << "]]);\n";

  return ret;
}


template <class NT, class CC>
typename CC::result_type pt3(Cartesian_moving_lifted_point_3<CGAL::Polynomial::Linear_polynomial<NT> > a,
			     Cartesian_moving_lifted_point_3<CGAL::Polynomial::Linear_polynomial<NT> > b,
			     Cartesian_moving_lifted_point_3<CGAL::Polynomial::Linear_polynomial<NT> > c,
			     Cartesian_moving_lifted_point_3<CGAL::Polynomial::Linear_polynomial<NT> > d,
			     Cartesian_moving_lifted_point_3<CGAL::Polynomial::Linear_polynomial<NT> > e,
			     CC cc)
{
  typedef typename CC::result_type FT;
  bool flip=false;
  bool warning_this_doesnt_really_work;
  if (!e.is_constant()) {
    flip=!flip;
    std::swap(a,e);
  }
  else if (!d.is_constant()) {
    flip=!flip;
    std::swap(d,a);
  }
  else if (!c.is_constant()) {
    std::swap(c,a);
    flip=!flip;
  }
  else if (!b.is_constant()) {
    std::swap(b,a);
    flip=!flip;
  }
  /*  if (a.is_constant()){
      std::swap(a,e);
      std::swap(b,c);
      }*/

  //typedef typename RT::NT FT;
  // We translate the points so that T becomes the origin.
  //FT m[4][4];
  FT a03 = cc(a.point().x()) - cc(e.point().x());
  FT a13 = cc(a.point().y()) - cc(e.point().y());
  FT a23 = cc(a.point().z()) - cc(e.point().z());
  FT a33 = cc(a.lifted()) - cc(e.lifted());
  FT a02 = cc(b.point().x()) - cc(e.point().x());
  FT a12 = cc(b.point().y()) - cc(e.point().y());
  FT a22 = cc(b.point().z()) - cc(e.point().z());
  FT a32 = cc(b.lifted()) - cc(e.lifted());
  FT a01 = cc(c.point().x()) - cc(e.point().x());
  FT a11 = cc(c.point().y()) - cc(e.point().y());
  FT a21 = cc(c.point().z()) - cc(e.point().z());
  FT a31 = cc(c.lifted()) - cc(e.lifted());
  FT a00 = cc(d.point().x()) - cc(e.point().x());
  FT a10 = cc(d.point().y()) - cc(e.point().y());
  FT a20 = cc(d.point().z()) - cc(e.point().z());
  FT a30 = cc(d.lifted()) - cc(e.lifted());

  FT ret;
  {
    // First compute the det2x2
    const NT m01 = a10[0]*a01[0] - a00[0]*a11[0];
    const NT m02 = a20[0]*a01[0] - a00[0]*a21[0];
    const NT m03 = a30[0]*a01[0] - a00[0]*a31[0];
    const NT m12 = a20[0]*a11[0] - a10[0]*a21[0];
    const NT m13 = a30[0]*a11[0] - a10[0]*a31[0];
    const NT m23 = a30[0]*a21[0] - a20[0]*a31[0];
    // Now compute the minors of rank 3
    const NT m012 = m12*a02[0] - m02*a12[0] + m01*a22[0];
    const NT m013 = m13*a02[0] - m03*a12[0] + m01*a32[0];
    const NT m023 = m23*a02[0] - m03*a22[0] + m02*a32[0];
    const NT m123 = m23*a12[0] - m13*a22[0] + m12*a32[0];
    // Now compute the minors of rank 4
    const FT m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
    ret= m0123;
  }
  if (flip) ret=-ret;
  /*FT ret = CGAL::determinant(dpx, dpy, dpz, dpt,
    dqx, dqy, dqz, dqt,
    drx, dry, drz, drt,
    dsx, dsy, dsz, dst);*/

#if 0
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << std::endl << std::endl;
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "pt3 lifted\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "m:=matrix(5,5,[[";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE)
    << a.point().x() << ", " << a.point().y() << ", " << a.point().z() << ", " << a.lifted() << ", 1], ["
    << b.point().x() << ", " << b.point().y() << ", " << b.point().z() << ", " << b.lifted() << ", 1], ["
    << c.point().x() << ", " << c.point().y() << ", " << c.point().z() << ", " << c.lifted() << ", 1], ["
    << d.point().x() << ", " << d.point().y() << ", " << d.point().z() << ", " << d.lifted() << ", 1], ["
    << e.point().x() << ", " << e.point().y() << ", " << e.point().z() << ", " << e.lifted() << ", 1]]);\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "det(m)-( " << ret << ");\n";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m2:= matrix(4,4,[[";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dpx << "," << dpy << "," << dpz << "," << dpt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dqx << "," << dqy << "," << dqz << "," << dqt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << drx << "," << dry << "," << drz << "," << drt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dsx << "," << dsy << "," << dsz << "," << dst << "]]);\n";
#endif
  return ret;
}


template <class C, class CC>
typename CC::result_type pt3(Cartesian_moving_lifted_point_3<C> a,
			     Cartesian_moving_lifted_point_3<C> b,
			     Cartesian_moving_lifted_point_3<C> c,
			     const Cartesian_moving_lifted_point_3<C> &d,
			     Cartesian_moving_lifted_point_3<C> e,
			     CC cc)
{
  typedef typename CC::result_type FT;

  if (a.is_constant()) {
    std::swap(a,e);
    std::swap(b,c);
  }

  //typedef typename RT::NT FT;
  // We translate the points so that T becomes the origin.
  FT dpx = cc(a.point().x()) - cc(e.point().x());
  FT dpy = cc(a.point().y()) - cc(e.point().y());
  FT dpz = cc(a.point().z()) - cc(e.point().z());
  FT dpt = cc(a.lifted()) - cc(e.lifted());
  FT dqx = cc(b.point().x()) - cc(e.point().x());
  FT dqy = cc(b.point().y()) - cc(e.point().y());
  FT dqz = cc(b.point().z()) - cc(e.point().z());
  FT dqt = cc(b.lifted()) - cc(e.lifted());
  FT drx = cc(c.point().x()) - cc(e.point().x());
  FT dry = cc(c.point().y()) - cc(e.point().y());
  FT drz = cc(c.point().z()) - cc(e.point().z());
  FT drt = cc(c.lifted()) - cc(e.lifted());
  FT dsx = cc(d.point().x()) - cc(e.point().x());
  FT dsy = cc(d.point().y()) - cc(e.point().y());
  FT dsz = cc(d.point().z()) - cc(e.point().z());
  FT dst = cc(d.lifted()) - cc(e.lifted());

  FT ret = CGAL::determinant(dpx, dpy, dpz, dpt,
				   dqx, dqy, dqz, dqt,
				   drx, dry, drz, drt,
				   dsx, dsy, dsz, dst);
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << std::endl << std::endl;
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "pt3 lifted\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "m:=matrix(5,5,[[";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE)
    << a.point().x() << ", " << a.point().y() << ", " << a.point().z() << ", " << a.lifted() << ", 1], ["
    << b.point().x() << ", " << b.point().y() << ", " << b.point().z() << ", " << b.lifted() << ", 1], ["
    << c.point().x() << ", " << c.point().y() << ", " << c.point().z() << ", " << c.lifted() << ", 1], ["
    << d.point().x() << ", " << d.point().y() << ", " << d.point().z() << ", " << d.lifted() << ", 1], ["
    << e.point().x() << ", " << e.point().y() << ", " << e.point().z() << ", " << e.lifted() << ", 1]]);\n";
  CGAL::Kinetic::log()->stream(CGAL::KDS::Log::MAPLE) << "det(m)-( " << ret << ");\n";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << "m2:= matrix(4,4,[[";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dpx << "," << dpy << "," << dpz << "," << dpt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dqx << "," << dqy << "," << dqz << "," << dqt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << drx << "," << dry << "," << drz << "," << drt << "], [";
  ::CGAL::Kinetic::log()->stream(::CGAL::KDS::Log::MAPLE) << dsx << "," << dsy << "," << dsz << "," << dst << "]]);\n";
  return ret;
}
#endif

} } } //namespace CGAL::Kinetic::internal
#endif
