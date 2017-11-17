// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_TYPE_SPHERE_H
#define CGAL_KD_TYPE_SPHERE_H
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
namespace CGAL {
template <class R_> class Sphere {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	Point_ c_;
	FT_ r2_;

	public:
	Sphere(Point_ const&p, FT_ const&r2): c_(p), r2_(r2) {}
	// TODO: Add a piecewise constructor?

	Point_ const& center()const{return c_;}
	FT_ const& squared_radius()const{return r2_;}
};

namespace CartesianDKernelFunctors {
template <class R_> struct Construct_sphere : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_sphere)
  typedef typename Get_type<R_, Sphere_tag>::type	result_type;
  typedef typename Get_type<R_, Point_tag>::type	Point;
  typedef typename Get_type<R_, FT_tag>::type FT;
  result_type operator()(Point const&a, FT const&b)const{
    return result_type(a,b);
  }
  // Not really needed
  result_type operator()()const{
    typename Get_functor<R_, Construct_ttag<Point_tag> >::type cp(this->kernel());
    return result_type(cp(),0);
  }
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    typename Get_functor<R_, Construct_circumcenter_tag>::type cc(this->kernel());
    typename Get_functor<R_, Squared_distance_tag>::type sd(this->kernel());

    // It should be possible to avoid copying the center by moving this code to a constructor.
    Point center = cc(f, e);
    FT const& r2 = sd(center, *f);
    return this->operator()(CGAL_MOVE(center), r2);
  }
};

template <class R_> struct Center_of_sphere : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Center_of_sphere)
  typedef typename Get_type<R_, Sphere_tag>::type	Sphere;
  // No reference because of the second overload
  typedef typename Get_type<R_, Point_tag>::type	result_type;

  result_type const& operator()(Sphere const&s)const{
    return s.center();
  }

  template<class Iter>
  result_type operator()(Iter b, Iter e)const{
    typename Get_functor<R_, Construct_ttag<Sphere_tag> >::type	cs(this->kernel());
    return operator()(cs(b,e)); // computes the radius needlessly
  }
};

template <class R_> struct Squared_radius {
  CGAL_FUNCTOR_INIT_IGNORE(Squared_radius)
  typedef typename Get_type<R_, Sphere_tag>::type	Sphere;
  typedef typename Get_type<R_, FT_tag>::type const&	result_type;
  // TODO: Is_exact?
  result_type operator()(Sphere const&s)const{
    return s.squared_radius();
  }
};

// FIXME: Move it to the generic functors, using the two above and conditional to the existence of sqrt(FT)
template<class R_> struct Point_of_sphere : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Point_of_sphere)
  typedef R_ R;
  typedef typename Get_type<R, FT_tag>::type FT;
  typedef typename Get_type<R, RT_tag>::type RT;
  typedef typename Get_type<R, Point_tag>::type Point;
  typedef typename Get_type<R, Sphere_tag>::type Sphere;
  typedef typename Get_functor<R, Construct_ttag<Point_tag> >::type CP;
  typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
  typedef typename Get_functor<R, Point_dimension_tag>::type PD;
  typedef Point result_type;
  typedef Sphere first_argument_type;
  typedef int second_argument_type;
  struct Trans : CGAL::binary_function<FT,int,FT> {
    FT const& r_; int idx; bool sgn;
    Trans (int n, FT const& r, bool b) : r_(r), idx(n), sgn(b) {}
    FT operator()(FT const&x, int i)const{
      return (i == idx) ? sgn ? x + r_ : x - r_ : x;
    }
  };
  result_type operator()(Sphere const&s, int i)const{
    CI ci(this->kernel());
    PD pd(this->kernel());
    typedef boost::counting_iterator<int,std::random_access_iterator_tag> Count;
    Point const&c = s.center();
    int d=pd(c);
    bool last = (i == d);
    FT r = sqrt(s.squared_radius());
    Trans t(last ? 0 : i, r, !last);
    return CP(this->kernel())(make_transforming_pair_iterator(ci(c,Begin_tag()),Count(0),t),make_transforming_pair_iterator(ci(c,End_tag()),Count(d),t));
  }
};
}
CGAL_KD_DEFAULT_TYPE(Sphere_tag,(CGAL::Sphere<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Sphere_tag>,(CartesianDKernelFunctors::Construct_sphere<K>),(Sphere_tag,Point_tag),(Construct_ttag<Point_tag>,Compute_point_cartesian_coordinate_tag,Squared_distance_tag,Squared_distance_to_origin_tag,Point_dimension_tag));
CGAL_KD_DEFAULT_FUNCTOR(Center_of_sphere_tag,(CartesianDKernelFunctors::Center_of_sphere<K>),(Sphere_tag,Point_tag),(Construct_ttag<Sphere_tag>));
CGAL_KD_DEFAULT_FUNCTOR(Squared_radius_tag,(CartesianDKernelFunctors::Squared_radius<K>),(Sphere_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_of_sphere_tag,(CartesianDKernelFunctors::Point_of_sphere<K>),(Sphere_tag,Point_tag),(Construct_ttag<Point_tag>, Construct_ttag<Point_cartesian_const_iterator_tag>));
} // namespace CGAL
#endif
