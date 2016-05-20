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
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_TYPE_WP_H
#define CGAL_KD_TYPE_WP_H
#include <CGAL/NewKernel_d/store_kernel.h>
#include <boost/iterator/counting_iterator.hpp>
namespace CGAL {
namespace KerD {
template <class R_> class Weighted_point {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	Point_ c_;
	FT_ w_;

	public:
	Weighted_point(Point_ const&p, FT_ const&w): c_(p), w_(w) {}
	// TODO: Add a piecewise constructor?

	Point_ const& point()const{return c_;}
	FT_ const& weight()const{return w_;}
};
}

namespace CartesianDKernelFunctors {
template <class R_> struct Construct_weighted_point : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_weighted_point)
  typedef typename Get_type<R_, Weighted_point_tag>::type	result_type;
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
};

template <class R_> struct Point_drop_weight {
  CGAL_FUNCTOR_INIT_IGNORE(Point_drop_weight)
  typedef typename Get_type<R_, Weighted_point_tag>::type	argument_type;
  typedef typename Get_type<R_, Point_tag>::type const&		result_type;
  // Returning a reference is fragile

  result_type operator()(argument_type const&s)const{
    return s.point();
  }
};

template <class R_> struct Point_weight {
  CGAL_FUNCTOR_INIT_IGNORE(Point_weight)
  typedef typename Get_type<R_, Weighted_point_tag>::type	argument_type;
  typedef typename Get_type<R_, FT_tag>::type			result_type;

  result_type operator()(argument_type const&s)const{
    return s.weight();
  }
};

template<class R_> struct Power_side_of_power_sphere : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Power_side_of_power_sphere)
  typedef R_ R;
  typedef typename Get_type<R, Oriented_side_tag>::type result_type;

  template<class Iter, class Pt>
    result_type operator()(Iter const& f, Iter const& e, Pt const& p0) const {
      typename Get_functor<R, Power_side_of_power_sphere_raw_tag>::type ptr(this->kernel());
      typename Get_functor<R, Point_drop_weight_tag>::type pdw(this->kernel());
      typename Get_functor<R, Point_weight_tag>::type pw(this->kernel());
      return ptr (
	  make_transforming_iterator (f, pdw),
	  make_transforming_iterator (e, pdw),
	  make_transforming_iterator (f, pw),
	  pdw (p0),
	  pw (p0));
    }
};

template<class R_> struct In_flat_power_side_of_power_sphere : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(In_flat_power_side_of_power_sphere)
  typedef R_ R;
  typedef typename Get_type<R, Oriented_side_tag>::type result_type;

  template<class Fo, class Iter, class Pt>
    result_type operator()(Fo const& fo, Iter const& f, Iter const& e, Pt const& p0) const {
      typename Get_functor<R, In_flat_power_side_of_power_sphere_raw_tag>::type ptr(this->kernel());
      typename Get_functor<R, Point_drop_weight_tag>::type pdw(this->kernel());
      typename Get_functor<R, Point_weight_tag>::type pw(this->kernel());
      return ptr (
	  fo,
	  make_transforming_iterator (f, pdw),
	  make_transforming_iterator (e, pdw),
	  make_transforming_iterator (f, pw),
	  pdw (p0),
	  pw (p0));
    }
};

}
CGAL_KD_DEFAULT_TYPE(Weighted_point_tag,(CGAL::KerD::Weighted_point<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Weighted_point_tag>,(CartesianDKernelFunctors::Construct_weighted_point<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_drop_weight_tag,(CartesianDKernelFunctors::Point_drop_weight<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_weight_tag,(CartesianDKernelFunctors::Point_weight<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Power_side_of_power_sphere_tag,(CartesianDKernelFunctors::Power_side_of_power_sphere<K>),(Weighted_point_tag),(Power_side_of_power_sphere_raw_tag,Point_drop_weight_tag,Point_weight_tag));
CGAL_KD_DEFAULT_FUNCTOR(In_flat_power_side_of_power_sphere_tag,(CartesianDKernelFunctors::In_flat_power_side_of_power_sphere<K>),(Weighted_point_tag),(In_flat_power_side_of_power_sphere_raw_tag,Point_drop_weight_tag,Point_weight_tag));
} // namespace CGAL
#endif
