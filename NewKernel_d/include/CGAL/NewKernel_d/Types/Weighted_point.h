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

namespace CartesianDKernelFunctors {
template <class R_> struct Construct_weighted_point : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_weighted_point)
  typedef typename Get_type<R_, Weighted_point_tag>::type	result_type;
  typedef typename Get_type<R_, Point_tag>::type	Point;
  typedef typename Get_type<R_, FT_tag>::type FT;
  result_type operator()(Point const&a, FT const&b)const{
    return result_type(a,b);
  }
};

template <class R_> struct Point_drop_weight {
  CGAL_FUNCTOR_INIT_IGNORE(Point_drop_weight)
  typedef typename Get_type<R_, Weighted_point_tag>::type	argument_type;
  typedef typename Get_type<R_, Point_tag>::type const&		result_type;

  result_type operator()(argument_type const&s)const{
    return s.point();
  }
};

template <class R_> struct Point_weight {
  CGAL_FUNCTOR_INIT_IGNORE(Point_weight)
  typedef typename Get_type<R_, Weighted_point_tag>::type	argument_type;
  typedef typename Get_type<R_, FT_tag>::type const&		result_type;

  result_type operator()(argument_type const&s)const{
    return s.weight();
  }
};

}
CGAL_KD_DEFAULT_TYPE(Weighted_point_tag,(CGAL::Weighted_point<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Weighted_point_tag>,(CartesianDKernelFunctors::Construct_weighted_point<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_drop_weight_tag,(CartesianDKernelFunctors::Point_drop_weight<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_weight_tag,(CartesianDKernelFunctors::Point_weight<K>),(Weighted_point_tag,Point_tag),());
} // namespace CGAL
#endif
