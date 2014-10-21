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

template <class R_> struct Power_distance : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Power_distance)
  typedef typename Get_type<R_, Weighted_point_tag>::type	first_argument_type;
  typedef first_argument_type					second_argument_type;
  typedef typename Get_type<R_, FT_tag>::type			result_type;

  result_type operator()(first_argument_type const&a, second_argument_type const&b)const{
    typename Get_functor<R_, Point_drop_weight_tag>::type pdw(this->kernel());
    typename Get_functor<R_, Point_weight_tag>::type pw(this->kernel());
    typename Get_functor<R_, Squared_distance_tag>::type sd(this->kernel());
    return sd(pdw(a),pdw(b))-pw(a)-pw(b);
  }
};
template <class R_> struct Power_distance_to_point : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Power_distance_to_point)
  typedef typename Get_type<R_, Weighted_point_tag>::type	first_argument_type;
  typedef typename Get_type<R_, Point_tag>::type		second_argument_type;
  typedef typename Get_type<R_, FT_tag>::type			result_type;

  result_type operator()(first_argument_type const&a, second_argument_type const&b)const{
    typename Get_functor<R_, Point_drop_weight_tag>::type pdw(this->kernel());
    typename Get_functor<R_, Point_weight_tag>::type pw(this->kernel());
    typename Get_functor<R_, Squared_distance_tag>::type sd(this->kernel());
    return sd(pdw(a),b)-pw(a);
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

// Construct a point at (weighted) distance 0 from all the input
template <class R_> struct Power_center : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Power_center)
  typedef typename Get_type<R_, Weighted_point_tag>::type	WPoint;
  typedef WPoint result_type;
  typedef typename Get_type<R_, Point_tag>::type	Point;
  typedef typename Get_type<R_, FT_tag>::type FT;
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    // 2*(x-y).c == (x^2-wx^2)-(y^2-wy^2)
    typedef typename R_::LA LA;
    typedef typename LA::Square_matrix Matrix;
    typedef typename LA::Vector Vec;
    typedef typename LA::Construct_vector CVec;
    typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
    typename Get_functor<R_, Construct_ttag<Point_tag> >::type cp(this->kernel());
    typename Get_functor<R_, Point_dimension_tag>::type pd(this->kernel());
    typename Get_functor<R_, Squared_distance_to_origin_tag>::type sdo(this->kernel());
    typename Get_functor<R_, Power_distance_to_point_tag>::type pdp(this->kernel());
    typename Get_functor<R_, Point_drop_weight_tag>::type pdw(this->kernel());
    typename Get_functor<R_, Point_weight_tag>::type pw(this->kernel());
    typename Get_functor<R_, Construct_ttag<Weighted_point_tag> >::type cwp(this->kernel());

    WPoint const& wp0 = *f;
    Point const& p0 = pdw(wp0);
    int d = pd(p0);
    FT const& n0 = sdo(p0) - pw(wp0);
    Matrix m(d,d);
    Vec b = typename CVec::Dimension()(d);
    // Write the point coordinates in lines.
    int i;
    for(i=0; ++f!=e; ++i) {
      WPoint const& wp=*f;
      Point const& p=pdw(wp);
      FT const& np = sdo(p) - pw(wp);
      for(int j=0;j<d;++j) {
	m(i,j)=2*(c(p,j)-c(p0,j));
	b[i] = np - n0;
      }
    }
    CGAL_assertion (i == d);
    Vec res = typename CVec::Dimension()(d);;
    //std::cout << "Mat: " << m << "\n Vec: " << one << std::endl;
    LA::solve(res, CGAL_MOVE(m), CGAL_MOVE(b));
    //std::cout << "Sol: " << res << std::endl;
    Point center = cp(d,LA::vector_begin(res),LA::vector_end(res));
    FT const& r2 = pdp (wp0, center);
    return cwp(CGAL_MOVE(center), r2);
  }
};
}
CGAL_KD_DEFAULT_TYPE(Weighted_point_tag,(CGAL::KerD::Weighted_point<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Weighted_point_tag>,(CartesianDKernelFunctors::Construct_weighted_point<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_drop_weight_tag,(CartesianDKernelFunctors::Point_drop_weight<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Point_weight_tag,(CartesianDKernelFunctors::Point_weight<K>),(Weighted_point_tag,Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Power_side_of_power_sphere_tag,(CartesianDKernelFunctors::Power_side_of_power_sphere<K>),(Weighted_point_tag),(Power_side_of_power_sphere_raw_tag,Point_drop_weight_tag,Point_weight_tag));
CGAL_KD_DEFAULT_FUNCTOR(In_flat_power_side_of_power_sphere_tag,(CartesianDKernelFunctors::In_flat_power_side_of_power_sphere<K>),(Weighted_point_tag),(In_flat_power_side_of_power_sphere_raw_tag,Point_drop_weight_tag,Point_weight_tag));
CGAL_KD_DEFAULT_FUNCTOR(Power_distance_tag,(CartesianDKernelFunctors::Power_distance<K>),(Weighted_point_tag,Point_tag),(Squared_distance_tag,Point_drop_weight_tag,Point_weight_tag));
CGAL_KD_DEFAULT_FUNCTOR(Power_distance_to_point_tag,(CartesianDKernelFunctors::Power_distance_to_point<K>),(Weighted_point_tag,Point_tag),(Squared_distance_tag,Point_drop_weight_tag,Point_weight_tag));
CGAL_KD_DEFAULT_FUNCTOR(Power_center_tag,(CartesianDKernelFunctors::Power_center<K>),(Weighted_point_tag,Point_tag),(Compute_point_cartesian_coordinate_tag,Construct_ttag<Point_tag>,Construct_ttag<Weighted_point_tag>,Point_dimension_tag,Squared_distance_to_origin_tag,Point_drop_weight_tag,Point_weight_tag,Power_distance_to_point_tag));
} // namespace CGAL
#endif
