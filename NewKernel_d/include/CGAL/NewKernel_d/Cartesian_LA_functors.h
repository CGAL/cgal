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

#ifndef CGAL_CARTESIAN_LA_FUNCTORS_H
#define CGAL_CARTESIAN_LA_FUNCTORS_H

#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/is_iterator.h>
#include <CGAL/argument_swaps.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace CartesianDVectorBase {
#ifndef CGAL_CXX11
namespace internal {
template<class R_,class Dim_> struct Construct_LA_vector_ {
	struct Never_use {};
	void operator()(Never_use)const;
};
#define CGAL_CODE(Z,N,_) template<class R> struct Construct_LA_vector_<R,Dimension_tag<N> > { \
	typedef typename R::Constructor Constructor; \
	typedef typename Get_type<R, RT_tag>::type RT; \
	typedef typename R::Vector_ result_type; \
	result_type operator() \
	(BOOST_PP_ENUM_PARAMS(N,RT const& t)) const { \
	return typename Constructor::Values()(BOOST_PP_ENUM_PARAMS(N,t)); \
	} \
	result_type operator() \
	(BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(N),RT const& t)) const { \
	return typename Constructor::Values_divide()(t##N,BOOST_PP_ENUM_PARAMS(N,t)); \
	} \
	};
BOOST_PP_REPEAT_FROM_TO(2, 11, CGAL_CODE, _ )
#undef CGAL_CODE
}
#endif

template<class R_,class Zero_> struct Construct_LA_vector
: private Store_kernel<R_>
#ifndef CGAL_CXX11
, public internal::Construct_LA_vector_<R_,typename R_::Default_ambient_dimension>
#endif
{
	//CGAL_FUNCTOR_INIT_IGNORE(Construct_LA_vector)
	CGAL_FUNCTOR_INIT_STORE(Construct_LA_vector)
	typedef R_ R;
	typedef typename R::Constructor Constructor;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, FT_tag>::type FT;
	typedef typename R::Vector_ result_type;
	typedef typename R_::Default_ambient_dimension Dimension;
	result_type operator()(int d)const{
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		return typename Constructor::Dimension()(d);
	}
	result_type operator()()const{
		return typename Constructor::Dimension()((std::max)(0,this->kernel().dimension()));
	}
	result_type operator()(int d, Zero_ const&)const{
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		return typename Constructor::Dimension()(d);
	}
	result_type operator()(Zero_ const&)const{
	  // Makes no sense for an unknown dimension.
		return typename Constructor::Dimension()(this->kernel().dimension());
	}
	result_type operator()(result_type const& v)const{
		return v;
	}
#ifdef CGAL_CXX11
	result_type operator()(result_type&& v)const{
		return std::move(v);
	}
#endif
#ifdef CGAL_CXX11
	template<class...U>
	typename std::enable_if<Constructible_from_each<RT,U...>::value &&
		boost::is_same<Dimension_tag<sizeof...(U)>, Dimension>::value,
	  result_type>::type
	operator()(U&&...u)const{
		return typename Constructor::Values()(std::forward<U>(u)...);
	}
	//template<class...U,class=typename std::enable_if<Constructible_from_each<RT,U...>::value>::type,class=typename std::enable_if<(sizeof...(U)==static_dim+1)>::type,class=void>
	template<class...U>
	typename std::enable_if<Constructible_from_each<RT,U...>::value &&
		boost::is_same<Dimension_tag<sizeof...(U)-1>, Dimension>::value,
	  result_type>::type
	operator()(U&&...u)const{
		return Apply_to_last_then_rest()(typename Constructor::Values_divide(),std::forward<U>(u)...);
	}
#else
	using internal::Construct_LA_vector_<R_,typename R::Default_ambient_dimension>::operator();
#endif
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,Cartesian_tag t)const
	{
		return this->operator()((int)std::distance(f,g),f,g,t);
	}
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(int d,Iter f,Iter g,Cartesian_tag)const
	{
		CGAL_assertion(d==std::distance(f,g));
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		return typename Constructor::Iterator()(d,f,g);
	}
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::bidirectional_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,Homogeneous_tag)const
	{
		--g;
		return this->operator()((int)std::distance(f,g),f,g,*g);
	}
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::bidirectional_iterator_tag>,result_type>::type operator()
		(int d,Iter f,Iter g,Homogeneous_tag)const
	{
		--g;
		return this->operator()(d,f,g,*g);
	}
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g)const
	{
	  // Shouldn't it try comparing dist(f,g) to the dimension if it is known?
		return this->operator()(f,g,typename R::Rep_tag());
	}
	template<class Iter> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(int d,Iter f,Iter g)const
	{
		return this->operator()(d,f,g,typename R::Rep_tag());
	}

	// Last homogeneous coordinate given separately
	template<class Iter,class NT> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(int d,Iter f,Iter g,NT const&l)const
	{
		CGAL_assertion(d==std::distance(f,g));
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		// RT? better be safe for now
		return typename Constructor::Iterator()(d,CGAL::make_transforming_iterator(f,Divide<FT,NT>(l)),CGAL::make_transforming_iterator(g,Divide<FT,NT>(l)));
	}
	template<class Iter,class NT> inline
	  typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,NT const&l)const
	{
		return this->operator()((int)std::distance(f,g),f,g,l);
	}
};

template<class R_> struct Compute_cartesian_coordinate {
	CGAL_FUNCTOR_INIT_IGNORE(Compute_cartesian_coordinate)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename R::Vector_ first_argument_type;
	typedef int second_argument_type;
	typedef Tag_true Is_exact;
#ifdef CGAL_CXX11
	typedef decltype(std::declval<const first_argument_type>()[0]) result_type;
#else
	typedef RT const& result_type;
	// RT const& doesn't work with some LA (Eigen2 for instance) so we
	// should use plain RT or find a way to detect this.
#endif

	result_type operator()(first_argument_type const& v,int i)const{
		return v[i];
	}
};

template<class R_> struct Construct_cartesian_const_iterator {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_cartesian_const_iterator)
	typedef R_ R;
	typedef typename R::Vector_ argument_type;
	typedef typename R::LA_vector S_;
	typedef typename R::Point_cartesian_const_iterator result_type;
	// same as Vector
	typedef Tag_true Is_exact;

	result_type operator()(argument_type const& v,Begin_tag)const{
		return S_::vector_begin(v);
	}
	result_type operator()(argument_type const& v,End_tag)const{
		return S_::vector_end(v);
	}
};

template<class R_> struct Midpoint {
	CGAL_FUNCTOR_INIT_IGNORE(Midpoint)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type first_argument_type;
	typedef typename Get_type<R, Point_tag>::type second_argument_type;
	typedef typename Get_type<R, Point_tag>::type result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return (a+b)/2;
	}
};

template<class R_> struct Sum_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Sum_of_vectors)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type first_argument_type;
	typedef typename Get_type<R, Vector_tag>::type second_argument_type;
	typedef typename Get_type<R, Vector_tag>::type result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a+b;
	}
};

template<class R_> struct Difference_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Difference_of_vectors)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type first_argument_type;
	typedef typename Get_type<R, Vector_tag>::type second_argument_type;
	typedef typename Get_type<R, Vector_tag>::type result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a-b;
	}
};

template<class R_> struct Opposite_vector {
	CGAL_FUNCTOR_INIT_IGNORE(Opposite_vector)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type result_type;
	typedef typename Get_type<R, Vector_tag>::type argument_type;

	result_type operator()(result_type const& v)const{
		return -v;
	}
};

template<class R_> struct Scalar_product {
	CGAL_FUNCTOR_INIT_IGNORE(Scalar_product)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename Get_type<R, RT_tag>::type result_type;
	typedef typename Get_type<R, Vector_tag>::type first_argument_type;
	typedef typename Get_type<R, Vector_tag>::type second_argument_type;

	result_type operator()(first_argument_type const& a, second_argument_type const& b)const{
		return LA::dot_product(a,b);
	}
};

template<class R_> struct Squared_distance_to_origin_stored {
	CGAL_FUNCTOR_INIT_IGNORE(Squared_distance_to_origin_stored)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename Get_type<R, RT_tag>::type result_type;
	typedef typename Get_type<R, Point_tag>::type argument_type;

	result_type operator()(argument_type const& a)const{
		return LA::squared_norm(a);
	}
};

template<class R_> struct Squared_distance_to_origin_via_dotprod {
	CGAL_FUNCTOR_INIT_IGNORE(Squared_distance_to_origin_via_dotprod)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename Get_type<R, RT_tag>::type result_type;
	typedef typename Get_type<R, Point_tag>::type argument_type;

	result_type operator()(argument_type const& a)const{
		return LA::dot_product(a,a);
	}
};

template<class R_> struct Orientation_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Orientation_of_vectors)
	typedef R_ R;
	typedef typename R::Vector_cartesian_const_iterator first_argument_type;
	typedef typename R::Vector_cartesian_const_iterator second_argument_type;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename R::LA_vector LA;

	template<class Iter>
	result_type operator()(Iter const& f, Iter const& e) const {
		return LA::determinant_of_iterators_to_vectors(f,e);
	}
};

template<class R_> struct Orientation_of_points {
	CGAL_FUNCTOR_INIT_IGNORE(Orientation_of_points)
	typedef R_ R;
	typedef typename R::Point_cartesian_const_iterator first_argument_type;
	typedef typename R::Point_cartesian_const_iterator second_argument_type;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename R::LA_vector LA;

	template<class Iter>
	result_type operator()(Iter const& f, Iter const& e) const {
		return LA::determinant_of_iterators_to_points(f,e);
	}
};

template<class R_> struct PV_dimension {
	CGAL_FUNCTOR_INIT_IGNORE(PV_dimension)
	typedef R_ R;
	typedef typename R::Vector_ argument_type;
	typedef int result_type;
	typedef typename R::LA_vector LA;
	typedef Tag_true Is_exact;

	template<class T>
	result_type operator()(T const& v) const {
		return LA::size_of_vector(v);
	}
};

template<class R_> struct Identity_functor {
  CGAL_FUNCTOR_INIT_IGNORE(Identity_functor)
  template<class T>
  T const& operator()(T const&t) const { return t; }
};

}
} // namespace CGAL
#endif // CGAL_CARTESIAN_LA_FUNCTORS_H
