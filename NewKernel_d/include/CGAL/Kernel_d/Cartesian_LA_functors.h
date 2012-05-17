#ifndef CGAL_CARTESIAN_LA_FUNCTORS_H
#define CGAL_CARTESIAN_LA_FUNCTORS_H

#include <CGAL/marcutils.h>
#include <CGAL/is_iterator.h>
#include <CGAL/argument_swaps.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/store_kernel.h>

namespace CGAL {
namespace CartesianDVectorBase {
#ifndef CGAL_CXX0X
namespace internal {
template<class R_,int dim_> struct Construct_LA_vector_ {
	struct Never_use {};
	void operator()(Never_use)const;
};
#define CODE(Z,N,_) template<class R> struct Construct_LA_vector_<R,N> { \
	typedef typename R::Constructor Constructor; \
	typedef typename R::RT RT; \
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
BOOST_PP_REPEAT_FROM_TO(2, 11, CODE, _ )
#undef CODE
}
#endif

template<class R_,class Zero_> struct Construct_LA_vector
: private Store_kernel<R_>
#ifndef CGAL_CXX0X
, public internal::Construct_LA_vector_<R_,R_::Default_ambient_dimension::value>
#endif
{
	//CGAL_FUNCTOR_INIT_IGNORE(Construct_LA_vector)
	CGAL_FUNCTOR_INIT_STORE(Construct_LA_vector)
	typedef R_ R;
	typedef typename R::Constructor Constructor;
	typedef typename R::RT RT;
	typedef typename R::FT FT;
	typedef typename R::Vector_ result_type;
	typedef typename R_::Default_ambient_dimension Dimension;
	static const int static_dim=Dimension::value;
	result_type operator()(int d)const{
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		return typename Constructor::Dimension()(d);
	}
	result_type operator()()const{
		return typename Constructor::Dimension()(this->kernel().dimension());
	}
	result_type operator()(Zero_ const&)const{
		return typename Constructor::Dimension()(this->kernel().dimension());
	}
	result_type operator()(result_type const& v)const{
		return v;
	}
#ifdef CGAL_CXX0X
	result_type operator()(result_type&& v)const{
		return std::move(v);
	}
#endif
#ifdef CGAL_CXX0X
	template<class...U>
	typename std::enable_if<Constructible_from_each<RT,U...>::value &&
		(sizeof...(U)==static_dim), result_type>::type
	operator()(U&&...u)const{
		return typename Constructor::Values()(std::forward<U>(u)...);
	}
	//template<class...U,class=typename std::enable_if<Constructible_from_each<RT,U...>::value>::type,class=typename std::enable_if<(sizeof...(U)==static_dim+1)>::type,class=void>
	template<class...U>
	typename std::enable_if<Constructible_from_each<RT,U...>::value &&
		(sizeof...(U)==static_dim+1), result_type>::type
	operator()(U&&...u)const{
		return Apply_to_last_then_rest()(typename Constructor::Values_divide(),std::forward<U>(u)...);
	}
#else
	using internal::Construct_LA_vector_<R_,R::Default_ambient_dimension::value>::operator();
#endif
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,Cartesian_tag)const
	{
		int d=std::distance(f,g);
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		return typename Constructor::Iterator()(d,f,g);
	}
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::bidirectional_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,Homogeneous_tag)const
	{
		--g;
		return this->operator()(f,g,*g);
	}
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g)const
	{
		return this->operator()(f,g,typename R::Rep_tag());
	}
	template<class Iter,class NT> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter f,Iter g,NT const&l)const
	{
		int d=std::distance(f,g);
		CGAL_assertion(check_dimension_eq(d,this->kernel().dimension()));
		// RT? better be safe for now
		return typename Constructor::Iterator()(d,CGAL::make_transforming_iterator(f,Divide<FT,NT>(l)),CGAL::make_transforming_iterator(g,Divide<FT,NT>(l)));
	}
};

template<class R_> struct Compute_cartesian_coordinate {
	CGAL_FUNCTOR_INIT_IGNORE(Compute_cartesian_coordinate)
	typedef R_ R;
	typedef typename R_::RT RT;
	typedef typename R::Vector_ first_argument_type;
	typedef int second_argument_type;
	typedef Tag_true Is_exact;
#ifdef CGAL_CXX0X
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
	typedef typename R::Point first_argument_type;
	typedef typename R::Point second_argument_type;
	typedef typename R::Point result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return (a+b)/2;
	}
};

template<class R_> struct Sum_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Sum_of_vectors)
	typedef R_ R;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;
	typedef typename R::Vector result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a+b;
	}
};

template<class R_> struct Difference_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Difference_of_vectors)
	typedef R_ R;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;
	typedef typename R::Vector result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a-b;
	}
};

template<class R_> struct Opposite_vector {
	CGAL_FUNCTOR_INIT_IGNORE(Opposite_vector)
	typedef R_ R;
	typedef typename R::Vector result_type;
	typedef typename R::Vector argument_type;

	result_type operator()(result_type const& v)const{
		return -v;
	}
};

template<class R_> struct Scalar_product {
	CGAL_FUNCTOR_INIT_IGNORE(Scalar_product)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename R::RT result_type;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;

	result_type operator()(first_argument_type const& a, second_argument_type const& b)const{
		return LA::dot_product(a,b);
	}
};

template<class R_> struct Squared_distance_to_origin_stored {
	CGAL_FUNCTOR_INIT_IGNORE(Squared_distance_to_origin_stored)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename R::RT result_type;
	typedef typename R::Point argument_type;

	result_type operator()(argument_type const& a)const{
		return LA::squared_norm(a);
	}
};

template<class R_> struct Squared_distance_to_origin_via_dotprod {
	CGAL_FUNCTOR_INIT_IGNORE(Squared_distance_to_origin_via_dotprod)
	typedef R_ R;
	typedef typename R::LA_vector LA;
	typedef typename R::RT result_type;
	typedef typename R::Point argument_type;

	result_type operator()(argument_type const& a)const{
		return LA::dot_product(a,a);
	}
};

template<class R_> struct Orientation_of_vectors {
	CGAL_FUNCTOR_INIT_IGNORE(Orientation_of_vectors)
	typedef R_ R;
	typedef typename R::Vector_cartesian_const_iterator first_argument_type;
	typedef typename R::Vector_cartesian_const_iterator second_argument_type;
	typedef typename R::Orientation result_type;
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
	typedef typename R::Orientation result_type;
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


}
} // namespace CGAL
#endif // CGAL_CARTESIAN_LA_FUNCTORS_H
