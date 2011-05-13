#ifndef CGAL_CARTESIAN_LA_FUNCTORS_H
#define CGAL_CARTESIAN_LA_FUNCTORS_H

#include <CGAL/marcutils.h>
#include <CGAL/is_iterator.h>
#include <CGAL/argument_swaps.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>

namespace CGAL {
namespace CartesianDVectorBase {
#ifndef CGAL_CXX0X
namespace internal {
template<class R_,int dim_> struct Construct_LA_vector_;
#define CODE(Z,N,_) template<class R> struct Construct_LA_vector_<R,N> { \
	typedef typename R::Constructor Constructor; \
	typedef typename R::FT FT; \
	typedef typename R::LA_vector result_type; \
	result_type operator() \
	(BOOST_PP_ENUM_PARAMS(N,FT const& t)) const { \
	return typename Constructor::Values()(BOOST_PP_ENUM_PARAMS(N,t)); \
	} \
	result_type operator() \
	(BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(N),FT const& t)) const { \
	return typename Constructor::Values_divide()(t##N,BOOST_PP_ENUM_PARAMS(N,t)); \
	} \
	};
BOOST_PP_REPEAT_FROM_TO(2, 11, CODE, _ )
#undef CODE
}
#endif

template<class R_,class Zero_> struct Construct_LA_vector
#ifndef CGAL_CXX0X
: internal::Construct_LA_vector_<R_,R_::Default_ambient_dimension::value>
#endif
{
	typedef R_ R;
	typedef typename R::Constructor Constructor;
	typedef typename R::FT FT;
	typedef typename R::LA_vector result_type;
	typedef typename R_::Default_ambient_dimension Dimension;
	static const int dim=Dimension::value;
	result_type operator()(int d)const{
		CGAL_assertion(d==dim);
		return typename Constructor::Dimension()(d);
	}
	result_type operator()()const{
		return typename Constructor::Dimension()(dim);
	}
	result_type operator()(Zero_ const&)const{
		return typename Constructor::Dimension()(dim);
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
	typename std::enable_if<Constructible_from_each<FT,U...>::value &&
		(sizeof...(U)==dim), result_type>::type
	operator()(U&&...u)const{
		return typename Constructor::Values()(std::forward<U>(u)...);
	}
	//template<class...U,class=typename std::enable_if<Constructible_from_each<FT,U...>::value>::type,class=typename std::enable_if<(sizeof...(U)==dim+1)>::type,class=void>
	template<class...U>
	typename std::enable_if<Constructible_from_each<FT,U...>::value &&
		(sizeof...(U)==dim+1), result_type>::type
	operator()(U&&...u)const{
		return Apply_to_last_then_rest()(typename Constructor::Values_divide(),std::forward<U>(u)...);
	}
#else
	using internal::Construct_LA_vector_<R_,R::Default_ambient_dimension::value>::operator();
#endif
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter const&f,Iter const&g,Cartesian_tag)const
	{
		int d=std::distance(f,g);
		CGAL_assertion(d==dim);
		return typename Constructor::Iterator()(dim,f,g);
	}
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::bidirectional_iterator_tag>,result_type>::type operator()
		(Iter const&f,Iter g,Homogeneous_tag)const
	{
		--g;
		return this->operator()(f,g,*g);
	}
	template<class Iter> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter const&f,Iter const&g)const
	{
		return this->operator()(f,g,typename R::Rep_tag());
	}
	template<class Iter,class NT> typename boost::enable_if<is_iterator_type<Iter,std::forward_iterator_tag>,result_type>::type operator()
		(Iter const&f,Iter const&g,NT const&l)const
	{
		int d=std::distance(f,g);
		CGAL_assertion(d==dim);
		return typename Constructor::Iterator()(dim,CGAL::make_transforming_iterator(f,Divide<FT,NT>(l)),CGAL::make_transforming_iterator(g,Divide<FT,NT>(l)));
	}
};

template<class R_> struct Compute_cartesian_coordinate {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::LA_vector first_argument_type;
	typedef int second_argument_type;
	typedef Tag_true Is_exact;
#ifdef CGAL_CXX0X
	typedef decltype(std::declval<const first_argument_type>()[0]) result_type;
#else
	typedef FT const& result_type;
	// FT const& doesn't work with some LA (Eigen2 for instance) so we
	// should use plain FT or find a way to detect this.
#endif

	result_type operator()(first_argument_type const& v,int i)const{
		return v[i];
	}
};

template<class R_> struct Construct_cartesian_const_iterator {
	typedef R_ R;
	typedef typename R::LA_vector argument_type;
	typedef typename R::LA_vector_selector S_;
	typedef typename R::Point_cartesian_const_iterator result_type;
	// same as Vector

	result_type begin(argument_type const& v)const{
		return S_::vector_begin(v);
	}
	result_type end(argument_type const& v)const{
		return S_::vector_end(v);
	}
};

template<class R_> struct Construct_midpoint {
	typedef R_ R;
	typedef typename R::Point first_argument_type;
	typedef typename R::Point second_argument_type;
	typedef typename R::Point result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return (a+b)/2;
	}
};

template<class R_> struct Construct_sum_of_vectors {
	typedef R_ R;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;
	typedef typename R::Vector result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a+b;
	}
};

template<class R_> struct Construct_difference_of_vectors {
	typedef R_ R;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;
	typedef typename R::Vector result_type;

	result_type operator()(result_type const& a, result_type const& b)const{
		return a-b;
	}
};

template<class R_> struct Construct_opposite_vector {
	typedef R_ R;
	typedef typename R::Vector result_type;
	typedef typename R::Vector argument_type;

	result_type operator()(result_type const& v)const{
		return -v;
	}
};

template<class R_> struct Compute_scalar_product {
	typedef R_ R;
	typedef typename R::LA LA;
	typedef typename R::FT result_type;
	typedef typename R::Vector first_argument_type;
	typedef typename R::Vector second_argument_type;

	result_type operator()(first_argument_type const& a, second_argument_type const& b)const{
		return LA::dot_product(a,b);
	}
};


}
} // namespace CGAL
#endif // CGAL_CARTESIAN_LA_FUNCTORS_H
