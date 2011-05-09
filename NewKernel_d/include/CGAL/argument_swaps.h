#ifndef CGAL_ARGUMENT_SWAPS_H
#define CGAL_ARGUMENT_SWAPS_H

#include <CGAL/config.h>
#include <utility>

#ifndef CGAL_CXX0X
#include <boost/preprocessor/repetition.hpp>
#include <boost/utility/result_of.hpp>
#endif

namespace CGAL {

#ifdef CGAL_CXX0X

namespace internal {

template<int,class...> struct Apply_to_last_then_rest_;

template<int d,class F,class T,class... U>
struct Apply_to_last_then_rest_<d,F,T,U...> {
	typedef typename Apply_to_last_then_rest_<d-1,F,U...,T>::result_type result_type;
	inline result_type operator()(F&&f,T&&t,U&&...u)const{
		return Apply_to_last_then_rest_<d-1,F,U...,T>()(
				std::forward<F>(f),
				std::forward<U>(u)...,
				std::forward<T>(t));
	}
};

template<class F,class T,class... U>
struct Apply_to_last_then_rest_<0,F,T,U...> {
	typedef decltype(std::declval<F>()(std::declval<T>(), std::declval<U>()...)) result_type;
	inline result_type operator()(F&&f,T&&t,U&&...u)const{
	return std::forward<F>(f)(std::forward<T>(t), std::forward<U>(u)...);
	}
};

} // namespace internal


struct Apply_to_last_then_rest {
	template<class F,class T,class...U> inline
	typename internal::Apply_to_last_then_rest_<sizeof...(U),F,T,U...>::result_type
	operator()(F&&f,T&&t,U&&...u)const{
	return internal::Apply_to_last_then_rest_<sizeof...(U),F,T,U...>()(
			std::forward<F>(f),
			std::forward<T>(t),
			std::forward<U>(u)...);
	}
};

#else // CGAL_CXX0X

struct Apply_to_last_then_rest {
#define CODE(Z,N,_) template<class F,class T,BOOST_PP_ENUM_PARAMS(N,class T)> \
	typename boost::result_of<F(T,BOOST_PP_ENUM_PARAMS(N,T))>::type \
	operator()(F const&f, BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t), T const&t) const { \
		return f(t,BOOST_PP_ENUM_PARAMS(N,t)); \
	}
	BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
};

#endif // CGAL_CXX0X

} // namespace CGAL

#endif // CGAL_ARGUMENT_SWAPS_H
