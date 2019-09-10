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

#ifndef CGAL_ARGUMENT_SWAPS_H
#define CGAL_ARGUMENT_SWAPS_H

#include <CGAL/config.h>
#include <utility>

#ifndef CGAL_CXX11
#include <boost/preprocessor/repetition.hpp>
#include <boost/utility/result_of.hpp>
#endif

namespace CGAL {

#ifdef CGAL_CXX11

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

#else // CGAL_CXX11

struct Apply_to_last_then_rest {
#define CGAL_CODE(Z,N,_) template<class F,class T,BOOST_PP_ENUM_PARAMS(N,class T)> \
	typename boost::result_of<F(T,BOOST_PP_ENUM_PARAMS(N,T))>::type \
	operator()(F const&f, BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t), T const&t) const { \
		return f(t,BOOST_PP_ENUM_PARAMS(N,t)); \
	}
	BOOST_PP_REPEAT_FROM_TO(1,11,CGAL_CODE,_)
#undef CGAL_CODE
};

#endif // CGAL_CXX11

} // namespace CGAL

#endif // CGAL_ARGUMENT_SWAPS_H
