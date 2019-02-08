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

} // namespace CGAL

#endif // CGAL_ARGUMENT_SWAPS_H
