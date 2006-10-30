// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias

#ifndef CGAL_ROOT_OF_2_FWD_H
#define CGAL_ROOT_OF_2_FWD_H

#include <CGAL/config.h>
#include <utility>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template < typename T > class Root_of_2;

// template < typename T >
// Comparison_result compare(const Root_of_2<T>&, const Root_of_2<T>&);

// template < typename T >
// Comparison_result compare(const Root_of_2<T> &, const T &);

// template < typename T >
// Comparison_result compare(const T &, const Root_of_2<T> &);

// template < typename RT >
// Comparison_result
// compare(const Root_of_2<RT> &a,
// 	const typename Root_of_traits< RT >::RootOf_1 &b)
// 
// template < typename RT >  inline
// Comparison_result
// compare(const typename Root_of_traits< RT >::RootOf_1 &a,
// 	const Root_of_2<RT> &b)
// 
// template < typename RT >
// Comparison_result
// compare(const Root_of_2<RT> &a, const CGAL_CK_int(RT) &b)
// 
// template < typename RT >
// Comparison_result
// compare(const CGAL_CK_int(RT) &a, const Root_of_2<RT> &b)

// template < typename T >
// Sign sign(const Root_of_2<T> &);

// template < typename T >
// double to_double(const Root_of_2<T>&);

// template < typename T >
// std::pair<double, double> to_interval (const Root_of_2<T>&);

template < typename T >
bool is_valid(const Root_of_2<T>&);

template < typename T >
bool is_finite(const Root_of_2<T>&);

template < typename T >
io_Operator io_tag(const Root_of_2<T>&);

template < typename T >
Root_of_2<T> square(const Root_of_2<T>&);

template < typename T >
Root_of_2<T> inverse(const Root_of_2<T>&);

template < typename T >
Root_of_2<T> make_sqrt(const T&);

namespace CGALi {
template < typename RT,
           typename Has_sqrt = typename Number_type_traits<RT>::Has_sqrt >
struct Make_root_of_2_helper;
} // CGALi

// Template default version generating a Root_of_2<>.
template < typename RT >
typename CGALi::Make_root_of_2_helper<RT>::result_type
make_root_of_2(const RT &a, const RT &b, const RT &c, bool smaller);

//template < typename RT >
//typename CGALi::Make_root_of_2_helper<RT>::result_type
//make_root_of_2(const typename Root_of_traits< RT >::RootOf_1 &a, 
//               const typename Root_of_traits< RT >::RootOf_1 &b, 
//               const typename Root_of_traits< RT >::RootOf_1 &c);

CGAL_END_NAMESPACE

#endif // CGAL_ROOT_OF_2_FWD_H
