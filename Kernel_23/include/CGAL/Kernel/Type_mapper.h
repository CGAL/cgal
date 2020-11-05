// Copyright (c) 2005, 2006
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_TYPE_MAPPER_H
#define CGAL_KERNEL_TYPE_MAPPER_H

#include <CGAL/basic.h>

#include <vector>

#include <boost/type_traits/remove_cv.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/remove.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <boost/preprocessor/facilities/expand.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

namespace CGAL {

namespace internal {

// the default implementation is required to catch the odd one-out
// object like Bbox
template<typename T, typename K1, typename K2 >
struct Type_mapper_impl {
  typedef T type;
};

template < typename T, typename K1, typename K2 >
struct Type_mapper_impl<std::vector< T >, K1, K2 > {
  typedef std::vector< typename Type_mapper_impl<T, K1, K2>::type > type;
};

template < typename T, typename K1, typename K2 >
struct Type_mapper_impl<boost::optional<T>, K1, K2 > {
  typedef boost::optional< typename Type_mapper_impl<T, K1, K2>::type > type;
};


/// The following code is equivalent to the one commented in CODE_TAG
/// except that with this one, the variant is really variant<A,B,C> and not
/// a internal obfuscated type
#define CGAL_TYPEMAP_TYPEDEFS(z, n, t) typedef typename Type_mapper_impl< t##n, K1, K2 >::type A##n;

#define CGAL_VARIANT_TYPEMAP(z, n, d) \
template< typename K1, typename K2, BOOST_PP_ENUM_PARAMS(n, class T) >            \
struct Type_mapper_impl<boost::variant<BOOST_PP_ENUM_PARAMS(n, T)>, K1, K2> {                    \
  BOOST_PP_REPEAT(n, CGAL_TYPEMAP_TYPEDEFS, T)                            \
  typedef boost::variant<BOOST_PP_ENUM_PARAMS(n, A)> type; \
};

BOOST_PP_REPEAT_FROM_TO(1, 10, CGAL_VARIANT_TYPEMAP, _)

#undef CGAL_TYPEMAP_TYPEDEFS
#undef CGAL_VARIANT_TYPEMAP

// CODE_TAG
//template<typename K1, typename K2, BOOST_VARIANT_ENUM_PARAMS(typename U) >
//struct Type_mapper_impl<boost::variant<BOOST_VARIANT_ENUM_PARAMS(U)>, K1, K2 > {
//  typedef typename boost::make_variant_over<
//    typename boost::mpl::transform<
//      typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(U)>::types,
//      Type_mapper_impl<boost::mpl::_1, K1, K2> >::type
//    >::type type;
//};

// Then we specialize for all kernel objects.
#define CGAL_Kernel_obj(X) \
  template < typename K1, typename K2 > \
  struct Type_mapper_impl < typename K1::X, K1, K2 > \
  { typedef typename K2::X type; };

#include <CGAL/Kernel/interface_macros.h>

template < typename K1, typename K2 >
struct Type_mapper_impl < typename K1::FT, K1, K2 >
{ typedef typename K2::FT type; };

} // internal

// This is a tool to obtain the K2::Point_2 from K1 and K1::Point_2.
// Similarly for other kernel types.

// TODO : add more specializations ?  Use a different mechanism ?

template < typename T, typename K1, typename K2 >
struct Type_mapper :
    internal::Type_mapper_impl< typename boost::remove_cv<
                                  typename boost::remove_reference < T >::type
                                  >::type, K1, K2 >
{ };

} //namespace CGAL

#endif // CGAL_KERNEL_TYPE_MAPPER_H
