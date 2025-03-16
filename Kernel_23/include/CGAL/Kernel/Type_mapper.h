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

namespace CGAL {

template < typename T, typename K1, typename K2 >
struct Type_mapper;

namespace internal {

template < typename T, typename K1, typename K2, typename = void >
struct Type_mapper_impl {
  typedef T type;
};

template < typename K1, typename K2>
struct Type_mapper_impl<K1, K1, K2> {
  typedef K2 type;
};

// 'Rep' gets a weird partial specialization because of Return_base_tag shenanigans.
// See https://github.com/CGAL/cgal/issues/3035#issuecomment-428721414
#define CGAL_Kernel_obj(X) \
  template < typename K1, typename K2 > \
  struct Type_mapper_impl < typename K1::X, K1, K2 > \
  { typedef typename K2::X type; }; \
  template < typename K1, typename K2 > \
  struct Type_mapper_impl < typename K1::X::Rep, K1, K2 > \
  { typedef typename K2::X type; };

#include <CGAL/Kernel/interface_macros.h>

template < typename K1, typename K2 >
struct Type_mapper_impl < typename K1::FT, K1, K2 >
{ typedef typename K2::FT type; };

// This matches about anything and recursively calls Type_mapper on the template parameters
// until reaching the other cases (kernel objects, K1, FT)
template <template <typename...> class T, typename... Params, typename K1, typename K2>
struct Type_mapper_impl<T<Params...>, K1, K2,
                        std::enable_if_t<
#define CGAL_Kernel_obj(X) !std::is_same<T<Params...>, typename K1::X::Rep>::value && \
                           !std::is_same<T<Params...>, typename K1::X>::value &&
#include <CGAL/Kernel/interface_macros.h>
                           !std::is_same<T<Params...>, K1>::value && // Matches K1 directly
                           !std::is_same<T<Params...>, typename K1::FT>::value // Matches K1::FT
                        >> {
  typedef T<typename Type_mapper<Params, K1, K2>::type...> type;
};

} // internal

// This is a tool to obtain e.g. K2::Point_2 from K1 and K1::Point_2.
template < typename T, typename K1, typename K2 >
struct Type_mapper
  : internal::Type_mapper_impl<CGAL::cpp20::remove_cvref_t<T>, K1, K2 >
{ };

} // namespace CGAL

#endif // CGAL_KERNEL_TYPE_MAPPER_H
