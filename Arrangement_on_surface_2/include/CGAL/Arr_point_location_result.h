// Copyright (c) 2013 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Efi Fogel   <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POINT_LOCATION_RESULT_H
#define CGAL_ARR_POINT_LOCATION_RESULT_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

// The macro CGAL_ARR_POINT_LOCATION_VERSION controls which version of the
// point location is used. Currently two values are supported:
// 1. Point location with CGAL::Object
// 2. Point location with boost::optional<boost::variant<...> >
// The default value is 2.

#if !defined(CGAL_ARR_POINT_LOCATION_VERSION)
#define CGAL_ARR_POINT_LOCATION_VERSION 2
#endif

#include <CGAL/Object.h>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#ifdef CGAL_CFG_BOOST_VARIANT_SWAP_BUG
#if CGAL_ARR_POINT_LOCATION_VERSION > 1
#include <CGAL/Arrangement_2/Arrangement_2_iterators.h>
// workaround for this bug:
// https://svn.boost.org/trac/boost/ticket/2839
namespace boost{ namespace detail { namespace variant {

template <class CI, class F, class MI, class V, class D, class C>
inline void move_swap(
  ::CGAL::I_Filtered_const_iterator<CI, F, MI, V, D, C>& lhs,
  ::CGAL::I_Filtered_const_iterator<CI, F, MI,  V, D, C>& rhs)
{
    ::CGAL::I_Filtered_const_iterator<CI, F, MI, V, D, C> tmp( boost::detail::variant::move(lhs) );
    lhs = boost::detail::variant::move(rhs);
    rhs = boost::detail::variant::move(tmp);
}
} } }
#endif
#endif

namespace CGAL {

template <typename Arrangement_>
struct Arr_point_location_result {
  typedef Arrangement_                                   Arrangement_2;

  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle      Face_const_handle;

#if CGAL_ARR_POINT_LOCATION_VERSION < 2
  typedef CGAL::Object                                   Type;
#else
  typedef typename boost::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>     Type;
#endif
  typedef Type                                           type;

  // This function returns either make_object() or a result_type constructor
  // to generate return values. The Object version takes a dummy template
  // argument, which is needed for the return of the other option, e.g.,
  // boost::optional<boost::variant> >.
  // In theory a one parameter variant could be returned, but this _could_
  // lead to conversion overhead, and so we rather go for the real type.
  // Overloads for empty returns are also provided.
#if CGAL_ARR_POINT_LOCATION_VERSION < 2
  template<typename T>
  static
  inline CGAL::Object make_result(T t) { return CGAL::make_object(t); }

  static
  inline CGAL::Object empty_optional_result() { return CGAL::Object(); }

  template<typename T>
  const T* assign(CGAL::Object obj) const { return CGAL::object_cast<T>(&obj); }
#else
  template<typename T>
  static
  inline Type make_result(T t) { return Type(t); }

  inline
  static
  boost::optional<Type> empty_optional_result() { return boost::optional<Type>(); }

  template<typename T>
  const T* assign(const Type& obj) const { return boost::get<T>(&obj); }
#endif // CGAL_ARR_POINT_LOCATION_VERSION < 2

  //this one is only to remove warnings in functions
  static
  inline Type default_result(){
    CGAL_error_msg("This functions should never have been called!");
    return Type();
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
