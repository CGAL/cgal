// Copyright (c) 2011  INRIA Saclay Ile-de-France (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Marc Glisse


#ifndef CGAL_TYPE_TRAITS_IS_ITERATOR_H
#define CGAL_TYPE_TRAITS_IS_ITERATOR_H

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/logical.hpp>

#include <iterator>

namespace CGAL {
namespace internal {

BOOST_MPL_HAS_XXX_TRAIT_DEF(iterator_category)
BOOST_MPL_HAS_XXX_TRAIT_DEF(value_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(difference_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(pointer)
BOOST_MPL_HAS_XXX_TRAIT_DEF(reference)

//We request the type to be either a pointer or to
//provide all 5 nested types provided by iterator_traits
template <class T>
struct is_iterator_
  : public boost::mpl::or_<
             boost::mpl::and_<
               has_iterator_category<T>,
               has_value_type<T>,
               has_difference_type<T>,
               has_pointer<T>,
               has_reference<T> >,
             boost::is_pointer<T> >
{ };

template <class T, class U, bool = is_iterator_<T>::value>
struct is_iterator_type_
  : public boost::mpl::false_
{ };

template <class T,class U>
struct is_iterator_type_<T, U, true>
  : public //boost::is_base_of<U,typename std::iterator_traits<T>::iterator_category>
           boost::is_convertible<typename std::iterator_traits<T>::iterator_category, U>
{ };

} // namespace internal

// NOTE: we don't want the real std::decay or functions are included
template <class T>
struct is_iterator
  : public internal::is_iterator_<typename boost::remove_cv<typename boost::remove_reference<T>::type>::type>
{ };

template <class T>
CGAL_CPP17_INLINE constexpr bool is_iterator_v = is_iterator<T>::value;

template <class T, class Tag>
struct is_iterator_type
  : public internal::is_iterator_type_<typename boost::remove_cv<typename boost::remove_reference<T>::type>::type, Tag>
{ };

template <class T, class Tag>
CGAL_CPP17_INLINE constexpr bool is_iterator_type_v = is_iterator_type<T,Tag>::value;


template <class T, class U, bool = is_iterator<T>::value>
struct is_iterator_to
  : public boost::mpl::false_
{ };

template <class T, class U>
struct is_iterator_to<T, U, true>
  : public boost::is_convertible<typename std::iterator_traits<T>::value_type, U>
{ };

template <class T, class U>
CGAL_CPP17_INLINE constexpr bool is_iterator_to_v = is_iterator_to<T,U>::value;


} // namespace CGAL

#endif // CGAL_IS_ITERATOR_H
