// Copyright (c) 2011  INRIA Saclay Ile-de-France (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Marc Glisse


#ifndef CGAL_IS_ITERATOR_H
#define CGAL_IS_ITERATOR_H

#include <iterator>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(iterator_category)
BOOST_MPL_HAS_XXX_TRAIT_DEF(value_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(difference_type)
BOOST_MPL_HAS_XXX_TRAIT_DEF(pointer)
BOOST_MPL_HAS_XXX_TRAIT_DEF(reference)
  
//We request the type to be either a pointer or to
//provide all 5 nested types provided by iterator_traits  
template <class T> struct is_iterator_ {
	enum { value = 
          ( has_iterator_category<T>::value &&
            has_value_type<T>::value &&
            has_difference_type<T>::value &&
            has_pointer<T>::value &&
            has_reference<T>::value
          )
	       	|| boost::is_pointer<T>::value }; 
};

template <class T,class U,bool=is_iterator_<T>::value>
struct is_iterator_type_ {
	enum { value=false };
};
template <class T,class U> struct is_iterator_type_<T,U,true> :
	//boost::is_base_of<U,typename std::iterator_traits<T>::iterator_category>
	boost::is_convertible<typename std::iterator_traits<T>::iterator_category,U>
	{};

}

// NOTE: we don't want the real std::decay or functions are included
template <class T> struct is_iterator :
	internal::is_iterator_<typename boost::remove_cv<typename boost::remove_reference<T>::type>::type> {};

template <class T,class Tag> struct is_iterator_type :
	internal::is_iterator_type_<typename boost::remove_cv<typename boost::remove_reference<T>::type>::type,Tag> {};

}

#endif // CGAL_IS_ITERATOR_H
