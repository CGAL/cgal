// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : 


#ifndef CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC
#include <iterator>
// #include <utility>

#define CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(T)                    \
namespace std {                                                        \
    template <>                                                        \
    struct iterator_traits<const T*> {                                 \
	typedef random_access_iterator_tag iterator_category;          \
	typedef T                          value_type;                 \
	typedef ptrdiff_t                  difference_type;            \
	typedef const T*                   pointer;                    \
	typedef const T&                   reference;                  \
    };                                                                 \
    template <>                                                        \
    struct iterator_traits<T*> {                                       \
	typedef random_access_iterator_tag iterator_category;          \
	typedef T                          value_type;                 \
	typedef ptrdiff_t                  difference_type;            \
	typedef T*                         pointer;                    \
	typedef T&                         reference;                  \
    };                                                                 \
}

// add more stuff accoring to taste...
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int)
namespace std {
    template <>                                                        \
    struct iterator_traits<const void*> {                                 \
	typedef random_access_iterator_tag iterator_category;          \
	typedef ptrdiff_t                  difference_type;            \
	typedef const void*                   pointer;                    \
    };                                                                 \
    template <>                                                        \
    struct iterator_traits<void*> {                                       \
	typedef random_access_iterator_tag iterator_category;          \
	typedef ptrdiff_t                  difference_type;            \
	typedef void*                         pointer;                    \
    };                                                                 \
}

  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(void*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short*)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(void**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(bool**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(float**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(double**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(int**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned int**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(signed char**)
  CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(unsigned short**)

#endif

