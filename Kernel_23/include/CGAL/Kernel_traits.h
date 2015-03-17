// Copyright (c) 2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_TRAITS_H
#define CGAL_KERNEL_TRAITS_H

#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/if.hpp>

namespace CGAL {

namespace internal_kernel_traits{
  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_R,R,false)

  /// In case no specialization of CGAL::Kernel_traits is provided and the type used
  /// is not coming from a CGAL Kernel, we provide a dummy Kernel just to break later
  /// and only when it is used. This was introduced to avoid issues with a free function
  /// return a nested type of the kernel through Kernel_traits. Even if the function
  /// was not the one to consider, its return type should be instanciable.
  template <class T>
  struct Dummy_kernel
  {
    struct FT{};
  };

  template <class T, bool is_cgal_kernel = Has_nested_R<T>::value >
  struct Kernel_traits{
    typedef typename T::R type;
  };

  template <class T>
  struct Kernel_traits<T, false>{
    typedef Dummy_kernel<T> type;
  };
} // end of namespace internal_kernel_traits

template <class T>
struct Kernel_traits
{
  typedef typename internal_kernel_traits::Kernel_traits<T>::type Kernel;
  typedef Kernel type;
};

} // end namespace CGAL

#endif // CGAL_KERNEL_TRAITS_H
