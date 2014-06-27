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
//
// Author(s)     : Marc Glisse

#ifndef CGAL_EXACTNESS_H
#define CGAL_EXACTNESS_H
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>
namespace CGAL {

#define CGAL_STRAWBERRY(Is_pretty) \
  namespace internal { \
    BOOST_MPL_HAS_XXX_TRAIT_DEF(Is_pretty) \
  } \
  template<class T,bool=internal::has_##Is_pretty<T>::value> \
  struct Is_pretty : boost::false_type {}; \
  template<class T> \
  struct Is_pretty<T,true> : T::Is_pretty {}

CGAL_STRAWBERRY(Is_exact);
CGAL_STRAWBERRY(Is_fast);
CGAL_STRAWBERRY(Is_stored);
#undef CGAL_STRAWBERRY
}
#endif // CGAL_EXACTNESS_H
