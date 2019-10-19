// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
