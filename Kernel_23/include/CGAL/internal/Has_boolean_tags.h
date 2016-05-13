// Copyright (c) 2016  GeometryFactory SARL(France).
// All rights reserved.
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERNAL_HAS_BOOLEAN_TAGS_H
#define CGAL_INTERNAL_HAS_BOOLEAN_TAGS_H

#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>

namespace CGAL{

namespace internal{

#define CGAL_HAS_XXX_MEMBER_NAMED_DEF(CLASS,MEMBER) \
template<typename T> struct CLASS { \
    struct Base { int MEMBER; }; \
    struct Derived : T, Base { }; \
\
    template<typename C, C> struct Check; \
\
    template<typename C> static char (&f(Check<int Base::*, &C::MEMBER>*))[1]; \
    template<typename C> static char (&f(...))[2]; \
\
    static bool const value = sizeof(f<Derived>(0)) == 2;\
};

CGAL_HAS_XXX_MEMBER_NAMED_DEF(Has_nested_type_Has_filtered_predicates,Has_filtered_predicates)
CGAL_HAS_XXX_MEMBER_NAMED_DEF(Has_nested_type_Has_static_filters,Has_static_filters)

#undef CGAL_HAS_XXX_MEMBER_NAMED_DEF

template<class Traits, bool has_tag = Has_nested_type_Has_filtered_predicates<Traits>::value >
struct Has_filtered_predicates;

template<class Traits>
struct Has_filtered_predicates<Traits,false>
{
  static const bool value = false;
};

template<class Traits>
struct Has_filtered_predicates<Traits,true>
{
  static const bool value = Traits::Has_filtered_predicates;
};

template<class Traits, bool has_tag = Has_nested_type_Has_static_filters<Traits>::value >
struct Has_static_filters;

template<class Traits>
struct Has_static_filters<Traits,false>
{
  static const bool value = false;
};

template<class Traits>
struct Has_static_filters<Traits,true>
{
  static const bool value = Traits::Has_static_filters;
};

} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_HAS_BOOLEAN_TAGS_H
