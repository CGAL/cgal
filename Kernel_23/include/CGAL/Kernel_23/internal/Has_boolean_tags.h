// Copyright (c) 2016  GeometryFactory SARL(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
