// Copyright (c) 1999
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
// Author(s)     : Stefan Schirra


#ifndef CGAL_TAGS_H
#define CGAL_TAGS_H

#include <CGAL/IO/io_tags.h>

namespace CGAL {

struct Void {};

template <bool b>
using Boolean_tag = std::bool_constant<b>;

typedef Boolean_tag<true>   Tag_true;
typedef Boolean_tag<false>  Tag_false;

// the function check_tag is deprecated since CGAL 3.3
inline bool check_tag( Tag_true)  {return true;}
inline bool check_tag( Tag_false) {return false;}

struct Null_tag {};

struct Null_functor {
  typedef Null_tag result_type;
  typedef Null_tag second_argument_type;
};

// For concurrency
struct Sequential_tag {};
struct Parallel_tag : public Sequential_tag {};

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Parallel_if_available_tag;
#else
typedef CGAL::Sequential_tag Parallel_if_available_tag;
#endif

// For Surface_mesher and Mesh_3
struct Non_manifold_tag {};
struct Manifold_tag {};
struct Manifold_with_boundary_tag {};

// A function that asserts a specific compile time tag
// forcing its two arguments to have equal type.
template <class Base>
struct Assert_tag_class
{
    void match_compile_time_tag( const Base&) const {}
};

template <class Tag, class Derived>
inline
void
Assert_compile_time_tag( const Tag&, const Derived& b)
{
  Assert_tag_class<Tag> x;
  x.match_compile_time_tag(b);
}

// To distinguish between kernel predicates for which a division-less FT is sufficient
template <typename T>
struct Needs_FT
{
  T value;
  Needs_FT(T v) : value(v) {}
  operator T() const { return value; }
};

template <typename T>
struct Remove_needs_FT
{
  using Type = T;
};

template <typename T>
struct Remove_needs_FT<Needs_FT<T> >
{
  using Type = T;
};

} // namespace CGAL

#endif // CGAL_TAGS_H
