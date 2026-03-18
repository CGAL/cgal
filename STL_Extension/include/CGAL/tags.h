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
#include <type_traits>

namespace CGAL {

struct Void {};

/// \ingroup PkgSTLExtensionUtilities
/// Depending on bool value the class `Boolean_tag` indicates that something is `true` or `false` respectively.
template <bool b>
using Boolean_tag = std::bool_constant<b>;

/// \ingroup PkgSTLExtensionUtilities
/// The typedef `Tag_true` is `Boolean_tag<true>`.
typedef Boolean_tag<true>   Tag_true;

/// \ingroup PkgSTLExtensionUtilities
/// The typedef `Tag_false` is `Boolean_tag<false>`.
typedef Boolean_tag<false>  Tag_false;

// the function check_tag is deprecated since CGAL 3.3
inline bool check_tag( Tag_true)  {return true;}
inline bool check_tag( Tag_false) {return false;}

/// \ingroup PkgSTLExtensionUtilities
/// General tag indicating that none of any other possible tags is valid.
struct Null_tag {};

/// \ingroup PkgSTLExtensionUtilities
/// Class indicating the absence of a functor.
struct Null_functor {
  typedef Null_tag result_type;
  typedef Null_tag second_argument_type;
};

// For concurrency

/// \ingroup PkgSTLExtensionUtilities
/// Tag used to disable concurrency.
///
/// For example, it may be used by a user to request the sequential version of an algorithm.
struct Sequential_tag { static constexpr bool is_parallel = false; };

/// \ingroup PkgSTLExtensionUtilities
/// Tag used to enable concurrency.
///
/// For example, it may be used by a user to request the parallel version of an algorithm.
struct Parallel_tag : public Sequential_tag { static constexpr bool is_parallel = true; };

#ifdef CGAL_LINKED_WITH_TBB
/// \ingroup PkgSTLExtensionUtilities
/// This tag is a convenience typedef to `Parallel_tag` if the third party library Intel TBB has been found and linked, and to `Sequential_tag` otherwise.
typedef CGAL::Parallel_tag Parallel_if_available_tag;
#else
/// \ingroup PkgSTLExtensionUtilities
typedef CGAL::Sequential_tag Parallel_if_available_tag;
#endif

// For Surface_mesher and Mesh_3

/// \ingroup PkgSTLExtensionUtilities
/// The class `Non_manifold_tag` is a tag class used to monitor the surface meshing algorithm.
struct Non_manifold_tag {};

/// \ingroup PkgSTLExtensionUtilities
/// The class `Manifold_tag` is a tag class used to monitor the surface meshing algorithm.
struct Manifold_tag {};

/// \ingroup PkgSTLExtensionUtilities
/// The class `Manifold_with_boundary_tag` is a tag class used to monitor the surface meshing algorithm.
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

#if __cpp_lib_execution >= 201603L && defined(CGAL_LINKED_WITH_TBB)
#  include <execution>

namespace CGAL {
  constexpr auto std_execution_policy_aux(Sequential_tag)
  {
    return std::execution::seq;
  }
  constexpr auto std_execution_policy_aux(Parallel_tag)
  {
    return std::execution::par;
  }

  template <typename Tag>
  inline constexpr auto std_execution_policy = std_execution_policy_aux(Tag{});
} // namespace CGAL

#  define CGAL_MAYBE_EXEC_POLICY(Tag) CGAL::std_execution_policy<Tag>, // with the comma
#else // not CGAL_LINKED_WITH_TBB
#  define CGAL_MAYBE_EXEC_POLICY(Tag)
#endif // not CGAL_LINKED_WITH_TBB

#endif // CGAL_TAGS_H