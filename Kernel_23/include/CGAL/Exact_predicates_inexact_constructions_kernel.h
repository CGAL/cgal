// Copyright (c) 2003
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
// Author(s)     : Menelaos Karavelas, Sylvain Pion

#ifndef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Converting_construction.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/double.h>
#include <CGAL/float.h>

namespace CGAL {

constexpr bool epick_use_static_filter =
#ifdef CGAL_NO_STATIC_FILTERS
    false;
#else
    true;
#endif

// Here Epick is a class, and Double_precision_epick an alias to it.
class Epick;
class Single_precision_epick;
using Double_precision_epick = Epick;

namespace internal {

  // Basic objects, constructions, and predicates, using the same base class as
  // Simple_cartesian<NT>: Cartesian_base without reference counting.
  template <typename NT, typename Kernel>
  using Epick_base =
      typename Simple_cartesian<NT>::template Base<Kernel>::Type;

  // Add the type equality property, by changing the objects types
  template <typename NT, typename Kernel>
  using Epick_base_with_type_equality =
      Type_equality_wrapper<Epick_base<NT, Kernel>, Kernel>;

  // Change the predicates to filtered predicates (with static filters or not)
  template <typename NT, typename Kernel>
  using Epick_with_filtered_predicates =
      Filtered_kernel_adaptor<Epick_base_with_type_equality<NT, Kernel>, epick_use_static_filter>;
};

// The following is equivalent to Filtered_kernel< Simple_cartesian<double> >,
// but it's shorter in terms of template name length (for error messages, mangling...).
class Epick : public internal::Epick_with_filtered_predicates<double, Epick>
{};

template <>
struct Triangulation_structural_filtering_traits<Double_precision_epick> {
  using Use_structural_filtering_tag = Tag_true;
};

using Exact_predicates_inexact_constructions_kernel = Epick;


// This kernel is Epick with the difference that its `FT` is `float`.
class Single_precision_epick
    : public Converting_constructions_kernel_adaptor<
          internal::Epick_with_filtered_predicates<float, Single_precision_epick>,
          Double_precision_epick>
{};

template <>
struct Triangulation_structural_filtering_traits<Single_precision_epick> {
  using Use_structural_filtering_tag = Tag_true;
};

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
