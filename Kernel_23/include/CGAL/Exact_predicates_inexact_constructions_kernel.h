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

namespace CGAL {

constexpr bool epick_use_static_filter =
#ifdef CGAL_NO_STATIC_FILTERS
    false;
#else
    true;
#endif

#ifndef CGAL_EPICK_SINGLE_PRECISION
// Here Epick is a class, and Double_precision_epick an alias to it.
class Epick;
class Single_precision_epick;
using Double_precision_epick = Epick;
#else // CGAL_EPICK_SINGLE_PRECISION
// Here, Epick is an alias to Single_precision_epick.
class Single_precision_epick;
class Double_precision_epick;
using Epick = Single_precision_epick;
#endif // CGAL_EPICK_SINGLE_PRECISION

namespace internal {

  // Basic objects, constructions, and predicates, using `double` for `FT`.
  using Double_precision_epick_base =
      Simple_cartesian<double>::Base<Double_precision_epick>::Type;

  // Add the type equality property, by changing the objects types
  using Double_precision_epick_base_with_type_equality =
      Type_equality_wrapper<Double_precision_epick_base, Double_precision_epick>;

  // Change the predicates to filtered predicates (with static filters or not)
  using Double_precision_epick_with_filtered_predicates =
      Filtered_kernel_adaptor<Double_precision_epick_base_with_type_equality, epick_use_static_filter>;

  // Same, but with objects using `float` for FT
  using Single_precision_epick_base =
      Simple_cartesian<float>::Base<Single_precision_epick>::Type;
  using Single_precision_epick_base_with_type_equality =
      Type_equality_wrapper<Single_precision_epick_base, Single_precision_epick>;
  using Single_precision_epick_with_filtered_predicates =
      Filtered_kernel_adaptor<Single_precision_epick_base_with_type_equality, epick_use_static_filter>;
};

#ifndef CGAL_EPICK_SINGLE_PRECISION
// The following is equivalent to Filtered_kernel< Simple_cartesian<double> >,
// but it's shorter in terms of template name length (for error messages, mangling...).
class Epick : public internal::Double_precision_epick_with_filtered_predicates
{};
#else // CGAL_EPICK_SINGLE_PRECISION
class Double_precision_epick : public internal::Double_precision_epick_with_filtered_predicates
{};
#endif // CGAL_EPICK_SINGLE_PRECISION

template <>
struct Triangulation_structural_filtering_traits<Double_precision_epick> {
  using Use_structural_filtering_tag = Tag_true;
};

using Exact_predicates_inexact_constructions_kernel = Epick;


// This kernel is Epick with the difference that its `FT` is `float`.
class Single_precision_epick
    : public Converting_constructions_kernel_adaptor<internal::Single_precision_epick_with_filtered_predicates,
                                                     Double_precision_epick>
{};

template <>
struct Triangulation_structural_filtering_traits<Single_precision_epick> {
  using Use_structural_filtering_tag = Tag_true;
};

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
