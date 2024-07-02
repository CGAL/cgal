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

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/Number_types/internal/Exact_type_selector.h>

#ifndef CGAL_DONT_USE_LAZY_KERNEL
#  include <CGAL/Lazy_kernel.h>
#endif

namespace CGAL {

constexpr bool epeck_use_static_filter =
#ifdef CGAL_NO_STATIC_FILTERS
    false;
#else
    true;
#endif

// Epeck_ft is either Gmpq, or leda_rational, or Quotient<MP_float>
using Epeck_ft = internal::Exact_field_selector<double>::Type;

// The following are redefined kernels instead of simple typedefs in order to shorten
// template name length (for error messages, mangling...).

#ifdef CGAL_DONT_USE_LAZY_KERNEL

// Equivalent to Filtered_kernel<Simple_cartesian<Lazy_exact_nt<Epeck_ft> > >
class Epeck
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper< Simple_cartesian<Lazy_exact_nt<Epeck_ft> >::Base<Epeck>::Type, Epeck >,
               epeck_use_static_filter >
{}; // end class Epeck

#else // no CGAL_DONT_USE_LAZY_KERNEL

namespace internal {
  template <typename FT>
  using Epeck_sc = Simple_cartesian<FT>;
  using Epeck_interval = Simple_cartesian<Interval_nt_advanced>;

  template <typename FT>
  using Epeck_converter = Cartesian_converter<Epeck_sc<FT>, Epeck_interval>;

  template <typename FT, typename Kernel>
  using Epeck_lazy_base = Lazy_kernel_base<Epeck_sc<FT>, Epeck_interval, Epeck_converter<FT>, Kernel>;

  template <typename FT, typename Kernel>
  using Epeck_lazy_base_with_type_equality = Type_equality_wrapper<Epeck_lazy_base<FT, Kernel>, Kernel>;
} // namespace internal

// Equivalent to Lazy_kernel<Simple_cartesian<Epeck_ft> >
class Epeck : public internal::Epeck_lazy_base_with_type_equality<Epeck_ft, Epeck> {};

#endif // no CGAL_DONT_USE_LAZY_KERNEL

using Exact_predicates_exact_constructions_kernel = Epeck;

template <>
struct Triangulation_structural_filtering_traits<Epeck> {
  using Use_structural_filtering_tag = Tag_true;
};

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
