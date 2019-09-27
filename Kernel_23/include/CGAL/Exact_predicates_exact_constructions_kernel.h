// Copyright (c) 2003  
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas, Sylvain Pion

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/internal/Exact_type_selector.h>

#ifndef CGAL_DONT_USE_LAZY_KERNEL
#  include <CGAL/Lazy_kernel.h>
#endif

namespace CGAL {

// Epeck_ft is either Gmpq, or leda_rational, or Quotient<MP_float>
typedef internal::Exact_field_selector<double>::Type Epeck_ft;

// The following are redefined kernels instead of simple typedefs in order to shorten
// template name length (for error messages, mangling...).

#ifdef CGAL_DONT_USE_LAZY_KERNEL

// Equivalent to Filtered_kernel<Simple_cartesian<Lazy_exact_nt<Epeck_ft> > >
class Epeck
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper< Simple_cartesian<Lazy_exact_nt<Epeck_ft> >::Base<Epeck>::Type, Epeck >,
#ifdef CGAL_NO_STATIC_FILTERS
               false >
#else
               true >
#endif
{}; // end class Epeck

class Atomic_ref_counted_epeck
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper< Simple_cartesian<Lazy_exact_nt<Epeck_ft, Atomic_reference_counting_tag> >::Base<Atomic_ref_counted_epeck>::Type, Atomic_ref_counted_epeck >,
#ifdef CGAL_NO_STATIC_FILTERS
               false >
#else
               true >
#endif
{}; // end class Atomic_ref_counted_epeck

class Thread_safe_epeck
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper< Simple_cartesian<Lazy_exact_nt<Epeck_ft, Thread_safe_tag> >::Base<Thread_safe_epeck>::Type, Thread_safe_epeck >,
#ifdef CGAL_NO_STATIC_FILTERS
               false >
#else
               true >
#endif
{}; // end class Thread_safe_epeck

#else // no CGAL_DONT_USE_LAZY_KERNEL

// Equivalent to Lazy_kernel<Simple_cartesian<Epeck_ft> >
class Epeck
  : public Type_equality_wrapper<
             Lazy_kernel_base< Simple_cartesian<Epeck_ft>,
                               Simple_cartesian<Interval_nt_advanced>,
                               Sequential_tag,
	                       Cartesian_converter< Simple_cartesian<Epeck_ft>,
                                                    Simple_cartesian<Interval_nt_advanced> >,
                               Epeck>,
             Epeck >
{};

class Atomic_ref_counted_epeck
  : public Type_equality_wrapper<
             Lazy_kernel_base< Simple_cartesian<Epeck_ft>,
                               Simple_cartesian<Interval_nt_advanced>,
                               Atomic_reference_counting_tag,
	                       Cartesian_converter< Simple_cartesian<Epeck_ft>,
                                                    Simple_cartesian<Interval_nt_advanced> >,
                               Atomic_ref_counted_epeck>,
             Atomic_ref_counted_epeck >
{};

class Thread_safe_epeck
  : public Type_equality_wrapper<
             Lazy_kernel_base< Simple_cartesian<Epeck_ft>,
                               Simple_cartesian<Interval_nt_advanced>,
                               Thread_safe_tag,
	                       Cartesian_converter< Simple_cartesian<Epeck_ft>,
                                                    Simple_cartesian<Interval_nt_advanced> >,
                               Thread_safe_epeck>,
             Thread_safe_epeck >
{};
#endif // no CGAL_DONT_USE_LAZY_KERNEL

typedef Epeck Exact_predicates_exact_constructions_kernel;

template <>
struct Triangulation_structural_filtering_traits<Epeck> {
  typedef Tag_true Use_structural_filtering_tag;
};
template <>
struct Triangulation_structural_filtering_traits<Atomic_ref_counted_epeck> {
  typedef Tag_true Use_structural_filtering_tag;
};
template <>
struct Triangulation_structural_filtering_traits<Thread_safe_epeck> {
  typedef Tag_true Use_structural_filtering_tag;
};

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
