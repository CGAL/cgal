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
#include <CGAL/Triangulation_structural_filtering_traits.h>

namespace CGAL {

// The following is equivalent to Filtered_kernel< Simple_cartesian<double> >,
// but it's shorter in terms of template name length (for error messages, mangling...).

class Epick
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper< Simple_cartesian<double>::Base<Epick>::Type, Epick >,
#ifdef CGAL_NO_STATIC_FILTERS
               false >
#else
               true >
#endif
{};

typedef Epick Exact_predicates_inexact_constructions_kernel;

template <>
struct Triangulation_structural_filtering_traits<Epick> {
  typedef Tag_true Use_structural_filtering_tag;
};

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
