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

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H

#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#elif defined CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#else
#  error "You need LEDA or CORE installed."
#endif

namespace CGAL {

#ifdef CGAL_USE_LEDA
typedef Simple_cartesian< leda_real >
        Exact_predicates_exact_constructions_kernel_with_sqrt;
#elif defined CGAL_USE_CORE
typedef Simple_cartesian< CORE::Expr >
        Exact_predicates_exact_constructions_kernel_with_sqrt;
#endif

} //namespace CGAL

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H
