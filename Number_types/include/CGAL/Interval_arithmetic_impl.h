// Copyright (c) 1999-2004
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
// Author(s)     : Sylvain Pion

namespace CGAL {

#ifdef CGAL_CFG_DENORMALS_COMPILE_BUG
// For compilers which bug on denormalized values at compile time.
// We generate CGAL_IA_MIN_DOUBLE at run time.
namespace {
double init_min_double()
{
    double d = 1;
    double e = 1;
    do {
        d = e;
        e = CGAL_IA_FORCE_TO_DOUBLE(e/2);
    } while (e != 0);
    return d;
}
} // anonymous namespace

#ifndef CGAL_HEADER_ONLY

namespace internal {
  double minimin = init_min_double();
  double& get_static_minimin()
  {
    return minimin;
  }
}

#else // CGAL_HEADER_ONLY

namespace internal {
  double& get_static_minimin()
  {
    static double minimin = init_min_double();
    return minimin;
  }
}
#endif // CGAL_HEADER_ONLY

#endif // CGAL_CFG_DENORMALS_COMPILE_BUG

#ifndef CGAL_HEADER_ONLY

#ifdef _MSC_VER
namespace {
int dummy_symbol_for_stopping_VC_linker_warning;
} // namespace
#endif

#endif // CGAL_HEADER_ONLY

} //namespace CGAL
