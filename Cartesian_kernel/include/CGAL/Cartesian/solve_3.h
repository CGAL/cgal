// Copyright (c) 2000
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_CARTESIAN_SOLVE_3_H
#define CGAL_CARTESIAN_CARTESIAN_SOLVE_3_H

#include <CGAL/Kernel/solve.h>
#include <CGAL/Cartesian/Vector_3.h>

namespace CGAL {

namespace Cartesian_internal {

template <class R>
void solve (const VectorC3<R> &v0,
            const VectorC3<R> &v1,
            const VectorC3<R> &v2,
            const VectorC3<R> &d,
            typename R::FT &alpha, typename R::FT &beta, typename R::FT &gamma)
{
  CGAL::solve(v0.x(), v0.y(), v0.z(),
        v1.x(), v1.y(), v1.z(),
        v2.x(), v2.y(), v2.z(),
        d.x(),  d.y(),  d.z(),
        alpha, beta, gamma);
}

} // namespace Cartesian_internal

} //namespace CGAL

#endif // CGAL_CARTESIAN_CARTESIAN_SOLVE_3_H
