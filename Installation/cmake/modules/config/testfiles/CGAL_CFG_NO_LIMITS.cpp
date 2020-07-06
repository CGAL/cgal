// Copyright (c) 1997
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
// Author(s)     : various

//| If a compiler has a bug in the implementation of
//| std::numeric_limits<>::denorm_min(), such as PGCC 7.1-2,
//| CGAL_CFG_NO_LIMITS is set.

#include <limits>

int main()
{
  double d = std::numeric_limits<double>::denorm_min();
  double e = (std::numeric_limits<double>::min)();
  // Note : denorm_min == min is actually not necessarily a bug.
  // So a better test should be found.
  if (d == 0 || d == e)
    return 1;
  return 0;
}
