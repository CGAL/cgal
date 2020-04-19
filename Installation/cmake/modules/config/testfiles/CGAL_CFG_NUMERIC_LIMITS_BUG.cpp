// Copyright (c) 2004
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

//| This flag is set if the compiler bugs with std::numeric_limits

#include <cmath>
#include <limits>

bool
is_finite(double d)
{ return (d != std::numeric_limits<double>::infinity()) && (-d != std::numeric_limits<double>::infinity()); }

int main()
{

  double zero = 0;
  double inf = 1. / zero;
  double nan = zero*inf;
  bool b = true;
  b = b && !is_finite(inf);

  (void) nan; // stop warning

  if (!b)
    return -1;
  return 0;
}
