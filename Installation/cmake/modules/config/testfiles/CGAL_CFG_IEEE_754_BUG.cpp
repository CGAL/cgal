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
// Author(s)     : Sylvain Pion

//| This flag is set if the compiler bugs with special features with IEEE 754
//| handling, concerning is_valid() and is_finite() testing of infinity and
//| nans.  The workaround is to use bitfield operations.
//| At least VC++, Borland and PGCC have the bug.

bool
is_valid(double d)
{ return (d == d); }

bool
is_finite(double d)
{ return (d == d) && is_valid(d-d); }

int main()
{
  double zero = 0;
  double inf = 1/zero;
  double nan = zero*inf;

  bool b = true;
  b = b &&  is_valid(inf);
  b = b && !is_valid(nan);
  b = b && !is_finite(inf);
  b = b && !is_finite(nan);

  if (!b)
    return -1;
  return 0;
}
