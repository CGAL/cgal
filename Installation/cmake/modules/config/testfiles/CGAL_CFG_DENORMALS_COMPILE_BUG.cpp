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

//| This flag is set if the compiler bugs when handling denormal values at
//| compile time.  At least PGCC 7.1-2 has the bug.
//|
//| Laurent Rineau, 2012/06/14: no supported platform has the bug now.

#undef NDEBUG
#include <cassert>

int main()
{
  double d = 5e-324;
  assert(d != 0);
  return 0;
}
