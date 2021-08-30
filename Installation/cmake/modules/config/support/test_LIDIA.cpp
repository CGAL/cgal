// Copyright (c) 2006
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
// Author(s)     : Andreas Meyer

// Tests if LIDIA is available.

#include <LiDIA/bigint.h>

int main() {
  ::LiDIA::bigint x = 123;
  x = x + x;
  return 0;
}
