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
// Author(s)     : Andreas Fabri

//| If a compiler doesn't know nextafter() (or only knows _nextafter as VC++ 7.1).
//| nextafter() is part of ISO C99, but not ISO C++98 (hence <math.h> instead of <cmath>).
//| CGAL_CFG_NO_NEXTAFTER is set.

#include <math.h>

int main()
{
  double d = nextafter(0,0);
  (void) d;
  return 0;
}
