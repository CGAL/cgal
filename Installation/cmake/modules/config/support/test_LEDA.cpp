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
// Author(s)     : various

// Tests if LEDA is available.

#include <iostream>

#ifndef HAVE_LEDA4X
  #include <LEDA/system/basic.h>
  #include <LEDA/numbers/integer.h>
#else
  #include <LEDA/basic.h>
  #include <LEDA/integer.h>
#endif

#if defined(_MSC_VER)
#ifndef LEDA_STD_HEADERS
#  error "LEDA_STD_HEADERS is not defined. Please consult the troubleshooting section of the installation manual."
#endif // ! LEDA_STD_HEADERS
#endif // MSVC

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

#ifndef HAVE_LEDA4X
 #include <LEDA/internal/REDEFINE_NAMES.h>
 typedef integer I;
 #include <LEDA/internal/UNDEFINE_NAMES.h>
#else
 #include <LEDA/REDEFINE_NAMES.h>
 typedef integer I;
 #include <LEDA/UNDEFINE_NAMES.h>
#endif

int main()
{
  I a(123456);
  I b(456789);

  if (a + b != I(580245))
    return 1;
  if (a - b != I(-333333))
    return 1;

  std::cout << "version=" << __LEDA__ << std::endl;

  return 0;
}
