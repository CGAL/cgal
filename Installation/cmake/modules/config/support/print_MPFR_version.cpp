// Copyright (c) 2005
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
// Author(s)     : various

// Tests if MPFR is available.

#include <iostream>
#include <cstdio>
#include <cstddef>
#include "mpfr.h"

int main()
{
#ifdef MPFR_VERSION
    std::cout << "version=" << MPFR_VERSION_MAJOR << "."
                            << MPFR_VERSION_MINOR << "."
                            << MPFR_VERSION_PATCHLEVEL << std::endl;
#else
    // MPFR versions < 2.2.0 did not have version strings
    std::cout << "version=unknown" << std::endl;
#endif

    return 0;
}
