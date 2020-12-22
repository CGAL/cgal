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
// Author(s)     : various

// Tests if GMP is available.

#include <iostream>
#include <cstdio>
#include <cstddef>
#include "gmp.h"

int main()
{
    std::cout << "version=" << __GNU_MP_VERSION << "."
                            << __GNU_MP_VERSION_MINOR << "."
                            << __GNU_MP_VERSION_PATCHLEVEL << std::endl;

    return 0;
}
