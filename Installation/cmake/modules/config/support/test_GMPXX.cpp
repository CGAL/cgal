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

// Tests if GMPXX is available.

#include <iostream>
#include <cstring> // needed by GMP 4.1.4 since <gmpxx.h> misses it.
#include <gmpxx.h>

int main()
{
    mpz_class a = 1;
    mpq_class b = 2/a;

    std::cout << a << std::endl; // test ABI of libgmpxx

    std::cout << "version=" << __GNU_MP_VERSION << "."
                            << __GNU_MP_VERSION_MINOR << "."
                            << __GNU_MP_VERSION_PATCHLEVEL << std::endl;

    return 0;
}
