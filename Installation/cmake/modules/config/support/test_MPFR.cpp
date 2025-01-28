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
//
// Author(s)     : various

// Tests if MPFR is available.

#include <iostream>
#include "gmp.h"
#include "mpfr.h"

int main()
{
    mpfr_t b, p;
    mpfr_init (p);
    mpfr_init_set_str (b, "31", 10, GMP_RNDN);
    mpfr_mul_ui (p, b, 75, GMP_RNDU);          /* generate product */

    char *str = new char[50]; // 50 should be enough
    mp_exp_t exp;
    str = mpfr_get_str(str, &exp, 10, 0, p, GMP_RNDU);
    std::cout << str << " E " << exp << std::endl;
    delete[] str;

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
