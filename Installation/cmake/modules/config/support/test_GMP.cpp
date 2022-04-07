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

// Tests if GMP is available.

#include <iostream>
#include "gmp.h"

int main()
{
    mpz_t b, p;
    mpz_init (p);
    mpz_init_set_str (b, "31", 0);
    mpz_mul_ui (p, b, 75);          /* generate product */

    char *str = new char[mpz_sizeinbase(p, 10) + 2];
    str = mpz_get_str(str, 10, p);
    std::cout << str << std::endl;
    delete[] str;

    std::cout << "version=" << __GNU_MP_VERSION << "."
                            << __GNU_MP_VERSION_MINOR << "."
                            << __GNU_MP_VERSION_PATCHLEVEL << std::endl;

    return 0;
}
