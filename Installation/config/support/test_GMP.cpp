// Copyright (c) 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
