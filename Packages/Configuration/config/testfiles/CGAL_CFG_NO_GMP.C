// Copyright (c) 2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : various

// CGAL_CFG_NO_GMP.C
// ---------------------------------------------------------------------
// A short test program to evaluate a machine architecture.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Tests if GMP is available.

#include <cstdio>
#ifdef __SUNPRO_CC
using std::printf;
#endif
#ifdef __BORLANDC__
#include <cstddef>
using std::size_t;
using std::printf;
#endif 

#include "gmp.h"

int main()
{
    mpz_t b, p;
    mpz_init (p);
    mpz_init_set_str (b, "31", 0);
    mpz_mul_ui (p, b, 75);          /* generate product */
    char *str = new char[mpz_sizeinbase(p, 10) + 2];
    str = mpz_get_str(str, 10, p);
    printf("%s\n", str);
    delete[] str;
    return 0;
}
