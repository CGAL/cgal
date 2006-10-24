// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : LiS
// File          : test/io.C
// LiS_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/IO/io.h>
#include <CGAL/Testsuite/assert.h>
#include <iostream>
#include <cstdlib>

int main() {
    std::cout << CGAL::oformat(5) << std::endl;
    std::cout << CGAL::oformat("Ok") << std::endl;
    int i;
    std::cin  >> CGAL::iformat(i);
    CGAL_test_assert( i == 42);
    return EXIT_SUCCESS;
}
