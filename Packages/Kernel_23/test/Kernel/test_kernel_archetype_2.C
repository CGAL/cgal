// Copyright (c) 1999,2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken
 

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#ifndef CGAL_NO_DEPRECATED_CODE
#  define CGAL_NO_DEPRECATED_CODE
#endif
#define CGAL_CONCEPT_ARCHETYPE_PROVIDE_CONSTRUCTORS

#include <CGAL/Kernel_archetype.h>

// needed in kernel testsuite ...
CGAL::Vector_2_archetype operator-(const CGAL::Vector_2_archetype& v)
{ return v; }  

#include "CGAL/_test_new_2.h"

typedef CGAL::Kernel_archetype   Kernel;

int main()
{
  test_new_2( Kernel() );
  return 0;
}

