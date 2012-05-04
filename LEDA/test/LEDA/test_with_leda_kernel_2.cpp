// Copyright (c) 1999  
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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 
#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/rat_leda.h>
#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/predicates_on_points_rat_leda_2.h>
#include "../Kernel_23/include/CGAL/_test_fct_point_2.h"

int
main()
{
  _test_fct_point_2( CGAL::use_rat_leda_kernel() );
  return 0;
}
#else
int main() { return 0; }
#endif // CGAL_USE_LEDA
