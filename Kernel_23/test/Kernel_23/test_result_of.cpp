// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Moeller

#include <CGAL/basic.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Simple_cartesian.h>

#include "CGAL/_Result_of_kernel.h"
#include "CGAL/_test_2.h"
#include "CGAL/_test_3.h"

#include <cassert>

#include "CGAL/Precise_numbers.h"

template<typename K>
bool test(const K& k) {
  return _test_2(k) &&  _test_3(k);
}

int main()
{
#if defined(CGAL_RESULT_OF_KERNEL)
  typedef CGAL::Result_of_cartesian< CGAL::Quotient<Precise_integer> > A;
  typedef CGAL::Result_of_homogeneous< Precise_integer, CGAL::Quotient<Precise_integer> > B;
  
  test( A() );
  test( B() );
#endif  
  return 0;
}
