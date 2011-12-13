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

#include <CGAL/_Result_of_kernel.h>
#include <CGAL/_test_2.h>
#include <CGAL/_test_3.h>

#include <cassert>


template<typename K>
bool test(const K& k) {
  return _test_2(k) || _test_3(k);
}

int main()
{
  typedef CGAL::Result_of_cartesian< double > A;
  A a;

  assert( test( a ) );
  return 0;
}
