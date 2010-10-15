// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Sylvain Pion

//| If a compiler does not support std::tr1::array<> (from TR1)
//| CGAL_CFG_NO_TR1_ARRAY is set. 

#undef NDEBUG
#include <cassert>
#include <tr1/array>

int main()
{
  std::tr1::array<int, 3> a = { {0, 2, 4} };
  assert(a[1] == 2);
  return 0;
}
