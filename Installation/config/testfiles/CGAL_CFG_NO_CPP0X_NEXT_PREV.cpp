// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
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
// Author(s)     : Philipp Möller

//| If a compiler does not support std::next and std::prev (from C++0x)
//| CGAL_CFG_NO_CPP0X_NEXT_PREV is set. 

#undef NDEBUG
#include <cassert>
#include <iterator>

int main()
{
  int i[] = {1, 2, 3, 4, 5};
  //single argument
  assert(*std::next(i) == 2);
  assert(*std::prev(i + 5) == 5);
  //two argument version
  assert(*std::next(i, 2) == 3);
  assert(*std::prev(i + 5, 2) == 4);
  return 0;
}
