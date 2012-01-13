// Copyright (c) 2011  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Phillip Möller

//| If a compiler does not support std::copy_n (from C++0x)
//| CGAL_CFG_NO_CPP0X_COPY_N is set. 

#undef NDEBUG
#include <cassert>
#include <algorithm>

int main()
{
  int arr[] = {1, 2, 3, 4, 5 };
  int arr2[] = {0, 0, 0, 0, 0 };
  std::copy_n(arr, 3, arr2);

  assert(arr2[0] == 1);
  assert(arr2[1] == 2);
  assert(arr2[2] == 3);
  assert(arr2[3] == 0);
  assert(arr2[4] == 0);
  return 0;
}
