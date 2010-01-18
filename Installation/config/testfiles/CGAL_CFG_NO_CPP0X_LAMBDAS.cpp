// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
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

//| If a compiler does not support C++0x lambdas
//| CGAL_CFG_NO_CPP0X_LAMBDAS is set. 

#include <algorithm> 
#include <cmath> 
#include <cassert> 

int main()
{
  float f[3] = {3, 1, -2};

  std::sort(f, f+3, [](float a, float b) { return std::abs(a) < std::abs(b); });

  assert(f[0] == 1);
  assert(f[1] == -2);
  assert(f[2] == 3);

  return 0;
}
