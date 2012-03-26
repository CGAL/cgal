// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

//| If a compiler does not support std::tr1::tuple<> (from TR1)
//| CGAL_CFG_NO_TR1_TUPLE is set. 

#undef NDEBUG
#include <cassert>
#include <tr1/tuple>

int main()
{
  std::tr1::tuple<int, double> a (2, 1.0);
  assert(std::tr1::get<0>(a) == 2);
  return 0;
}
