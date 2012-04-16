// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Installation/config/testfiles/CGAL_CFG_NO_CPP0X_DECLTYPE.cpp $
// $Id: CGAL_CFG_NO_CPP0X_DECLTYPE.cpp 43247 2008-05-21 15:34:36Z spion $
//
// Author(s)     : Sylvain Pion

//| If a compiler does not support std::isfinite() (from C++0x)
//| CGAL_CFG_NO_CPP0X_ISFINITE is set. 

#undef NDEBUG
#include <cassert>
#include <cmath>

#ifdef isfinite
#  error isfinite cannot be a macro if one want to use C++0x std::isfinite
// On Intel Compiler 12, isfinite is a macro in <cmath.h> :-(
#endif

template < typename T >
void use(T) {}

int main()
{
  double d = 1.0;
  bool b = std::isfinite(d);
  assert(b);
  return 0;
}
