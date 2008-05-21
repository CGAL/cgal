// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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

// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by install_cgal.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler does not support variadic templates (from C++0x)
//| CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES is set. 

#undef NDEBUG
#include <cassert>
//#include <iostream> // commented to speed up compilation.

double total = 0.0;

// increment
template < typename T >
T inc(const T& i)
{
  total += i;
  return i+T(1);
}

// print all args on std::cout
void print() {}

template < typename T, typename... Args >
void print(const T&t, const Args&... args)
{
  //std::cout << t << std::endl;
  (void) t;
  print(args...);
}

// A mixture for testing
void f() {}

template < typename... Args >
void f(const Args&... args)
{
  print(inc(args)...);
  //std::cout << std::endl;
}

int main()
{
  f();
  f(1);
  f(2,3.5);
  f(2,3.5,1u);
  assert(total == 13);
  return 0;
}
