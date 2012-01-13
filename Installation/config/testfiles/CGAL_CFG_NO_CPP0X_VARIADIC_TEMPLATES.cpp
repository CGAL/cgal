// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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

//| If a compiler does not support variadic templates (from C++0x)
//| CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES is set. 

#undef NDEBUG
#include <cassert>

// It is annoying that the test passes in non-std=c++0x mode, hence
// triggering warnings all over the place.
// If GCC's non-c++0x mode finally rejects variadic templates at some point
// in some future release, we will be able to refine the version check.
#if defined __GNUC__ && (__GNUC__ == 4) && (__GNUC_MINOR__ >= 4) \
    && !defined __GXX_EXPERIMENTAL_CXX0X__
#  error GCC needs -std=c++0x to enable variadic templates without warnings
#endif

double total = 0.0;

template < typename T >
T inc(const T& i)
{
  total += i;
  return i+T(1);
}

void print() {}

template < typename T, typename... Args >
void print(const T&t, const Args&... args)
{
  (void) t;
  print(args...);
}

void f() {}

template < typename... Args >
void f(const Args&... args)
{
  print(inc(args)...);
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
