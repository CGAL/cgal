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
// $URL$
// $Id$
//
// Author(s)     : Sylvain Pion

//| If a compiler does not support statement expressions (a GCC extension)
//| CGAL_CFG_NO_STATEMENT_EXPRESSIONS is set. 

#undef NDEBUG
#include <cassert>

struct A {
  int* p;

  A(int i) : p(new int(i)) {}
  ~A() { delete p; }
  int value() const { return *p;}
};

int main()
{
  int i = __extension__ ({ int j = 2; j+j; });
  assert(i == 4);

  // The Intel Compiler complains with the following error:
  // "error: destructible entities are not allowed inside of a statement
  // expression"
  // See http://software.intel.com/en-us/articles/cdiag1487/
  i = __extension__ ({ A a(2); A b(3); a.value() + b.value(); });

  assert(i == 5);
  return 0;
}
