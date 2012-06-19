// Copyright (c) 2012 Inria Saclay (France). All rights reserved.
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
// Author: Marc Glisse

#include <CGAL/basic.h>

int main(){
  unsigned int  a = 42;
  (void)CGAL::compare(a,a);
  unsigned long b = 42;
  (void)CGAL::compare(b,b);
#ifdef CGAL_USE_LONG_LONG
  unsigned long long c = 42;
  (void)CGAL::compare(c,c);
#endif
}
