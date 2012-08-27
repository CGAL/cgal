// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <CGAL/basic.h>

#if defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFI) && defined(CGAL_USE_RS)

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polynomial.h>
#include <CGAL/RS/isolator_1.h>
#include "include/CGAL/_test_real_root_isolator.h"

int main(){
  typedef CGAL::Polynomial<CGAL::Gmpz>                          Polynomial_1;
  typedef CGAL::Gmpfr                                           Bound;
  typedef ::CGAL::internal::RS_real_root_isolator<Polynomial_1,Bound>
                                                                Isolator;

  // general test of concept RealRootIsolator
  CGAL::internal::test_real_root_isolator<Isolator>();
  return 0;
}
#else
int main(){
        return 0;
}
#endif
