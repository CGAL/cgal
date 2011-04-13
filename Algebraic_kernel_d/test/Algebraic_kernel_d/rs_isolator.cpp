// TODO: add copyright
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
  typedef CGAL::Polynomial<CGAL::Gmpz>                  Polynomial_1;
  typedef CGAL::Gmpfr                                   Bound;
  typedef ::CGAL::internal::RS_real_root_isolator<Polynomial_1,Bound> Isolator;
  // general test of concept RealRootIsolator
  CGAL::internal::test_real_root_isolator<Isolator>();
  return 0;
}
#else
int main(){
        return 0;
}
#endif
