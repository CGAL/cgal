// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
//
//
// Author(s)     : Michael Hemmer <Michael.Hemmer@sophia.inria.fr>
//
// ============================================================================

// It seems that this does not matter. However, I keep it for further testing
// #define  CGAL_GMPFR_NO_REFCOUNT

#include <CGAL/config.h>

#ifdef CGAL_USE_MPFI

#include <cassert> // for the assert macro
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>

int main(){

#ifdef CGAL_GMPFR_NO_REFCOUNT
  std::cout << "CGAL_GMPFR_NO_REFCOUNT defined "<< std::endl ;
#else
  std::cout << "CGAL_GMPFR_NO_REFCOUNT undefined "<< std::endl ;
#endif

  CGAL::Gmpq A("6420587669/17179869184");
  std::cout << A << std::endl;


  {
    CGAL::Gmpfr a = CGAL::Gmpfi(A).sup();
    assert( a == CGAL::Gmpfi(A).sup());

    CGAL::Gmpz z(0);
    mpfr_get_z_exp(z.mpz(),CGAL::Gmpfi(A).sup().fr()); // this usage of a does not cause a bug


    CGAL::Gmpfr b = CGAL::Gmpfi(A).sup();
    assert( b == CGAL::Gmpfi(A).sup());

    std::cout << (a == CGAL::Gmpfi(A).sup()) << " "
              << (b == CGAL::Gmpfi(A).sup()) << std::endl;

    z+=0;

    CGAL::Gmpfr c = CGAL::Gmpfi(A).sup();
    assert( c == CGAL::Gmpfi(A).sup());

    std::cout << (a == CGAL::Gmpfi(A).sup()) << " "
              << (b == CGAL::Gmpfi(A).sup()) << " "
              << (c == CGAL::Gmpfi(A).sup()) << std::endl;

    assert( a == CGAL::Gmpfi(A).sup());
    assert( b == CGAL::Gmpfi(A).sup());
    assert( c == CGAL::Gmpfi(A).sup());

    assert( a == b);
    assert( b == c);
    assert( a == c);
  }

  return 0;
}
#else
int main(){
        return 0;
}
#endif
