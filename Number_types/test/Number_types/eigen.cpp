#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Lazy_exact_nt.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpzf.h>
#endif
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Testsuite/use.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Dense>

// Just check that it all compiles.
template <class NT, int s>
void check_(){
  Eigen::Matrix<NT,s,s> m(3,3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 10;
  NT d=(m+m).determinant();
  Eigen::Matrix<NT,s,1> v(3);
  v << 1, 2, 3;
  NT t=v.dot(v);
  v+=d*m*(t*v);
  int si=v.size();
  CGAL_USE(si);
}
template <class NT>
void check(){
  check_<NT,3>();
  check_<NT,Eigen::Dynamic>();
}

int main(){
#ifdef CGAL_USE_GMP
  check<CGAL::Gmpz>();
  check<CGAL::Gmpq>();
#endif
  check<CGAL::MP_Float>();
  check<CGAL::Quotient<CGAL::MP_Float> >();
  check<CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > >();
  {
    typedef CGAL::Interval_nt<true> I1;
    I1::Protector p1;
    check<I1>();
  }
  {
    typedef CGAL::Interval_nt<false> I2;
    I2::Protector p2;
    check<I2>();
  }
}

#else
#include <iostream>
int main(){
  std::cerr << "Eigen is not configured!\n";
}
#endif
