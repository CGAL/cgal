#include <CGAL/basic.h>
#include <CGAL/mpzf.h>
#include <iostream>
#include <stdlib.h>
#include "checked_NT.h"
#include <CGAL/Gmpq.h>
#ifdef CGAL_HAS_MPZF
using CGAL::mpzf;
template<class NT,class IT>
void test1(){
  NT a=IT(3);
  NT b=4.5;
  NT c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);
  NT d=b/a;
  NT e=.0003;
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(CGAL::square(b)*4==81);
  assert(CGAL::is_zero(c));
  assert(!CGAL::is_zero(a));
  assert(!CGAL::is_one(a));
  assert(CGAL::is_one(a-2));
  assert(CGAL::is_square(a+1));
  assert(!CGAL::is_square(a-1));
  assert(!CGAL::is_square(b));
  assert(e-e==0);
}
template<class NT>
void test2(){
  NT a=3;
  NT b=4.5;
  NT c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);
  NT d=b/a;
  try {
    NT e=1/a; // Will throw
    exit(42);
  } catch(...) {}
  try {
    NT e=a/0; // Will throw
    exit(43);
  } catch(...) {}
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(CGAL::square(b)*4==81);
  try {
    NT e=CGAL::sqrt(a); // Will throw
    exit(44);
  } catch(...) {}
  assert(CGAL::sqrt(NT(.25))*2-1==0);
  assert(CGAL::is_zero(c));
  assert(!CGAL::is_zero(a));
  assert(!CGAL::is_one(a));
  assert(CGAL::is_one(a-2));
  NT f=CGAL::gcd(a,b);
  assert(f/3+3/f>0);
  assert(CGAL::is_square(a+1));
  assert(!CGAL::is_square(a-1));
  assert(!CGAL::is_square(b));
  assert(a.is_canonical() && b.is_canonical() && c.is_canonical() && d.is_canonical() && f.is_canonical());
}
int main(){
#ifdef CGAL_USE_GMPXX
  test1<CGAL::checked_NT<mpzf,mpq_class>,mpz_class>();
#else
  test1<CGAL::checked_NT<mpzf,CGAL::Gmpq>,CGAL::Gmpz>();
#endif
  test2<mpzf>();
}
#else
int main(){
  std::cout << "mpzf is not supported on this platform.\n";
}
#endif
