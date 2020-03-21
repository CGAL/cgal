#include <CGAL/config.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Mpzf.h>
#endif
#include <iostream>
#include <stdlib.h>
#include "checked_NT.h"
#ifdef CGAL_HAS_MPZF
#include <CGAL/Gmpq.h>
using CGAL::Mpzf;
template<class NT,class IT>
void test1(){
  NT z;
  NT a=IT(3);
  NT b=4.5;
  NT c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);
  NT d=b/a;
  NT e=.0003;
  NT f=1e-90;
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(z<f && CGAL::compare(z,f)<0 && CGAL::compare(f,z)>0);
  assert(z==z && CGAL::compare(z,z)==0);
  assert(CGAL::square(b)*4==81);
  assert(CGAL::is_zero(c));
  assert(!CGAL::is_zero(a));
  assert(!CGAL::is_one(a));
  assert(CGAL::is_one(a-2));
  assert(CGAL::is_square(a+1));
  assert(!CGAL::is_square(a-1));
  assert(!CGAL::is_square(b));
  assert(e-e==0);
  assert(NT(1)/16*16==1);
}
template<class NT>
void test2(){
  NT a=3;
  NT b=4.5;
  NT c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);
  NT d=b/a;
#ifndef CGAL_NO_ASSERTIONS
  try {
    NT e=1/a; // Will throw
    exit(42);
  } catch(...) {}
#endif
  try {
    NT e=a/0; // Will throw
    exit(43);
  } catch(...) {}
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(CGAL::square(b)*4==81);
#ifndef CGAL_NO_ASSERTIONS
  try {
    NT e=CGAL::sqrt(a); // Will throw
    exit(44);
  } catch(...) {}
#endif
  assert(CGAL::sqrt(NT(4))==2);
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
  test1<CGAL::checked_NT<Mpzf,mpq_class>,mpz_class>();
#else
  test1<CGAL::checked_NT<Mpzf,CGAL::Gmpq>,CGAL::Gmpz>();
#endif
  test2<Mpzf>();
}
#else
int main(){
  std::cout << "Mpzf is not supported on this platform.\n";
}
#endif
