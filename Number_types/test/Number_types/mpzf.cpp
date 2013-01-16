#include <CGAL/basic.h>
#include <CGAL/mpzf.h>

using CGAL::mpzf;
int main(){
  mpzf a=3;
  mpzf b=4.5;
  mpzf c=2*(a+b)+-a*5;
  assert(CGAL::sign(c)==0);
  mpzf d=b/a;
  try {
    mpzf e=1/a; // Will throw
    return 42;
  } catch(...) {}
  try {
    mpzf e=a/0; // Will throw
    return 42;
  } catch(...) {}
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
  assert(a<b && CGAL::compare(b,a)>0);
  assert(CGAL::square(b)*4==81);
  try {
    mpzf e=CGAL::sqrt(a); // Will throw
    return 42;
  } catch(...) {}
  assert(CGAL::sqrt(mpzf(.25))*2-1==0);
  assert(CGAL::is_zero(c));
  assert(!CGAL::is_zero(a));
  assert(!CGAL::is_one(a));
  assert(CGAL::is_one(a-2));
  mpzf f=CGAL::gcd(a,b);
  assert(f/3+3/f>0);
  assert(CGAL::is_square(a+1));
  assert(!CGAL::is_square(a-1));
  assert(!CGAL::is_square(b));
  assert(a.is_canonical() && b.is_canonical() && c.is_canonical() && d.is_canonical() && f.is_canonical());
}
