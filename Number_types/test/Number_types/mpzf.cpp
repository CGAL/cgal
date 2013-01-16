#include <CGAL/basic.h>
#include <CGAL/mpzf.h>

using CGAL::mpzf;
int main(){
  mpzf a=3;
  mpzf b=4.5;
  mpzf c=2*(a+b)-5*a;
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
  b.print();
  assert(CGAL::to_double(b) == 4.5);
  std::pair<double,double> p=CGAL::to_interval(b);
  assert(p.first<=4.5 && p.second >= 4.5);
}
