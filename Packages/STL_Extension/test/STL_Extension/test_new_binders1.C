#include <CGAL/basic.h>
#include <CGAL/functional.h>

using CGAL::bind_1;
using CGAL::bind_2;
using CGAL::bind_3;
using CGAL::bind_4;
using CGAL::bind_5;
using CGAL::compose;
using std::cout;
using std::endl;

struct F0 {
  typedef CGAL::Arity_tag< 0 > Arity;
  typedef int result_type;
  int operator()() const { return 2; }
};

struct F1 {
  typedef CGAL::Arity_tag< 1 > Arity;
  typedef int result_type;
  int operator()(int x) const { return x*2; }
};

struct F2 {
  typedef CGAL::Arity_tag< 2 > Arity;
  typedef int result_type;
  int operator()(int x, int y) const { return x+y; }
};

struct F3 {
  typedef CGAL::Arity_tag< 3 > Arity;
  typedef int result_type;
  int operator()(int x, int y, int z) const { return (x-z)*y; }
};

struct F4 {
  typedef CGAL::Arity_tag< 4 > Arity;
  typedef int result_type;
  int operator()(int x, int y, int z, int k) const
  { return x*z-y*k; }
};

struct F5 {
  typedef CGAL::Arity_tag< 5 > Arity;
  typedef int result_type;
  int operator()(int x, int y, int z, int k, int l) const
  { return (x-k)*(l-z)*y; }
};

int main()
{
  F0 f0;
  F1 f1;
  F2 f2;
  F3 f3;
  F4 f4;
  F5 f5;

  // test binders
  cout << bind_1(f1, 2)() << endl;
  cout << bind_1(f2, 2)(7) << endl;
  cout << bind_1(f3, 2)(3, 4) << endl;
  cout << bind_1(f4, 2)(3, 4, 5) << endl;
  cout << bind_1(f5, 2)(3, 4, 5, 6) << endl;
  cout << bind_2(f2, 2)(7) << endl;
  cout << bind_2(f3, 2)(3, 4) << endl;
  cout << bind_2(f4, 2)(3, 4, 5) << endl;
  cout << bind_2(f5, 2)(3, 4, 5, 6) << endl;
  cout << bind_3(f3, 2)(3, 4) << endl;
  cout << bind_3(f4, 2)(3, 4, 5) << endl;
  cout << bind_3(f5, 2)(3, 4, 5, 6) << endl;
  cout << bind_4(f4, 2)(3, 4, 5) << endl;
  cout << bind_4(f5, 2)(3, 4, 5, 6) << endl;
  cout << bind_5(f5, 2)(3, 4, 5, 6) << endl;

  // test composers
  cout << compose(f1, f0)() << endl;
  cout << compose(f1, f1)(2) << endl;
  cout << compose(f1, f2)(2, 3) << endl;
  cout << compose(f1, f3)(2, 3, 4) << endl;
  cout << compose(f1, f4)(2, 3, 4, 5) << endl;
  cout << compose(f1, f5)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f0, f0)() << endl;
  cout << compose(f2, f0, f1)(2) << endl;
  cout << compose(f2, f1, f0)(2) << endl;
  cout << compose(f2, f0, f1)(2) << endl;
  cout << compose(f2, f0, f2)(2, 3) << endl;
  cout << compose(f2, f2, f0)(2, 3) << endl;
  cout << compose(f2, f1, f1)(2, 3) << endl;
  cout << compose(f2, f3, f0)(2, 3, 4) << endl;
  cout << compose(f2, f0, f3)(2, 3, 4) << endl;
  cout << compose(f2, f2, f1)(2, 3, 4) << endl;
  cout << compose(f2, f1, f2)(2, 3, 4) << endl;
  cout << compose(f2, f4, f0)(2, 3, 4, 5) << endl;
  cout << compose(f2, f0, f4)(2, 3, 4, 5) << endl;
  cout << compose(f2, f3, f1)(2, 3, 4, 5) << endl;
  cout << compose(f2, f1, f3)(2, 3, 4, 5) << endl;
  cout << compose(f2, f2, f2)(2, 3, 4, 5) << endl;
  cout << compose(f2, f5, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f0, f5)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f4, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f1, f4)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f3, f2)(2, 3, 4, 5, 6) << endl;
  cout << compose(f2, f2, f3)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f0, f0)() << endl;
  cout << compose(f3, f1, f0, f0)(2) << endl;
  cout << compose(f3, f0, f1, f0)(2) << endl;
  cout << compose(f3, f0, f0, f1)(2) << endl;
  cout << compose(f3, f2, f0, f0)(2, 3) << endl;
  cout << compose(f3, f0, f2, f0)(2, 3) << endl;
  cout << compose(f3, f0, f0, f2)(2, 3) << endl;
  cout << compose(f3, f1, f1, f0)(2, 3) << endl;
  cout << compose(f3, f1, f0, f1)(2, 3) << endl;
  cout << compose(f3, f0, f1, f1)(2, 3) << endl;
  cout << compose(f3, f3, f0, f0)(2, 3, 4) << endl;
  cout << compose(f3, f0, f3, f0)(2, 3, 4) << endl;
  cout << compose(f3, f0, f0, f3)(2, 3, 4) << endl;
  cout << compose(f3, f2, f1, f0)(2, 3, 4) << endl;
  cout << compose(f3, f2, f0, f1)(2, 3, 4) << endl;
  cout << compose(f3, f1, f2, f0)(2, 3, 4) << endl;
  cout << compose(f3, f0, f2, f1)(2, 3, 4) << endl;
  cout << compose(f3, f1, f0, f2)(2, 3, 4) << endl;
  cout << compose(f3, f0, f1, f2)(2, 3, 4) << endl;
  cout << compose(f3, f1, f1, f1)(2, 3, 4) << endl;
  cout << compose(f3, f4, f0, f0)(2, 3, 4, 5) << endl;
  cout << compose(f3, f0, f4, f0)(2, 3, 4, 5) << endl;
  cout << compose(f3, f0, f0, f4)(2, 3, 4, 5) << endl;
  cout << compose(f3, f3, f1, f0)(2, 3, 4, 5) << endl;
  cout << compose(f3, f3, f0, f1)(2, 3, 4, 5) << endl;
  cout << compose(f3, f1, f3, f0)(2, 3, 4, 5) << endl;
  cout << compose(f3, f0, f3, f1)(2, 3, 4, 5) << endl;
  cout << compose(f3, f1, f0, f3)(2, 3, 4, 5) << endl;
  cout << compose(f3, f0, f1, f3)(2, 3, 4, 5) << endl;
  cout << compose(f3, f2, f2, f0)(2, 3, 4, 5) << endl;
  cout << compose(f3, f2, f0, f2)(2, 3, 4, 5) << endl;
  cout << compose(f3, f0, f2, f2)(2, 3, 4, 5) << endl;
  cout << compose(f3, f2, f1, f1)(2, 3, 4, 5) << endl;
  cout << compose(f3, f1, f2, f1)(2, 3, 4, 5) << endl;
  cout << compose(f3, f1, f1, f2)(2, 3, 4, 5) << endl;
  cout << compose(f3, f5, f0, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f5, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f0, f5)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f4, f1, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f4, f0, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f4, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f1, f4, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f1, f4)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f1, f0, f4)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f3, f2, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f3, f0, f2)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f3, f2)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f2, f3, f0)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f0, f2, f3)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f2, f0, f3)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f3, f1, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f1, f3, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f1, f1, f3)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f2, f2, f1)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f2, f1, f2)(2, 3, 4, 5, 6) << endl;
  cout << compose(f3, f1, f2, f2)(2, 3, 4, 5, 6) << endl;

  return 0;
}

