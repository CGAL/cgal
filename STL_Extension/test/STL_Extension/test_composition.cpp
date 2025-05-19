#include <CGAL/function_objects.h>
#include <CGAL/int.h>
#include <CGAL/number_utils_classes.h>
#include <boost/functional.hpp>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>

using CGAL::compose1_1;
using CGAL::compose1_2;
using CGAL::compose2_1;
using CGAL::compose2_2;
using CGAL::compare_to_less;
using std::accumulate;
using std::plus;
using std::multiplies;
using std::divides;
using std::transform;
using std::equal;

struct Myf {
  typedef double result_type;
  result_type operator()(int a, int b, int c, int d, int e) const
  { return a + (b * c - d) * e; }
};

int main()
{
  plus< int >        pl;
  std::function<int(int)> op1 = [](int i){ return i+1; };
  std::function<int(int)> op2 = [](int i){ return i * 2; };

  // compose1_2:
  int a[] = {3,5,7,2,4};
  int r = accumulate(a, a + 5, 0, compose1_2(op2, pl));
  assert(r == 248);

  // compose2_2:
  r = accumulate(a, a + 5, 0, compose2_2(pl, op1, op2));
  assert(r == 47);

  // compose1_1:
  transform(a, a + 5, a, compose1_1(op1, op2));
  int b[] = {7,11,15,5,9};
  assert(equal(a, a + 5, b));

  // compose2_1:
  transform(b, b + 5, a, compose2_1(pl, op2, op2));
  transform(b, b + 5, b, [](int i) { return 4 * i;} );
  assert(equal(a, a + 5, b));

  // compare_to_less
  bool boo = compare_to_less(CGAL::Compare<int>())(3,5);
  (void) boo;

  return 0;
}
