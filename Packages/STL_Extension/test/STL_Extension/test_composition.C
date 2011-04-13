#include <CGAL/basic.h>
#include <CGAL/function_objects.h>
#include <functional>
#include <algorithm>
#include <numeric>

using CGAL::compose1_1;
using CGAL::compose1_2;
using CGAL::compose2_1;
using CGAL::compose2_2;
using std::binder1st;
using std::bind1st;
using std::accumulate;
using std::plus;
using std::multiplies;
using std::transform;
using std::equal;

int main()
{
  plus< int >        pl;
  multiplies< int >  mu;
  binder1st< plus< int > >       op1 = bind1st(pl, 1);
  binder1st< multiplies< int > > op2 = bind1st(mu, 2);

  // compose1_2:
  int a[] = {3,5,7,2,4};
  int r = accumulate(a, a + 5, 0, compose1_2(op2, pl));
  CGAL_assertion(r == 248);

  // compose2_2:
  r = accumulate(a, a + 5, 0, compose2_2(pl, op1, op2));
  CGAL_assertion(r == 47);

  // compose1_1:
  transform(a, a + 5, a, compose1_1(op1, op2));
  int b[] = {7,11,15,5,9};
  CGAL_assertion(equal(a, a + 5, b));

  // compose2_1:
  transform(b, b + 5, a, compose2_1(pl, op2, op2));
  transform(b, b + 5, b, bind1st(mu, 4));
  CGAL_assertion(equal(a, a + 5, b));

  return 0;
}
