
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using boost::tuple;
using boost::make_tuple;
using boost::tie;
using boost::get;

int main()
{
  tuple<int, double> a, b, c;
  a = tuple<int, double>();
  b = tuple<int, double>(1);
  c = tuple<int, double>(1, 3.14);
  a = make_tuple(1, 2.57);

  int i;
  double d;
  tie(i, d) = a;
  i = a.get<0>();
  d = a.get<1>();

  return 0;
}
