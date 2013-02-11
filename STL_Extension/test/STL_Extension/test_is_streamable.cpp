#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/is_streamable.h>
#include <utility> // std::pair
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

struct A {};
struct B {};
struct C {};
struct D {};

using std::ostream;
using std::istream;

ostream& operator<<(ostream& os, const B&) { return os; }
istream& operator>>(istream& is, const B&) { return is; }

ostream& operator<<(ostream& os, const C&) { return os; }
istream& operator>>(istream& is, const D&) { return is; }

int main() {
  using CGAL::is_streamable;

  CGAL_static_assertion(!is_streamable<A>::value);
  CGAL_static_assertion(is_streamable<B>::value);
  CGAL_static_assertion(!is_streamable<C>::value);
  CGAL_static_assertion(!is_streamable<D>::value);
  CGAL_static_assertion(is_streamable<int>::value);
  CGAL_static_assertion(is_streamable<double>::value);
  CGAL_static_assertion(! (is_streamable<std::pair<int, int> >::value) );
  CGAL_static_assertion(  (is_streamable<boost::tuple<int, int> >::value) );
}
