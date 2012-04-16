#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/is_streamable.h>
#include <utility> // std::pair

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
  typedef std::vector<int>::const_iterator vector_it;
  typedef int* int_p;
  using CGAL::is_streamable;

  CGAL_static_assertion(!is_streamable<A>::value);
  CGAL_static_assertion(is_streamable<B>::value);
  CGAL_static_assertion(!is_streamable<C>::value);
  CGAL_static_assertion(!is_streamable<D>::value);
  CGAL_static_assertion(is_streamable<int>::value);
  CGAL_static_assertion(is_streamable<double>::value);
  CGAL_static_assertion(! (is_streamable<std::pair<int, int> >::value) );
}
