#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/is_streamable.h>
#include <utility> // std::pair

struct A {};
struct B {};

using std::ostream;
using std::istream;

ostream& operator<<(ostream& os, const B&) { return os; }
istream& operator>>(istream& is, const B&) { return is; }

int main() {
  typedef std::vector<int>::const_iterator vector_it;
  typedef int* int_p;
  using CGAL::is_streamable;

  CGAL_static_assertion(!is_streamable<A>::value);
  CGAL_static_assertion(is_streamable<B>::value);
  CGAL_static_assertion(is_streamable<int>::value);
  CGAL_static_assertion(is_streamable<double>::value);
  CGAL_static_assertion(! (is_streamable<std::pair<int, int> >::value) );
}
