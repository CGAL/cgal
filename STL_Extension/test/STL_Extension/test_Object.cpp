#include <iostream>

#include <CGAL/Object.h>
#include <CGAL/assertions.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

void test_object() {
  int i = 0;
  double j = 0.0;
  boost::optional< boost::variant<int, char, double> > v(23);
  CGAL::Object o = v;
  CGAL_assertion(!o.empty());
  CGAL_assertion(CGAL::assign(i, o));
  CGAL_assertion(i == 23);
  //reassign the variant and assign it again
  v = 2.0;
  o = v;
  CGAL_assertion(!CGAL::assign(i, o));
  CGAL_assertion(CGAL::assign(j, o));
  CGAL_assertion(j == 2.0);
  //empty optional
  boost::optional< boost::variant<int, char, double> > v2;
  CGAL::Object o2 = v2;
  CGAL_assertion(o2.empty());
}

int main() {
  test_object();
}
