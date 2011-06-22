#include <iostream>

#include <CGAL/Object.h>
#include <CGAL/assertions.h>

void test_object() {
  int i = 0;
  double j = 0.0;
  boost::variant<int, char, double> v = 23;
  CGAL::Object o = v;
  CGAL_assertion(CGAL::assign(i, o));
  CGAL_assertion(i == 23);
  //reassign the variant and assign it again
  v = 2.0;
  o = v;
  CGAL_assertion(!CGAL::assign(i, o));
  CGAL_assertion(CGAL::assign(j, o));
  CGAL_assertion(j == 2.0);
}

int main() {
  test_object();
}
