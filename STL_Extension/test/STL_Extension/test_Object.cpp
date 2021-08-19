#include <iostream>

#include <CGAL/Object.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

void from_opt_var() {
  int i = 0;
  double j = 0.0;
  CGAL_USE(i);   CGAL_USE(j);
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

void from_var() {
  int i = 0;
  CGAL_USE(i);
  boost::variant<int, char, double> v(23);
  CGAL::Object o = v;
  CGAL_assertion(!o.empty());
  CGAL_assertion(CGAL::assign(i, o));
  CGAL_assertion(i == 23);
}

struct Foo {
};

void make_object_and_assign() {
  int i = 23, j = 0;
  CGAL_USE(j);
  CGAL::Object o = CGAL::make_object(i);
  CGAL_assertion(CGAL::assign(j, o));
  CGAL_assertion(j == i);
  CGAL_assertion(CGAL::object_cast<Foo>(&o) == nullptr);
  CGAL_assertion(CGAL::object_cast<int>(&o) != nullptr);
}

void safe_bool() {
  CGAL::Object o;
  CGAL_assertion(!o);
  CGAL::Object o2 = CGAL::make_object(23);
  CGAL_assertion(o2);

  // dummy code, we want to bork on this
  // if(o == o2) ;
  // if(o < 0) ;
}

void test_object() {
  make_object_and_assign();
  from_opt_var();
  from_var();
}

int main() {
  test_object();
}
