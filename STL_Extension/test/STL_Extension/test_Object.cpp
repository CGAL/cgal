#include <iostream>

#include <CGAL/Object.h>
#include <CGAL/assertions.h>

#include <variant>
#include <optional>
#include <cassert>

void from_opt_var() {
  int i = 0;
  double j = 0.0;

  std::optional< std::variant<int, char, double> > v(23);
  CGAL::Object o = v;
  assert(!o.empty());
  assert(CGAL::assign(i, o));
  assert(i == 23);
  //reassign the variant and assign it again
  v = 2.0;
  o = v;
  assert(!CGAL::assign(i, o));
  assert(CGAL::assign(j, o));
  assert(j == 2.0);
  //empty optional
  std::optional< std::variant<int, char, double> > v2;
  CGAL::Object o2 = v2;
  assert(o2.empty());
}

void from_var() {
  int i = 0;

  std::variant<int, char, double> v(23);
  CGAL::Object o = v;
  assert(!o.empty());
  assert(CGAL::assign(i, o));
  assert(i == 23);
}

struct Foo {
};

void make_object_and_assign() {
  int i = 23, j = 0;

  CGAL::Object o = CGAL::make_object(i);
  assert(CGAL::assign(j, o));
  assert(j == i);
  assert(CGAL::object_cast<Foo>(&o) == nullptr);
  assert(CGAL::object_cast<int>(&o) != nullptr);
}

void safe_bool() {
  CGAL::Object o;
  assert(!o);
  CGAL::Object o2 = CGAL::make_object(23);
  assert(o2);

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
