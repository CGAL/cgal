#include <CGAL/compiler_config.h>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#include <iostream>
#include <tuple>
#include <functional>
#include <CGAL/Overload.h>
#include <CGAL/assertions.h>
#include <boost/variant.hpp>


void test_overload() {
  using namespace CGAL;
  using std::function;
  
  auto t = std::make_tuple(function<int(int)>([](int) {  return 1; }),
                           function<int(char)>([](char) {  return 2; }), 
                           function<int(double)>([](double) { return 3; }));

  // manual construction
  Overload<function<int(int)>, function<int(char)>, function<int(double)> > o2(t);
  
  // basic overload checking
  auto o = make_overload(function<int(int)>([](int) {  return 1; }),
                         function<int(char)>([](char) {  return 2; }), 
                         function<int(double)>([](double) { return 3; }));

  CGAL_assertion(o(1) == 1);
  CGAL_assertion(o('a') == 2);
  CGAL_assertion(o(2.0) == 3);

  // check for variants
  boost::variant<int, char, double> v1 = 1;
  boost::variant<int, char, double> v2 = 'a';
  boost::variant<int, char, double> v3 = 2.0;

  CGAL_assertion(boost::apply_visitor(o, v1) == 1);
  CGAL_assertion(boost::apply_visitor(o, v2) == 2);
  CGAL_assertion(boost::apply_visitor(o, v3) == 3);
}

#endif

int main() {
#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)
  test_overload();
#endif
}
