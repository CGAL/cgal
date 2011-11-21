#include <iostream>
#include <CGAL/Overload.h>
#include <CGAL/compiler_config.h>

// This file requires a compiler with support for the C++0x features auto and lambda
int main()
{
#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
  auto overload = make_overload(
    function<int(int)>([](int) {  return 1; }),
    function<int(char)>([](char) {  return 2; }), 
    function<int(double)>([](double) { return 3; })
    );

  std::cout << o(1) << o('a') << " " << o(2.0) << std::endl;
#endif
  return 0;
}

