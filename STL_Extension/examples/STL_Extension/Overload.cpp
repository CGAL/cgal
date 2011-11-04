#include <iostream>
#include <CGAL/Overload.h>

// This file requires a compiler with support for the C++0x features auto and lambda
int main()
{
  auto overload = make_overload(
    function<int(int)>([](int) {  return 1; }),
    function<int(char)>([](char) {  return 2; }), 
    function<int(double)>([](double) { return 3; })
    );

  std::cout << o(1) << o('a') << " " << o(2.0) << std::endl;

  return 0;
}

