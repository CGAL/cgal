#include <CGAL/tl/expected.hpp>
#include <iostream>

using CGAL::cpp23::expected;
using CGAL::cpp23::unexpected;

expected<int, const char*> divide(int a, int b)
{
  if (b == 0)
    return unexpected("division by zero");
  return a / b;
}

int main()
{
  auto r = divide(10, 2);
  if (r)
    std::cout << "Result: " << *r << std::endl;
  else
    std::cout << "Error: " << r.error() << std::endl;

  return 0;
}
