#include <iostream>
#include <CGAL/array.h>
#include <CGAL/use.h>

struct B {};

struct A // move-only class, move-constructible from B
{
  A(B&&) {}

  A(const A&) = delete;
  A(A&&) = default;
};

int main()
{
  // this test requires C++17 mandatory return-value optimization (RVO)
  std::array<A, 1> a = CGAL::make_array<A>(B());
  auto b = CGAL::make_array<double>(1u);
  static_assert(std::is_same_v<decltype(b), std::array<double, 1>>);
  CGAL_USE(a);
  CGAL_USE(b);
  return 0;
}
