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
#if ! defined(_MSC_VER) || (_MSC_VER > 1916)
  // This test requires C++17 mandatory return-value optimization (RVO).
  //
  // MSVC-2017 does not implement C++17 correctly
  // See https://godbolt.org/z/7Y34Y1c53
  // and commit 15349f0bdafe60b85697f9d142c2652200d968e8
  // where we introduced a workaround for MSVC-2017: disable the correct
  // forwarding of arguments in the `make_array` function.
  // For that reason we disable this test for MSVC-2017)
  std::array<A, 1> a = CGAL::make_array<A>(B());
  CGAL_USE(a);
#endif
  auto b = CGAL::make_array<double>(1u);
  static_assert(std::is_same_v<decltype(b), std::array<double, 1>>);
  CGAL_USE(b);
  return 0;
}
