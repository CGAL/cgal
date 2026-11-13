#include <CGAL/Interval_nt.h>
#include <CGAL/FPU.h>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cassert>
#include <cmath>
#include <limits>

template <bool Protected>
void check_ldexp()
{
  using Interval = CGAL::Interval_nt<Protected>;
  typename Interval::Protector protector;
  const double infinity = std::numeric_limits<double>::infinity();
  const double largest = std::numeric_limits<double>::max();
  const double smallest = std::numeric_limits<double>::denorm_min();

  const auto overflow = CGAL::ldexp(Interval(largest), 1);
  assert(overflow.inf() == largest);
  assert(overflow.sup() == infinity);
  const auto negative_overflow = CGAL::ldexp(Interval(-largest), 1);
  assert(negative_overflow.inf() == -infinity);
  assert(negative_overflow.sup() == -largest);

  const auto underflow = CGAL::ldexp(Interval(smallest), -1);
  assert(underflow.inf() <= 0);
  assert(underflow.sup() >= smallest);
  const auto negative_underflow = CGAL::ldexp(Interval(-smallest), -1);
  assert(negative_underflow.inf() <= -smallest);
  assert(negative_underflow.sup() >= 0);

  const auto scaled_up = CGAL::ldexp(Interval(smallest), 1074);
  assert(scaled_up.inf() <= 1 && scaled_up.sup() >= 1);
  assert(scaled_up.sup() < 2);
  const auto scaled_down = CGAL::ldexp(Interval(largest), -1075);
  assert(scaled_down.inf() > 0);
  assert(scaled_down.sup() < 1);

  for(const int exponent : {std::numeric_limits<int>::min(), -1075, -1074, -1022,
                            -1, 0, 1, 1023, 1024, 1074, std::numeric_limits<int>::max()})
  {
    const auto zero = CGAL::ldexp(Interval(0), exponent);
    assert(zero.inf() <= 0 && zero.sup() >= 0);
    const auto interval = CGAL::ldexp(Interval(-1, 1), exponent);
    assert(interval.inf() <= 0 && interval.sup() >= 0);
    assert(interval.inf() <= interval.sup());
  }

  const auto ordinary = CGAL::ldexp(Interval(-3, 5), 4);
  assert(ordinary.inf() <= -48 && ordinary.sup() >= 80);

  using Reference = boost::multiprecision::cpp_bin_float_100;
  for(const double value : {-largest, -1.0, -smallest, 0.0, smallest,
                            std::numeric_limits<double>::min(), 0.1, 1.0, largest})
    for(const int exponent : {-2200, -1075, -1074, -1023, -1022, -1, 0,
                              1, 1023, 1024, 1074, 2200})
    {
      const Reference expected = boost::multiprecision::ldexp(Reference(value), exponent);
      const auto result = CGAL::ldexp(Interval(value), exponent);
      assert(Reference(result.inf()) <= expected);
      assert(Reference(result.sup()) >= expected);
    }
}

int main()
{
  const auto rounding_mode = CGAL::FPU_get_cw();
  check_ldexp<true>();
  assert(CGAL::FPU_get_cw() == rounding_mode);
  check_ldexp<false>();
  assert(CGAL::FPU_get_cw() == rounding_mode);
}
