#include <CGAL/MP_Float.h>
#include <CGAL/number_utils.h>
#include <cassert>
#include <cmath>
#include <limits>
#include <CGAL/Gmpq.h>

CGAL::Gmpq exact_value(const CGAL::MP_Float& value)
{
  CGAL::Gmpq result(0);
  for(auto exponent = value.max_exp(); exponent > value.min_exp(); )
    result = result * 65536 + value.of_exp(--exponent);
  const int exponent = static_cast<int>(value.min_exp() * 16);
  const CGAL::Gmpz scale = CGAL::Gmpz(1) << std::abs(exponent);
  return exponent >= 0 ? result * CGAL::Gmpq(scale) : result / CGAL::Gmpq(scale);
}

void check(double d)
{
  const CGAL::MP_Float value(d);
  assert(exact_value(value) == CGAL::Gmpq(d));
  assert(exact_value(value - CGAL::MP_Float(1)) == CGAL::Gmpq(d) - 1);
  assert(exact_value(-value) == -CGAL::Gmpq(d));
}

int main()
{
  for(double boundary : {CGAL::INTERN_MP_FLOAT::trunc_min,
                          CGAL::INTERN_MP_FLOAT::trunc_max})
    for(int exponent : {-64, -16, 0, 16, 64})
    {
      const double d = std::ldexp(boundary, exponent);
      check(d);
      check(std::nextafter(d, -std::numeric_limits<double>::infinity()));
      check(std::nextafter(d, std::numeric_limits<double>::infinity()));
    }

  for(double d : {0.0, 0.5, -0.5, 32767.5, -32768.5,
                   std::numeric_limits<double>::denorm_min(),
                   (std::numeric_limits<double>::min)(),
                   (std::numeric_limits<double>::max)()})
    check(d);

  const float f = static_cast<float>(CGAL::INTERN_MP_FLOAT::trunc_max);
  assert(CGAL::to_double(CGAL::MP_Float(f)) == f);
  const long double ld = CGAL::INTERN_MP_FLOAT::trunc_max;
  assert(CGAL::MP_Float(ld) == CGAL::MP_Float(static_cast<double>(ld)));
}
