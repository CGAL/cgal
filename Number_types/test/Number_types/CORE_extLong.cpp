#include <CGAL/config.h>

#ifdef CGAL_USE_CORE
#include <CGAL/CORE/extLong.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <cassert>
#include <limits>
#include <vector>

void check_product(long a, long b)
{
  using boost::multiprecision::cpp_int;
  const cpp_int expected = cpp_int(a) * b;
  const CORE::extLong product = CORE::extLong(a) * CORE::extLong(b);
  const long limit = (std::numeric_limits<long>::max)();
  assert(!product.isNaN());
  if(expected >= limit || expected <= -cpp_int(limit))
  {
    assert(expected > 0 ? product.isInfty() : product.isTiny());
    assert(product.sign() == (expected > 0 ? 1 : -1));
  }
  else
  {
    assert(!product.isInfty() && !product.isTiny());
    assert(product.asLong() == expected.convert_to<long>());
  }
}
#endif

int main()
{
#ifdef CGAL_USE_CORE
  const long limit = (std::numeric_limits<long>::max)();
  std::vector<long> values = {0, 1, -1, 2, -2, 3, -3,
                            limit - 1, 1 - limit,
                            limit / 2, -limit / 2,
                            limit / 2 + 1, -limit / 2 - 1};
  for(long power = 1; power <= limit / 2; power *= 2)
  {
    values.push_back(power);
    values.push_back(-power);
    values.push_back(power - 1);
    values.push_back(1 - power);
  }
  for(long a : values)
    for(long b : values)
      check_product(a, b);

  CORE::extLong square(limit / 2 + 1);
  square *= square;
  assert(square.isInfty() && square.sign() > 0);
  assert((CORE::extLong::getPosInfty() * CORE::extLong(-2)).sign() < 0);
  assert((CORE::extLong::getNegInfty() * CORE::extLong(-2)).sign() > 0);
  assert((CORE::extLong(true) * CORE::extLong(2)).isNaN());
#endif
}
