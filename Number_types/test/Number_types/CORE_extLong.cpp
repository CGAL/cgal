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
template <class Integer>
void check_unsigned_constructors()
{
  const long limit = (std::numeric_limits<long>::max)();
  auto check = [limit](Integer value)
  {
    const CORE::extLong number(value);
    if(boost::multiprecision::cpp_int(value) >= limit)
      assert(number.isInfty());
    else
      assert(!number.isInfty() && number.asLong() == static_cast<long>(value));
  };
  check(0);
  check((std::numeric_limits<Integer>::max)());
  for(int bit = 0; bit < std::numeric_limits<Integer>::digits; ++bit)
    check(Integer(1) << bit);
}

#endif

int main()
{
#ifdef CGAL_USE_CORE
  check_unsigned_constructors<unsigned int>();
  check_unsigned_constructors<unsigned long>();
  check_unsigned_constructors<std::size_t>();
  const long limit = (std::numeric_limits<long>::max)();
  std::vector<long> values = {0, 1, -1, 2, -2, 3, -3,
                            limit - 1, 1 - limit,
                            limit / 2, -limit / 2,
                            limit / 2 + 1, -limit / 2 - 1};
  if(limit % 7 == 0)
  {
    values.push_back(limit / 7);
    values.push_back(-limit / 7);
  }
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
  for(const CORE::extLong infinity : {CORE::extLong::getPosInfty(),
                                     CORE::extLong::getNegInfty()})
  {
    for(bool reverse : {false, true})
    {
      CORE::extLong product = reverse ? CORE::extLong(0) : infinity;
      product *= reverse ? infinity : CORE::extLong(0);
      assert(!product.isNaN() && !product.isInfty() && !product.isTiny());
      assert(product.asLong() == 0);
    }
    for(int sign : {-1, 1})
    {
      const CORE::extLong product = infinity * CORE::extLong(sign);
      assert(product.sign() == infinity.sign() * sign);
    }
  }
  assert((CORE::extLong(true) * CORE::extLong(0)).isNaN());
  assert((CORE::extLong(0) * CORE::extLong(true)).isNaN());
  CORE::extLong zero(0);
  zero *= zero;
  assert(!zero.isNaN() && zero.asLong() == 0);
#endif
}
