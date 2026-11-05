#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <cassert>
#include <limits>
#include <string>

using MPF = CGAL::MP_Float;
using QMPF = CGAL::Quotient<MPF>;

template<typename Integer>
void check(Integer value)
{
  const std::string digits = std::to_string(value);
  MPF expected(0);
  for(char digit : digits)
    if(digit != '-')
      expected = expected * MPF(10) + MPF(digit - '0');
  if(digits.front() == '-')
    expected = -expected;

  const MPF actual(value);
  assert(actual == expected);
  assert(QMPF(value) == QMPF(expected));
}

template<typename Integer>
void check_type()
{
  check(Integer(0));
  check(Integer(1));
  check((std::numeric_limits<Integer>::max)());
  check((std::numeric_limits<Integer>::min)());
}

int main()
{
  check_type<signed char>();
  check_type<unsigned char>();
  check_type<short>();
  check_type<unsigned short>();
  check_type<int>();
  check_type<unsigned int>();
  check_type<long>();
  check_type<unsigned long>();
  check_type<long long>();
  check_type<unsigned long long>();
  check_type<bool>();

  for(long long value : {32767LL, 32768LL, 65535LL, 65536LL, 9007199254740993LL})
  {
    check(value);
    check(-value);
  }
  assert(MPF(0.5) + MPF(0.25) == MPF(0.75));
}
