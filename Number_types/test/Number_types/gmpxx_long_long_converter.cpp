#include <CGAL/config.h>

#ifdef CGAL_USE_GMPXX
#include <CGAL/gmpxx.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <limits>
#include <random>
#include <string>

template <class Integer>
void check(Integer value)
{
  const mpz_class expected(std::to_string(value));
  assert((CGAL::NT_converter<Integer, mpz_class>()(value) == expected));
  assert((CGAL::NT_converter<Integer, mpq_class>()(value) == mpq_class(expected)));
}
#endif

int main()
{
#ifdef CGAL_USE_GMPXX
  for(long long value : {0LL, 1LL, -1LL,
                         (1LL << 53) + 1, -(1LL << 53) - 1,
                         (std::numeric_limits<long long>::min)(),
                         (std::numeric_limits<long long>::max)()})
    check(value);
  for(unsigned long long value : {0ULL, 1ULL, (1ULL << 63) + 1,
                                  (std::numeric_limits<unsigned long long>::max)()})
    check(value);
  std::mt19937_64 random(5700);
  for(int i = 0; i < 1000; ++i)
  {
    const unsigned long long value = random();
    check(value);
    const long long signed_value = static_cast<long long>(value >> 1);
    check(signed_value);
    check(-signed_value);
  }

  using K = CGAL::Filtered_kernel<CGAL::Simple_cartesian<long long>>;
  const K::Point_2 p((1LL << 60) + 1, 0);
  const K::Point_2 q((1LL << 60) + 2, 0);
  assert(K().equal_2_object()(p, p));
  assert(!K().equal_2_object()(p, q));
  assert(K().orientation_2_object()(p, q, K::Point_2((1LL << 60) + 1, 1)) == CGAL::LEFT_TURN);
#endif
}
