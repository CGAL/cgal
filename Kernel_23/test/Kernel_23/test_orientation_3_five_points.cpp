// Checks that the 5-point overload of Orientation_3, which returns
// { orientation(p,q,r,s), orientation(p,q,r,t) }, agrees with two calls
// to the 4-point overload, in particular in ranges where the static filter
// must not be used (underflow / overflow protection).

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

template <class K>
bool check(const typename K::Point_3& p, const typename K::Point_3& q,
           const typename K::Point_3& r, const typename K::Point_3& s,
           const typename K::Point_3& t)
{
  typename K::Orientation_3 orientation = K().orientation_3_object();
  const std::array<CGAL::Orientation, 2> o5 = orientation(p, q, r, s, t);
  const CGAL::Orientation os = orientation(p, q, r, s);
  const CGAL::Orientation ot = orientation(p, q, r, t);
  if(o5[0] != os || o5[1] != ot) {
    std::cerr << std::hexfloat << "Mismatch for\n p=" << p << "\n q=" << q << "\n r=" << r
              << "\n s=" << s << "\n t=" << t << "\n 5-point: " << int(o5[0]) << " " << int(o5[1])
              << "\n 4-point: " << int(os) << " " << int(ot) << std::defaultfloat << std::endl;
    return false;
  }
  return true;
}

template <class K>
int random_tests(int exp_xy, int exp_z, int n, std::mt19937& rng)
{
  typedef typename K::Point_3 P;
  std::uniform_int_distribution<int> coord(-7, 7);
  auto xy = [&]{ return std::ldexp(double(coord(rng)), exp_xy); };
  auto z  = [&]{ return std::ldexp(double(coord(rng)), exp_z); };

  int failures = 0;
  const P p(0, 0, 0);
  for(int i = 0; i < n; ++i) {
    const P q(xy(), xy(), z()), r(xy(), xy(), z()), s(xy(), xy(), z()), t(xy(), xy(), z());
    if(!check<K>(p, q, r, s, t))
      ++failures;
  }
  return failures;
}

template <class K>
int test_kernel()
{
  typedef typename K::Point_3 P;
  int failures = 0;

  // Simple configurations
  const P o(0,0,0), x(1,0,0), y(0,1,0), up(0,0,1), down(0,0,-1), in_plane(3,4,0);
  std::array<CGAL::Orientation, 2> res = K().orientation_3_object()(o, x, y, up, down);
  assert(res[0] == CGAL::POSITIVE && res[1] == CGAL::NEGATIVE);
  res = K().orientation_3_object()(o, x, y, in_plane, up);
  assert(res[0] == CGAL::COPLANAR && res[1] == CGAL::POSITIVE);

  // A configuration that used to be misclassified by the static filter: the x and y
  // extents are below the underflow threshold of the filter (1e-97), but not zero.
  failures += !check<K>(P(0, 0, 0),
                        P(-0x1.cp-536, -0x1p-538, -0x1p+1),
                        P( 0x1p-537,    0x0p+0,    0x1.4p+2),
                        P( 0x1.4p-536,  0x0p+0,   -0x1p+0),
                        P( 0x1.cp-536,  0x1p-538,  0x1.4p+2));

  std::mt19937 rng(42);
  failures += random_tests<K>(0, 0, 20000, rng);       // regular range
  failures += random_tests<K>(-538, 0, 50000, rng);    // tiny x/y extents: underflow guard
  failures += random_tests<K>(-536, -2, 20000, rng);
  failures += random_tests<K>(531, -340, 20000, rng);  // huge x/y extents: overflow guard
  failures += random_tests<K>(520, -320, 20000, rng);
  return failures;
}

int main()
{
  int failures = test_kernel<CGAL::Exact_predicates_inexact_constructions_kernel>();
  failures = test_kernel<CGAL::Exact_predicates_exact_constructions_kernel>();
  std::cout << failures << " failure(s)" << std::endl;
  assert(failures == 0);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

