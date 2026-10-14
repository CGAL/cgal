#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <array>
#include <cassert>

using Nt_traits = CGAL::CORE_algebraic_number_traits;
using Rational = Nt_traits::Rational;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Nt_traits::Algebraic>;
using Traits = CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

void test_ellipse(const std::array<Rational, 6>& coefficients)
{
  Traits traits;
  auto construct = traits.construct_curve_2_object();
  const auto& c = coefficients;
  auto automatic = construct(c[0], c[1], c[2], c[3], c[4], c[5]);
  auto reversed = construct(-c[0], -c[1], -c[2], -c[3], -c[4], -c[5]);
  assert(automatic.is_valid());
  assert(automatic.is_full_conic());
  assert(automatic.orientation() == CGAL::CLOCKWISE);
  assert(reversed.orientation() == CGAL::COUNTERCLOCKWISE);

  for(int sign : {-1, 1})
  {
    for(auto orientation : {CGAL::CLOCKWISE, CGAL::COUNTERCLOCKWISE})
    {
      auto curve = construct(sign*c[0], sign*c[1], sign*c[2],
                             sign*c[3], sign*c[4], sign*c[5], orientation);
      assert(curve.is_valid());
      assert(curve.is_full_conic());
      assert(curve.orientation() == orientation);

      const auto& expected = (orientation == CGAL::CLOCKWISE) ? automatic : reversed;
      assert(curve.r() == expected.r());
      assert(curve.s() == expected.s());
      assert(curve.t() == expected.t());
      assert(curve.u() == expected.u());
      assert(curve.v() == expected.v());
      assert(curve.w() == expected.w());

      CGAL::Arrangement_2<Traits> arrangement(&traits);
      CGAL::insert(arrangement, curve);
      assert(arrangement.is_valid());
      assert(arrangement.number_of_vertices() == 2);
      assert(arrangement.number_of_edges() == 2);
      assert(arrangement.number_of_faces() == 2);
    }
  }
}

int main()
{
  test_ellipse({1, 1, 0, 0, 0, -1});
  test_ellipse({2, 3, 1, -5, 4, -10});
  test_ellipse({Rational(1, 2), Rational(1, 3), Rational(1, 4),
                Rational(-1, 5), Rational(1, 6), Rational(-1, 7)});
}
